"""WSI patch extraction service.

Pipeline responsibilities:
  1. Open a Whole Slide Image (WSI). Uses ``openslide`` when available for real
     pyramidal formats (.svs/.ndpi/.tiff); FALLS BACK to PIL/numpy for plain
     raster images (.png/.jpg) so the service runs in a demo without a real slide.
  2. Tile the image into ``patch_size`` x ``patch_size`` patches (default 256).
  3. Otsu thresholding + tissue-ratio filter to DROP mostly-background tiles
     (glass / white space). This part is implemented for real.
  4. Macenko stain normalization -- a documented SIMPLIFIED version here
     (see ``_macenko_normalize_simplified``).

Heavy / optional libraries (openslide) are import-guarded so a missing lib
produces a clean warning rather than crashing at module import.
"""
import logging
from typing import List, Dict, Any

import numpy as np
from PIL import Image

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional dependency: openslide (real WSI reader). Import-guarded.
# ---------------------------------------------------------------------------
try:
    import openslide  # type: ignore

    _OPENSLIDE_AVAILABLE = True
except Exception as e:  # pragma: no cover - depends on host libs
    openslide = None
    _OPENSLIDE_AVAILABLE = False
    logger.warning(
        "openslide not available (%s); falling back to PIL for plain raster images.", e
    )

# ---------------------------------------------------------------------------
# Optional dependency: scikit-image for a robust Otsu implementation.
# We keep a pure-numpy Otsu fallback so the tissue filter always works.
# ---------------------------------------------------------------------------
try:
    from skimage.filters import threshold_otsu  # type: ignore

    _SKIMAGE_AVAILABLE = True
except Exception as e:  # pragma: no cover
    threshold_otsu = None
    _SKIMAGE_AVAILABLE = False
    logger.warning("scikit-image not available (%s); using numpy Otsu fallback.", e)


# Minimum fraction of a patch that must be tissue (non-background) to keep it.
TISSUE_RATIO_THRESHOLD = 0.25


def _otsu_threshold(gray: np.ndarray) -> float:
    """Compute an Otsu threshold for a grayscale array (uint8, 0..255).

    Prefers scikit-image; otherwise uses a pure-numpy histogram implementation
    of Otsu's between-class variance maximization.
    """
    if _SKIMAGE_AVAILABLE:
        try:
            return float(threshold_otsu(gray))
        except Exception as e:  # e.g. constant image -> skimage raises
            logger.warning("skimage Otsu failed (%s); using numpy fallback.", e)

    # --- Pure numpy Otsu (real implementation) ---
    hist, _ = np.histogram(gray.ravel(), bins=256, range=(0, 256))
    hist = hist.astype(np.float64)
    total = hist.sum()
    if total == 0:
        return 128.0
    prob = hist / total
    omega = np.cumsum(prob)                     # class 0 weight (cumulative)
    mu = np.cumsum(prob * np.arange(256))       # class 0 cumulative mean
    mu_t = mu[-1]                               # global mean
    # Between-class variance for each threshold; guard div-by-zero.
    denom = omega * (1.0 - omega)
    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_b2 = np.where(denom > 0, (mu_t * omega - mu) ** 2 / denom, 0.0)
    return float(np.argmax(sigma_b2))


def _tissue_ratio(patch_rgb: np.ndarray, otsu_thr: float) -> float:
    """Fraction of pixels in the patch that are tissue (darker than Otsu thr).

    Background in H&E WSIs is bright (near white); tissue is darker/more saturated.
    """
    gray = patch_rgb.mean(axis=2)
    tissue_mask = gray < otsu_thr
    return float(tissue_mask.mean())


def _macenko_normalize_simplified(patch_rgb: np.ndarray) -> np.ndarray:
    """SIMPLIFIED Macenko stain normalization.

    SCAFFOLD: The full Macenko (2009) method converts RGB -> optical density,
    estimates the two principal stain vectors (H&E) via SVD of the OD matrix,
    projects pixels onto them, and rescales to a fixed reference stain matrix.
    That requires per-slide eigen-decomposition and reference constants.

    Here we apply a lightweight, deterministic stand-in that performs the
    *spirit* of normalization -- per-channel mean/contrast alignment to a fixed
    target -- so downstream features are less sensitive to global staining/
    illumination differences. Replace with a full Macenko implementation
    (e.g. ``staintools`` / ``torchstain``) for production use.
    """
    img = patch_rgb.astype(np.float32)
    # Target per-channel mean (typical normalized H&E look) and std.
    target_mean = np.array([180.0, 140.0, 190.0], dtype=np.float32)
    target_std = np.array([40.0, 40.0, 40.0], dtype=np.float32)
    out = np.empty_like(img)
    for c in range(3):
        ch = img[..., c]
        std = ch.std()
        if std < 1e-3:
            out[..., c] = ch  # near-constant channel: leave as-is
        else:
            out[..., c] = (ch - ch.mean()) / std * target_std[c] + target_mean[c]
    return np.clip(out, 0, 255).astype(np.uint8)


def _load_image_array(wsi_image_path: str, magnification: int) -> np.ndarray:
    """Load the WSI/tile as an RGB numpy array.

    Uses openslide for real pyramidal slides when available; otherwise PIL.
    The ``magnification`` argument is accepted for API parity; with openslide it
    could be used to pick the matching pyramid level. For the simplified reader
    we read the full (or thumbnail) image and tile it directly.
    """
    if _OPENSLIDE_AVAILABLE and wsi_image_path.lower().endswith(
        (".svs", ".ndpi", ".tif", ".tiff", ".mrxs", ".vms", ".scn")
    ):
        # SCAFFOLD: a production reader would select the pyramid level matching
        # `magnification` and use read_region over a tissue bounding box rather
        # than loading a single downsampled thumbnail.
        slide = openslide.OpenSlide(wsi_image_path)
        try:
            level = slide.level_count - 1  # lowest-res level for a safe demo read
            dims = slide.level_dimensions[level]
            region = slide.read_region((0, 0), level, dims).convert("RGB")
            return np.asarray(region)
        finally:
            slide.close()

    # Fallback: plain raster image via PIL.
    img = Image.open(wsi_image_path).convert("RGB")
    return np.asarray(img)


def extract_patches(
    wsi_image_path: str,
    patch_size: int = 256,
    magnification: int = 20,
) -> Dict[str, Any]:
    """Extract tissue patches from a WSI.

    Returns a dict:
        {
          "patches": List[np.ndarray]  # kept, stain-normalized RGB patches (uint8)
          "num_patches": int           # candidate tiles found (grid count)
          "patches_kept": int          # tiles passing the Otsu tissue filter
        }

    Raises on unrecoverable I/O errors so the caller (main.py) can degrade the
    whole request gracefully.
    """
    image = _load_image_array(wsi_image_path, magnification)
    h, w = image.shape[:2]

    # Global Otsu threshold computed once on the whole image for stable tiling.
    gray_full = image.mean(axis=2).astype(np.uint8)
    otsu_thr = _otsu_threshold(gray_full)
    logger.info("Otsu threshold for %s = %.1f", wsi_image_path, otsu_thr)

    kept_patches: List[np.ndarray] = []
    num_candidates = 0

    for y in range(0, h - patch_size + 1, patch_size):
        for x in range(0, w - patch_size + 1, patch_size):
            num_candidates += 1
            patch = image[y : y + patch_size, x : x + patch_size, :]

            # Otsu + tissue-ratio filter: drop mostly-background tiles.
            if _tissue_ratio(patch, otsu_thr) < TISSUE_RATIO_THRESHOLD:
                continue

            # Macenko (simplified) stain normalization on kept tiles.
            norm = _macenko_normalize_simplified(patch)
            kept_patches.append(norm)

    logger.info(
        "Extracted %d candidate tiles, kept %d after tissue filtering.",
        num_candidates,
        len(kept_patches),
    )
    return {
        "patches": kept_patches,
        "num_patches": num_candidates,
        "patches_kept": len(kept_patches),
    }
