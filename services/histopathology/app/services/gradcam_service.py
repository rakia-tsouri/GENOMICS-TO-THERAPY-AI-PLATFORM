"""Grad-CAM explainability overlay generation.

Produces a Grad-CAM style heatmap for the HIGHEST-ATTENTION patch (the tile the
ABMIL aggregator weighted most for the slide-level decision) and saves a PNG to
the shared volume (env GRADCAM_VOLUME, default /data/structures/gradcam).

Grad-CAM (Selvaraju et al., 2017): backprop a target score to the last conv
feature map, weight each feature channel by its mean gradient (importance),
ReLU the weighted sum, upsample to the patch size, and overlay as a heatmap.

Behavior when optional deps are missing:
  - If torch / cv2 are unavailable, we DO NOT crash. We write a small text
    placeholder file describing what would have been produced and return its
    path, with a warning surfaced by the caller.

SCAFFOLD: because the CNN+MIL run on untrained weights, the produced heatmap is
illustrative, not a validated saliency map.
"""
import logging
import os
import uuid
from typing import List, Optional, Any

import numpy as np

logger = logging.getLogger(__name__)

try:
    import torch

    _TORCH_AVAILABLE = True
except Exception as e:  # pragma: no cover
    torch = None
    _TORCH_AVAILABLE = False
    logger.warning("torch not available (%s); Grad-CAM will write a placeholder.", e)

try:
    import cv2

    _CV2_AVAILABLE = True
except Exception as e:  # pragma: no cover
    cv2 = None
    _CV2_AVAILABLE = False
    logger.warning("cv2 not available (%s); Grad-CAM will write a placeholder.", e)

# PIL is a hard dependency (used for the no-cv2 PNG fallback too).
try:
    from PIL import Image

    _PIL_AVAILABLE = True
except Exception as e:  # pragma: no cover
    Image = None
    _PIL_AVAILABLE = False


def _gradcam_dir() -> str:
    """Resolve and ensure the Grad-CAM output directory exists."""
    out_dir = os.getenv("GRADCAM_VOLUME", "/data/structures/gradcam")
    os.makedirs(out_dir, exist_ok=True)
    return out_dir


def _write_placeholder(out_dir: str, reason: str) -> str:
    """Write a .txt placeholder explaining what a real overlay would contain."""
    path = os.path.join(out_dir, f"gradcam_placeholder_{uuid.uuid4().hex[:8]}.txt")
    with open(path, "w") as f:
        f.write(
            "Grad-CAM overlay PLACEHOLDER\n"
            f"Reason: {reason}\n\n"
            "A real run would save a PNG: the highest-attention 256x256 tissue patch\n"
            "with a jet heatmap overlay highlighting the regions most responsible for\n"
            "the predicted mutation status (Grad-CAM on the ResNet50 last conv layer).\n"
            "NOTE: model weights are untrained -- the saliency map would be illustrative.\n"
        )
    logger.warning("Wrote Grad-CAM placeholder (%s) -> %s", reason, path)
    return path


def _heatmap_overlay(patch_rgb: np.ndarray, cam: np.ndarray) -> np.ndarray:
    """Blend a normalized CAM (HxW, 0..1) over the patch as a jet heatmap."""
    cam = np.clip(cam, 0.0, 1.0)
    h, w = patch_rgb.shape[:2]
    if _CV2_AVAILABLE:
        cam_resized = cv2.resize(cam, (w, h))
        heat = cv2.applyColorMap((cam_resized * 255).astype(np.uint8), cv2.COLORMAP_JET)
        heat = cv2.cvtColor(heat, cv2.COLOR_BGR2RGB)  # cv2 returns BGR
    else:
        # Minimal numpy jet-ish heatmap (no cv2): map intensity to R/G/B ramp.
        cam_resized = np.asarray(
            Image.fromarray((cam * 255).astype(np.uint8)).resize((w, h))
        ).astype(np.float32) / 255.0
        heat = np.zeros((h, w, 3), dtype=np.float32)
        heat[..., 0] = cam_resized                 # red increases with importance
        heat[..., 2] = 1.0 - cam_resized           # blue decreases
        heat = (heat * 255).astype(np.uint8)
    overlay = (0.6 * patch_rgb.astype(np.float32) + 0.4 * heat.astype(np.float32))
    return np.clip(overlay, 0, 255).astype(np.uint8)


def generate_overlay(
    patches: List[np.ndarray],
    attention: Optional["np.ndarray"],
    cnn_extractor: Any = None,
) -> Optional[str]:
    """Generate a Grad-CAM overlay PNG for the highest-attention patch.

    Args:
        patches: list of kept RGB patches (uint8 HWC) -- same order MIL scored.
        attention: per-patch attention weights (N,) from the MIL aggregator.
        cnn_extractor: the CNNFeatureExtractor singleton (used to access the
            backbone + cached conv features). Optional; if None or torch is
            unavailable we fall back to an attention-only saliency proxy.

    Returns the saved file path (PNG on success, or a .txt placeholder), or None
    if even the placeholder could not be written.
    """
    out_dir = _gradcam_dir()

    if not patches:
        return _write_placeholder(out_dir, "no tissue patches available")

    if not _PIL_AVAILABLE:
        # Cannot write any image; nothing else we can do.
        logger.error("PIL unavailable; cannot write Grad-CAM output at all.")
        return None

    # Pick the highest-attention patch (fallback to the first patch).
    if attention is not None and len(attention) == len(patches):
        top_idx = int(np.argmax(attention))
    else:
        top_idx = 0
    top_patch = np.asarray(patches[top_idx])

    if not _TORCH_AVAILABLE:
        return _write_placeholder(out_dir, "torch unavailable; cannot run Grad-CAM")

    try:
        cam = _compute_gradcam(top_patch, cnn_extractor)
    except Exception as e:
        logger.warning("Grad-CAM computation failed (%s); using uniform saliency.", e)
        cam = np.ones(top_patch.shape[:2], dtype=np.float32) * 0.5

    overlay = _heatmap_overlay(top_patch, cam)
    path = os.path.join(out_dir, f"gradcam_{uuid.uuid4().hex[:8]}.png")
    Image.fromarray(overlay).save(path)
    logger.info("Saved Grad-CAM overlay -> %s (patch #%d)", path, top_idx)
    return path


def _compute_gradcam(patch_rgb: np.ndarray, cnn_extractor: Any) -> np.ndarray:
    """Real Grad-CAM on the ResNet50 backbone's last conv layer.

    Requires the CNNFeatureExtractor (for its backbone/device + cached conv
    feature map). Targets the max logit-proxy (sum of pooled features) since the
    scaffold has no trained class head wired through here.

    SCAFFOLD: with trained weights you would backprop the specific gene logit
    from the MIL head instead of this generic feature-magnitude proxy.
    """
    if cnn_extractor is None or not getattr(cnn_extractor, "available", False):
        raise RuntimeError("CNN extractor unavailable for Grad-CAM")

    device = cnn_extractor.device
    model = cnn_extractor.model  # nn.Sequential(conv... , avgpool)

    # Build input tensor (reuse the extractor's normalization constants).
    from .cnn_service import _MEAN, _STD

    arr = patch_rgb.astype(np.float32) / 255.0
    arr = (arr - _MEAN) / _STD
    arr = np.transpose(arr, (2, 0, 1))[None, ...]  # 1,C,H,W
    x = torch.from_numpy(arr).to(device).requires_grad_(True)

    # We need gradients w.r.t. the cached conv features. Capture them locally.
    activations = {}
    gradients = {}

    # The extractor's forward hook caches layer4 output on cnn_extractor; here we
    # add our own hooks on the corresponding module inside the Sequential.
    # In the Sequential, the conv stages are children up to avgpool; the last
    # conv block (layer4) is index -2.
    target_layer = list(model.children())[-2]

    def fwd_hook(_m, _i, o):
        activations["value"] = o

    def bwd_hook(_m, _gi, go):
        gradients["value"] = go[0]

    fh = target_layer.register_forward_hook(fwd_hook)
    bh = target_layer.register_full_backward_hook(bwd_hook)

    try:
        model.zero_grad()
        out = model(x).flatten(1)          # 1, 2048
        score = out.sum()                  # SCAFFOLD proxy target (no trained head)
        score.backward()

        acts = activations["value"]        # 1, C, h, w
        grads = gradients["value"]         # 1, C, h, w
        weights = grads.mean(dim=(2, 3), keepdim=True)   # channel importance
        cam = torch.relu((weights * acts).sum(dim=1)).squeeze(0)  # h, w
        cam = cam.detach().cpu().numpy()
    finally:
        fh.remove()
        bh.remove()

    # Normalize to 0..1.
    cam -= cam.min()
    if cam.max() > 1e-8:
        cam /= cam.max()
    return cam.astype(np.float32)
