"""CNN patch feature extractor (ResNet50 backbone).

Wraps a torchvision ResNet50 used as a frozen feature extractor: it maps each
256x256 RGB patch to a 2048-d feature vector (global-average-pooled conv
features). These per-patch features are then aggregated by the MIL module.

Key conventions (mirroring drug_discovery/app/services/gnn_service.py):
  - CPU/GPU detection via torch.cuda.is_available().
  - pretrained=False (weights=None) so NO weights are downloaded at build/runtime.
    SCAFFOLD: in production load ImageNet (or, better, a pathology SSL backbone
    such as CTransPath / UNI) weights here.
  - torch is import-guarded so a missing install yields a clean error to the
    caller rather than crashing at module import.

The extractor keeps the final-conv feature map accessible (via a forward hook)
so the Grad-CAM service can reuse the SAME backbone for explainability.
"""
import logging
from typing import List

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Heavy optional deps: torch / torchvision. Import-guarded.
# ---------------------------------------------------------------------------
try:
    import torch
    import torch.nn as nn
    from torchvision.models import resnet50

    _TORCH_AVAILABLE = True
except Exception as e:  # pragma: no cover
    torch = None
    nn = None
    resnet50 = None
    _TORCH_AVAILABLE = False
    logger.warning("torch/torchvision not available (%s); CNN features disabled.", e)


# Deterministic seed so demos are stable across runs.
_SEED = 42

# ImageNet normalization constants (used even with random weights for sane scaling).
_MEAN = np.array([0.485, 0.456, 0.406], dtype=np.float32)
_STD = np.array([0.229, 0.224, 0.225], dtype=np.float32)

FEATURE_DIM = 2048  # ResNet50 penultimate feature width.


class CNNFeatureExtractor:
    """ResNet50-based per-patch feature extractor."""

    def __init__(self):
        if not _TORCH_AVAILABLE:
            self.available = False
            self.device = None
            self.model = None
            self.last_conv_features = None
            return

        self.available = True
        # Deterministic init for stable demo outputs.
        torch.manual_seed(_SEED)

        self.device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        logger.info("CNNFeatureExtractor using device: %s", self.device)

        # SCAFFOLD: weights=None -> random init, NO network download at build time.
        # Replace with a trained pathology backbone for real predictions.
        backbone = resnet50(weights=None)

        # Use everything up to (and including) the global average pool; drop the
        # 1000-class fc head so we emit a 2048-d feature vector per patch.
        self.model = nn.Sequential(*list(backbone.children())[:-1])
        self.model.to(self.device)
        self.model.eval()

        # Storage for the last conv feature map (consumed by Grad-CAM). The
        # final conv block ("layer4") output is captured via a forward hook.
        self.last_conv_features = None
        # children()[:-2] indexing: layer4 is the block right before avgpool+fc.
        # In nn.Sequential above, avgpool is the last module, layer4 the one before.
        self._register_hook(backbone)

    def _register_hook(self, backbone):
        """Register a forward hook on layer4 to cache conv features for Grad-CAM."""

        def _hook(_module, _inp, out):
            self.last_conv_features = out

        # layer4 is the deepest conv stage of ResNet50.
        backbone.layer4.register_forward_hook(_hook)

    def _to_tensor(self, patches: List[np.ndarray]):
        """Stack uint8 HWC RGB patches into a normalized NCHW float tensor."""
        arr = np.stack(patches).astype(np.float32) / 255.0  # N,H,W,C
        arr = (arr - _MEAN) / _STD
        arr = np.transpose(arr, (0, 3, 1, 2))               # N,C,H,W
        return torch.from_numpy(arr).to(self.device)

    def extract_features(self, patches: List[np.ndarray], batch_size: int = 16) -> np.ndarray:
        """Return an (N, FEATURE_DIM) numpy array of per-patch features.

        Raises RuntimeError if torch is unavailable so the caller can degrade
        gracefully. Returns an empty (0, FEATURE_DIM) array when there are no
        patches.
        """
        if not self.available:
            raise RuntimeError("torch/torchvision unavailable; cannot extract CNN features.")

        if not patches:
            return np.zeros((0, FEATURE_DIM), dtype=np.float32)

        feats: List[np.ndarray] = []
        with torch.no_grad():
            for i in range(0, len(patches), batch_size):
                batch = patches[i : i + batch_size]
                x = self._to_tensor(batch)
                out = self.model(x)               # N,2048,1,1
                out = out.flatten(1)              # N,2048
                feats.append(out.cpu().numpy())
        return np.concatenate(feats, axis=0).astype(np.float32)


# Module-level singleton (matches `gnn_service = GNNService()` convention).
cnn_feature_extractor = CNNFeatureExtractor()
