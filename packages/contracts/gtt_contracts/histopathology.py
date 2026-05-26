"""Contract for the Histopathology Service (Track B).

Endpoint: POST /api/v1/predict-mutation

Pipeline: WSI -> patch extraction (Otsu + stain norm) -> CNN features
(ResNet/EfficientNet) -> MIL aggregation -> per-gene mutation prediction
-> Grad-CAM explainability overlay.
"""
from __future__ import annotations

from pydantic import BaseModel
from typing import List, Optional


class HistopathologyRequest(BaseModel):
    wsi_image_path: str                       # path to uploaded WSI/tile image (shared volume)
    cancer_type: Optional[str] = None         # e.g. "LUAD"
    target_mutations: List[str] = ["TP53", "IDH1", "KRAS"]
    magnification: int = 20
    patch_size: int = 256


class MutationPrediction(BaseModel):
    gene: str
    mutated: bool
    probability: float                        # P(mutated), 0..1
    confidence: str                           # "high" | "medium" | "low"


class HistopathologyResponse(BaseModel):
    status: str = "success"                   # "success" | "failed"
    num_patches: int = 0                      # candidate patches found
    patches_kept: int = 0                     # after Otsu/background filtering
    mutations: List[MutationPrediction] = []
    gradcam_overlay_path: Optional[str] = None
    model: str = "ResNet50+ABMIL"
    processing_time_sec: float = 0.0
    warnings: List[str] = []
    reason: Optional[str] = None              # set when status == "failed"
