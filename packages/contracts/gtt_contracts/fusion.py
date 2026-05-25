"""Contract for the cross-modal Fusion Layer (gateway-side).

The fusion layer reconciles mutation evidence from:
  * Track A (genomic): sequence-level signal derived from validation/annotation
  * Track B (visual):  morphology-based mutation prediction from histopathology

Concordant signals raise confidence; discordant signals are flagged for review.
"""
from __future__ import annotations

from pydantic import BaseModel
from typing import List, Optional


class GeneFusion(BaseModel):
    gene: str
    genomic_signal: Optional[bool] = None     # Track A evidence (None = not assessed)
    visual_signal: Optional[bool] = None       # Track B evidence (None = not assessed)
    visual_probability: Optional[float] = None
    agreement: str                             # "concordant" | "discordant" | "single_modality"
    combined_confidence: float                 # 0..1
    flag_for_review: bool = False
    note: Optional[str] = None


class FusionResult(BaseModel):
    genes: List[GeneFusion] = []
    overall_agreement: str = "single_modality"  # worst-case across genes
    cohen_kappa: Optional[float] = None          # cross-modal concordance (when computable)
    flagged_genes: List[str] = []
