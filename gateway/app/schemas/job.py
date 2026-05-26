from datetime import datetime
from typing import Any, List, Optional

from pydantic import BaseModel, ConfigDict, model_validator


class JobCreate(BaseModel):
    """Create a pipeline run. Supply DNA input (sequence or gene id) and/or a WSI path.

    At least one modality must be provided:
      * Track A (genomics): dna_sequence OR gene_id
      * Track B (histopathology): wsi_image_path
    """
    name: str = "Untitled analysis"
    project_id: Optional[int] = None

    dna_sequence: Optional[str] = None
    gene_id: Optional[str] = None
    wsi_image_path: Optional[str] = None
    target_mutations: List[str] = ["TP53", "IDH1", "KRAS"]

    @model_validator(mode="after")
    def _at_least_one_input(self):
        if not (self.dna_sequence or self.gene_id or self.wsi_image_path):
            raise ValueError(
                "Provide at least one input: dna_sequence, gene_id, or wsi_image_path."
            )
        return self


class JobSummary(BaseModel):
    """Lightweight row for list views."""
    model_config = ConfigDict(from_attributes=True)

    id: int
    name: str
    status: str
    current_step: str
    project_id: Optional[int]
    gene_id: Optional[str]
    created_at: datetime
    completed_at: Optional[datetime]


class JobRead(JobSummary):
    """Full job with all stage results."""
    error: str
    dna_sequence: Optional[str]
    wsi_image_path: Optional[str]
    target_mutations: Optional[List[str]]
    validation_result: Optional[dict[str, Any]]
    analysis_result: Optional[dict[str, Any]]
    drug_result: Optional[dict[str, Any]]
    histopathology_result: Optional[dict[str, Any]]
    fusion_result: Optional[dict[str, Any]]
