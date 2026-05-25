"""Contract for the Protein Analysis Service (Track A, step 2).

Endpoint: POST /api/v1/analyze
Input is one selected ORF from the validation step.
"""
from __future__ import annotations

from pydantic import BaseModel
from typing import List, Optional


class AnalyzeRequest(BaseModel):
    """Mirror of the validation output for a single chosen ORF."""
    valid: bool = True
    gene_id: str
    dna_length: int = 0
    protein_seq: str
    protein_length: int
    gc_percent: float = 0.0
    warnings: List[str] = []


class BlastHit(BaseModel):
    name: str
    organism: str
    identity_percent: float
    e_value: float
    coverage_percent: float


class BlastResult(BaseModel):
    top_hit_name: Optional[str] = None
    top_hit_organism: Optional[str] = None
    identity_percent: Optional[float] = None
    e_value: Optional[float] = None
    coverage_percent: Optional[float] = None
    uniprot_id: Optional[str] = None
    protein_status: str = "unknown"  # "known", "novel", "unknown"


class ActiveSite(BaseModel):
    position: int
    description: str


class BindingSite(BaseModel):
    positions: List[int]
    ligand: str


class ProteinAnnotation(BaseModel):
    uniprot_id: Optional[str] = None
    function: Optional[str] = None
    domains: List[str] = []
    active_sites: List[ActiveSite] = []
    binding_sites: List[BindingSite] = []
    diseases: List[str] = []


class FoldCheck(BaseModel):
    foldindex_score: Optional[float] = None
    foldable: bool
    disordered_percent: Optional[float] = None
    disordered_regions: List[List[int]] = []
    method: str
    warning: Optional[str] = None


class Structure3D(BaseModel):
    source: str  # "AlphaFold DB", "RCSB PDB", "ESMFold"
    status: str = "success"
    uniprot_id: Optional[str] = None
    pdb_file_path: Optional[str] = None
    plddt_mean: Optional[float] = None
    plddt_per_residue: List[float] = []
    confidence_level: Optional[str] = None  # "very high", "confident", "low", "failed"
    reason: Optional[str] = None
    fallback: bool = False


class AnalysisResponse(BaseModel):
    gene_id: str
    protein_seq: str
    protein_length: int
    blast: BlastResult
    annotation: Optional[ProteinAnnotation] = None
    fold_check: FoldCheck
    structure_3d: Structure3D
    processing_time_sec: float
