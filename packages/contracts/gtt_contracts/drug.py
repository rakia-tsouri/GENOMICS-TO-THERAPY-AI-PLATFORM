"""Contract for the Drug Discovery Service (shared final output).

Endpoint: POST /api/v1/predict
"""
from __future__ import annotations

from pydantic import BaseModel
from typing import Dict, List, Optional


class DrugStructure3D(BaseModel):
    pdb_file_path: str
    plddt_mean: float = 0.0
    plddt_per_residue: List[float] = []


class DrugBindingSite(BaseModel):
    positions: List[int]
    ligand: str


class DrugAnnotation(BaseModel):
    binding_sites: Optional[List[DrugBindingSite]] = []
    uniprot_id: Optional[str] = None


class PredictRequest(BaseModel):
    gene_id: str
    protein_seq: str
    structure_3d: DrugStructure3D
    annotation: Optional[DrugAnnotation] = None


class Pocket(BaseModel):
    residues: List[int]
    center: List[float]
    volume: float
    score: float


class DrugCandidate(BaseModel):
    chembl_id: Optional[str] = None
    smiles: str
    binding_score: float
    pchembl_value: Optional[float] = None
    lipinski_pass: bool
    admet: Dict[str, float] = {}
    toxicity_risk: str


class PredictResponse(BaseModel):
    gene_id: str
    pocket: Pocket
    top_candidates: List[DrugCandidate]
    visualization_html_path: str
