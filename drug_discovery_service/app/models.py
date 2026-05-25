from pydantic import BaseModel
from typing import List, Dict, Optional

class Structure3D(BaseModel):
    pdb_file_path: str
    plddt_mean: float
    plddt_per_residue: List[float]

class BindingSite(BaseModel):
    positions: List[int]
    ligand: str

class Annotation(BaseModel):
    binding_sites: Optional[List[BindingSite]] = []
    uniprot_id: Optional[str] = None

class PredictRequest(BaseModel):
    gene_id: str
    protein_seq: str
    structure_3d: Structure3D
    annotation: Optional[Annotation] = None

class Pocket(BaseModel):
    residues: List[int]
    center: List[float]
    volume: float
    score: float

class DrugCandidate(BaseModel):
    chembl_id: Optional[str]
    smiles: str
    binding_score: float
    pchembl_value: Optional[float]
    lipinski_pass: bool
    admet: Dict[str, float]
    toxicity_risk: str

class PredictResponse(BaseModel):
    gene_id: str
    pocket: Pocket
    top_candidates: List[DrugCandidate]
    visualization_html_path: str
