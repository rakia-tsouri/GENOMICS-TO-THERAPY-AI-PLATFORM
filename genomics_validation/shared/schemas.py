from pydantic import BaseModel, Field
from typing import List, Optional

class ORFRecord(BaseModel):
    start: int
    end: int
    length: int
    dna_seq: str
    protein_seq: str
    protein_length: int

class ValidationRequest(BaseModel):
    dna_sequence: Optional[str] = None
    gene_id: Optional[str] = None

class ValidationResponse(BaseModel):
    valid: bool
    gene_id: Optional[str] = None
    dna_length: int = 0
    gc_percent: float = 0.0
    orfs: List[ORFRecord] = []
    warnings: List[str] = []
    errors: List[str] = []
