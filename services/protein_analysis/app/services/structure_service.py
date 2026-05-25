import httpx
import os
import asyncio
from typing import Optional, Dict, Any
from app.models import Structure3D, BlastResult
from app.utils.parsers import extract_plddt_from_pdb
import logging

logger = logging.getLogger(__name__)

DATA_PATH = os.getenv("DATA_PATH", "/data/structures")
os.makedirs(DATA_PATH, exist_ok=True)

async def get_3d_structure(gene_id: str, sequence: str, blast_result: BlastResult) -> Structure3D:
    """
    Acquire 3D structure: Retrieval for known proteins, prediction for novel ones.
    """
    uniprot_id = blast_result.uniprot_id
    
    # CASE 1: Known Protein (Identity > 90%)
    if blast_result.protein_status == "known" and uniprot_id:
        # Try AlphaFold DB
        af_url = f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"
        try:
            async with httpx.AsyncClient(timeout=20.0) as client:
                response = await client.get(af_url)
                if response.status_code == 200:
                    data = response.json()
                    if data and len(data) > 0:
                        pdb_url = data[0].get("pdbUrl")
                        if pdb_url:
                            pdb_res = await client.get(pdb_url)
                            if pdb_res.status_code == 200:
                                pdb_content = pdb_res.text
                                return save_and_create_structure(gene_id, pdb_content, "AlphaFold DB", uniprot_id)
        except Exception as e:
            logger.warning(f"AlphaFold DB retrieval failed: {e}")

        # Fallback to RCSB PDB if not in AlphaFold
        # (Simplified: check RCSB search API could be complex, we'll try a common pattern for now)
        # In a real app, we'd use RCSB Search API.

    # CASE 2: Novel Protein (Identity < 30%) or Fallback for others
    # ESMFold Prediction
    status, pdb_content, reason = await predict_with_esmfold(sequence)
    if status == "success" and pdb_content:
        return save_and_create_structure(gene_id, pdb_content, "ESMFold", uniprot_id)
    
    return Structure3D(
        source="ESMFold",
        status="failed",
        reason=reason,
        fallback=True,
        confidence_level="failed"
    )

async def predict_with_esmfold(sequence: str, max_retries: int = 3) -> tuple:
    """
    Call Meta's ESMFold API with retries and timeout.
    """
    url = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    
    if len(sequence) > 400:
        logger.warning(f"Sequence length ({len(sequence)}) exceeds ESMFold ideal limit (400).")

    for attempt in range(max_retries):
        try:
            async with httpx.AsyncClient(timeout=60.0) as client:
                # Meta API accepts raw sequence in body
                response = await client.post(url, data=sequence)
                if response.status_code == 200:
                    return "success", response.text, None
                else:
                    logger.warning(f"ESMFold attempt {attempt+1} failed with status {response.status_code}")
        except asyncio.TimeoutError:
            logger.warning(f"ESMFold attempt {attempt+1} timed out.")
        except Exception as e:
            logger.warning(f"ESMFold attempt {attempt+1} error: {e}")
        
        if attempt < max_retries - 1:
            await asyncio.sleep(2 * (attempt + 1)) # Exponential backoff

    return "failed", None, "timeout or api error after 3 attempts"

def save_and_create_structure(gene_id: str, pdb_content: str, source: str, uniprot_id: Optional[str]) -> Structure3D:
    """Save PDB to shared volume and extract metadata."""
    file_path = os.path.join(DATA_PATH, f"{gene_id}.pdb")
    try:
        with open(file_path, "w") as f:
            f.write(pdb_content)
        
        plddt_info = extract_plddt_from_pdb(pdb_content, source=source)
        
        return Structure3D(
            source=source,
            status="success",
            uniprot_id=uniprot_id,
            pdb_file_path=file_path,
            plddt_mean=plddt_info["mean"],
            plddt_per_residue=plddt_info["per_residue"],
            confidence_level=plddt_info["confidence"]
        )
    except Exception as e:
        logger.error(f"Failed to save PDB file for {gene_id}: {e}")
        return Structure3D(source=source, status="failed", reason="file_save_error")
