import time
import logging
from fastapi import FastAPI, HTTPException
from app.models import Service1Input, AnalysisResponse, BlastResult, FoldCheck, Structure3D
from app.services.blast_service import run_blast_search
from app.services.uniprot_service import get_uniprot_annotation
from app.services.fold_check_service import check_foldability
from app.services.structure_service import get_3d_structure

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

app = FastAPI(
    title="Protein Analysis Service",
    description="Microservice 2 for biomedical AI pipeline: BLAST, UniProt, FoldCheck, and 3D Structure.",
    version="1.0.0"
)

@app.post("/api/v1/analyze", response_model=AnalysisResponse)
async def analyze_protein(input_data: Service1Input):
    """
    Complete protein analysis pipeline.
    """
    start_time = time.time()
    logger.info(f"Starting analysis for Gene ID: {input_data.gene_id}")
    
    sequence = input_data.protein_seq
    
    # Step A: BLAST Search
    try:
        blast_result = await run_blast_search(sequence)
    except Exception as e:
        logger.error(f"Error in BLAST step: {e}")
        blast_result = BlastResult(protein_status="unknown")

    # Step B: UniProt Annotation
    annotation = None
    if blast_result.uniprot_id:
        try:
            annotation = await get_uniprot_annotation(blast_result.uniprot_id)
        except Exception as e:
            logger.error(f"Error in UniProt step: {e}")

    # Step C: Fold Check
    try:
        fold_check = await check_foldability(sequence)
    except Exception as e:
        logger.error(f"Error in Fold Check step: {e}")
        # Fallback to a safe fail state
        fold_check = FoldCheck(foldable=False, method="failed", warning=str(e))

    # Step D: 3D Structure
    try:
        structure_3d = await get_3d_structure(input_data.gene_id, sequence, blast_result)
    except Exception as e:
        logger.error(f"Error in Structure step: {e}")
        structure_3d = Structure3D(source="failed", status="failed", reason=str(e))

    processing_time = time.time() - start_time
    logger.info(f"Analysis completed in {processing_time:.2f} seconds.")

    return AnalysisResponse(
        gene_id=input_data.gene_id,
        protein_seq=sequence,
        protein_length=input_data.protein_length,
        blast=blast_result,
        annotation=annotation,
        fold_check=fold_check,
        structure_3d=structure_3d,
        processing_time_sec=round(processing_time, 2)
    )

@app.get("/health")
async def health_check():
    return {"status": "healthy"}

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8001)
