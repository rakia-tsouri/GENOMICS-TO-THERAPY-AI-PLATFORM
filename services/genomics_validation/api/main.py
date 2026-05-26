from fastapi import FastAPI, HTTPException
from shared.schemas import ValidationRequest, ValidationResponse, ORFRecord
from .logic import clean_sequence, validate_dna, calculate_gc_content, find_all_orfs
from .ncbi import fetch_sequence_by_id
from .cache import cache
import logging

app = FastAPI(title="Genomics Validation Service", version="1.0.0")

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

@app.post("/api/v1/validate", response_model=ValidationResponse)
async def validate_genomics(request: ValidationRequest):
    warnings = []
    errors = []
    active_gene_id = request.gene_id
    raw_dna = request.dna_sequence

    # 1. Handle Gene ID fetch if provided
    if request.gene_id and not raw_dna:
        cached_dna = cache.get_sequence(request.gene_id)
        if cached_dna:
            logger.info(f"Cache hit for gene_id: {request.gene_id}")
            raw_dna = cached_dna
        else:
            logger.info(f"Cache miss for gene_id: {request.gene_id}. Fetching from NCBI...")
            fetched_id, fetched_dna, fetch_error = fetch_sequence_by_id(request.gene_id)
            if fetch_error:
                return ValidationResponse(valid=False, errors=[fetch_error])
            
            raw_dna = fetched_dna
            active_gene_id = fetched_id
            cache.set_sequence(request.gene_id, raw_dna)

    if not raw_dna:
        return ValidationResponse(valid=False, errors=["No DNA sequence or Gene ID provided."])

    # 2. Parsing & Cleaning
    gene_id_from_fasta, cleaned_dna, cleaning_warnings = clean_sequence(raw_dna)
    warnings.extend(cleaning_warnings)
    if not active_gene_id:
        active_gene_id = gene_id_from_fasta

    # 3. Biological Validation
    validation_errors = validate_dna(cleaned_dna)
    if validation_errors:
        return ValidationResponse(valid=False, errors=validation_errors, gene_id=active_gene_id)

    # 4. GC Content Assessment
    gc_percent = calculate_gc_content(cleaned_dna)
    if gc_percent < 35 or gc_percent > 65:
        warnings.append(f"Atypical GC content: {gc_percent}% (typical range 35-65%)")

    # 5. ORF Discovery (Multiple ORFs)
    orfs = find_all_orfs(cleaned_dna)
    
    if not orfs:
        errors.append("No valid ORF found (Start ATG + In-frame Stop).")
        return ValidationResponse(valid=False, errors=errors, gene_id=active_gene_id, dna_length=len(cleaned_dna))

    # 6. Success Response
    return ValidationResponse(
        valid=True,
        gene_id=active_gene_id,
        dna_length=len(cleaned_dna),
        gc_percent=gc_percent,
        orfs=orfs,
        warnings=warnings
    )

@app.get("/health")
async def health_check():
    return {"status": "healthy"}
