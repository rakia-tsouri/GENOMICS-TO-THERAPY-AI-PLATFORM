import os
import time
from Bio.Blast import NCBIWWW
from app.utils.cache import cache, TTL_24H
from app.utils.parsers import parse_blast_xml
from app.models import BlastResult
import logging

logger = logging.getLogger(__name__)

USE_BLAST = os.getenv("USE_BLAST", "true").lower() == "true"

async def run_blast_search(sequence: str) -> BlastResult:
    """
    Run a blastp search against SwissProt using NCBI's public API.
    Includes Redis caching and rate limit management.
    """
    if not USE_BLAST:
        return BlastResult(protein_status="unknown")

    # Check cache
    cache_key = f"blast:{sequence}"
    cached_result = cache.get(cache_key)
    if cached_result:
        logger.info("Retrieved BLAST results from cache.")
        return BlastResult(**cached_result)

    try:
        logger.info("Starting NCBI BLAST search (this may take 30-60s)...")
        # NCBIWWW uses public API. qblast(program, database, query)
        # SwissProt is better for dev as it's smaller/faster than 'nr'
        result_handle = NCBIWWW.qblast("blastp", "swissprot", sequence)
        xml_data = result_handle.read()
        result_handle.close()

        hits = parse_blast_xml(xml_data)
        
        if not hits:
            result = BlastResult(protein_status="novel")
        else:
            top_hit = hits[0]
            status = "known" if top_hit["identity_percent"] > 90 else "novel"
            if top_hit["identity_percent"] < 30:
                status = "novel"
            
            result = BlastResult(
                top_hit_name=top_hit["name"],
                top_hit_organism=top_hit["organism"],
                identity_percent=round(top_hit["identity_percent"], 2),
                e_value=top_hit["e_value"],
                coverage_percent=round(top_hit["coverage_percent"], 2),
                uniprot_id=top_hit["uniprot_id"],
                protein_status=status
            )

        # Save to cache
        cache.set(cache_key, result.model_dump(), ttl=TTL_24H)
        return result

    except Exception as e:
        logger.error(f"BLAST search failed: {e}")
        return BlastResult(protein_status="unknown")
