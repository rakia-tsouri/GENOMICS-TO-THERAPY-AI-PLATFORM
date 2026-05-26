from Bio import Entrez, SeqIO
from typing import Tuple, Optional
import io

# Default settings. In production, these should be from env vars.
Entrez.email = "ai-bio-project@example.com"  # NCBI requires an email

def fetch_sequence_by_id(gene_id: str) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """
    Fetches sequence and metadata from NCBI Entrez API.
    Returns (gene_id, sequence, error_message).
    """
    try:
        # Step 1: Search for the ID
        # Note: gene_id could be an accession like NM_000546
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="fasta", retmode="text")
        record_content = handle.read()
        handle.close()
        
        if not record_content:
            return None, None, f"No sequence found for ID: {gene_id}"
            
        # Parse the fasta string
        fasta_io = io.StringIO(record_content)
        record = SeqIO.read(fasta_io, "fasta")
        
        return record.id, str(record.seq), None
        
    except Exception as e:
        return None, None, f"NCBI Fetch Error: {str(e)}"
