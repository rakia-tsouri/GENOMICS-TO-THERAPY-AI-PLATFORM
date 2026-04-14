import re
from Bio.Seq import Seq
from typing import List, Tuple, Dict
from shared.schemas import ORFRecord

def clean_sequence(raw_data: str) -> Tuple[str, str, List[str]]:
    """
    Parses FASTA/Raw string and returns (gene_id, cleaned_sequence, warnings).
    """
    lines = raw_data.strip().split('\n')
    gene_id = "Unknown"
    sequence_parts = []
    warnings = []

    if lines and lines[0].startswith('>'):
        gene_id = lines[0][1:].strip().split()[0]
        sequence_parts = lines[1:]
    else:
        sequence_parts = lines

    raw_sequence = "".join(sequence_parts)
    # Remove whitespace, digits, etc.
    cleaned = re.sub(r'[^a-zA-Z]', '', raw_sequence).upper()
    
    return gene_id, cleaned, warnings

def validate_dna(sequence: str) -> List[str]:
    """
    Validates the DNA sequence for length and nucleotide composition.
    Returns a list of errors.
    """
    errors = []
    
    # 1. Length validation (150 - 100k bp)
    if len(sequence) < 150:
        errors.append(f"Sequence too short: {len(sequence)} bp (min 150 bp)")
    elif len(sequence) > 100000:
        errors.append(f"Sequence too long: {len(sequence)} bp (max 100,000 bp)")
        
    # 2. Nucleotide composition
    if not re.match(r'^[ATCGN]*$', sequence):
        errors.append("Invalid nucleotides: Only A, T, C, G, N allowed.")
        
    return errors

def calculate_gc_content(sequence: str) -> float:
    """Calculates GC% of the sequence."""
    if not sequence:
        return 0.0
    g = sequence.count('G')
    c = sequence.count('C')
    return round((g + c) / len(sequence) * 100, 2)

def find_all_orfs(sequence: str) -> List[ORFRecord]:
    """
    Scans the sequence in all 3 forward frames for ORFs.
    An ORF starts with ATG and ends with a Stop codon (TAA, TAG, TGA).
    """
    orfs = []
    stop_codons = {"TAA", "TAG", "TGA"}
    
    seq_obj = Seq(sequence)
    seq_len = len(sequence)

    for frame in range(3):
        for i in range(frame, seq_len - 2, 3):
            codon = sequence[i:i+3]
            if codon == "ATG":
                # Found a potential start. Look for an in-frame stop.
                for j in range(i + 3, seq_len - 2, 3):
                    stop_codon = sequence[j:j+3]
                    if stop_codon in stop_codons:
                        orf_dna = sequence[i:j+3]
                        # Translate the ORF
                        protein = str(Seq(orf_dna).translate(table=1))
                        # Remove trailing stop '*' if present or handle it
                        # Seq.translate usually adds a * at the end if it's a stop
                        protein_seq = protein.rstrip('*')
                        
                        orfs.append(ORFRecord(
                            start=i,
                            end=j+3,
                            length=len(orf_dna),
                            dna_seq=orf_dna,
                            protein_seq=protein_seq,
                            protein_length=len(protein_seq)
                        ))
                        break # Stop at the first in-frame stop codon for this ATG

    # Sort ORFs by length (descending) as requested to "keep the most long" or just return all
    # The user recommendation was "Return MULTIPLE ORFs".
    orfs.sort(key=lambda x: x.length, reverse=True)
    return orfs
