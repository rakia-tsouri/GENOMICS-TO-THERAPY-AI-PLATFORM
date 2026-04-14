import pytest
from api.logic import clean_sequence, validate_dna, calculate_gc_content, find_all_orfs

def test_clean_sequence_fasta():
    fasta = ">TP53\nATGGAGGAGCCG\nCAGTCAGATCC"
    gene_id, cleaned, warnings = clean_sequence(fasta)
    assert gene_id == "TP53"
    assert cleaned == "ATGGAGGAGCCGCAGTCAGATCC"

def test_validate_dna_valid():
    # Length 150+ sequence
    dna = "ATG" * 50 + "TAA"
    errors = validate_dna(dna)
    assert len(errors) == 0

def test_validate_dna_invalid_length():
    dna = "ATGC"
    errors = validate_dna(dna)
    assert any("too short" in e for e in errors)

def test_validate_dna_invalid_chars():
    dna = "ATGCNXYZ" * 20
    errors = validate_dna(dna)
    assert any("Invalid nucleotides" in e for e in errors)

def test_gc_content():
    dna = "GGCC" + "AA"
    # 4/6 = 66.67%
    assert calculate_gc_content(dna) == 66.67

def test_find_all_orfs():
    # Sequence with 2 ORFs
    # 1. ATGGCC TAA (6 bp)
    # 2. ATG TGA (6 bp) - in a different frame or later
    dna = "ATGGCC" + "TAA" + "NNN" + "ATGAAA" + "TGA"
    # Pad to 150 for validation if needed, but logic works regardless
    orfs = find_all_orfs(dna)
    assert len(orfs) >= 2
    assert orfs[0].dna_seq == "ATGGCC" + "TAA" or orfs[0].dna_seq == "ATGAAA" + "TGA"
