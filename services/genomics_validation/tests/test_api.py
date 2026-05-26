import pytest
from api.logic import clean_sequence, validate_dna, calculate_gc_content, find_all_orfs

def test_find_all_orfs_detailed():
    sequence = "ATGGCCAAATAGNNNATGAAATGA"
    results = find_all_orfs(sequence)
    
    print("\n--- Analyse des ORFs trouvés ---")
    for orf in results:
        print(f"Position: {orf.start}-{orf.end} | Longueur: {orf.length}bp")
        print(f"  ADN: {orf.dna_seq}")
        print(f"  PROT: {orf.protein_seq} (Taille: {orf.protein_length})")
        
    assert results[0].protein_seq == "MAK"  # ATGGCCAAATAG -> MAK
    assert results[1].protein_seq == "MK"   # ATGAAATGA -> MK