import pytest
from app.services.fold_check_service import calculate_hydrophobicity_fallback

def test_hydrophobicity_calculation():
    # Sequence with mixed residues
    sequence = "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDE"
    mean_score, regions = calculate_hydrophobicity_fallback(sequence)
    
    assert isinstance(mean_score, float)
    assert isinstance(regions, list)
    # Check if regions overlap or are out of bounds
    for r in regions:
        assert r[0] >= 1
        assert r[1] <= len(sequence)
        assert r[0] < r[1]

def test_hydrophobicity_length_residue():
    # Very short sequence
    sequence = "AAA"
    mean_score, regions = calculate_hydrophobicity_fallback(sequence, window_size=3)
    assert mean_score == 1.8
    assert regions == [] # Ala is hydrophobic, so no disordered regions expected (threshold -1.0)

def test_hydrophobicity_disordered():
    # Very hydrophilic sequence (Arg/Lys)
    sequence = "RRRRRRRRRRRR"
    mean_score, regions = calculate_hydrophobicity_fallback(sequence, window_size=5)
    assert mean_score < -4.0
    assert len(regions) > 0
