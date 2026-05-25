import pytest
from httpx import AsyncClient
from app.main import app
from unittest.mock import patch
import respx
from httpx import Response

@pytest.mark.asyncio
async def test_analyze_endpoint_mocked():
    """
    Test the /analyze endpoint with mocked internal services.
    """
    input_payload = {
        "valid": True,
        "gene_id": "TP53",
        "dna_length": 1182,
        "protein_seq": "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYPQGLNGTVNLFRNLNKSSPQPKKKPLDGEYFTLKIRGRERFEMFRE",
        "protein_length": 393,
        "gc_percent": 51.2,
        "warnings": []
    }

    # Mocking BLAST, UniProt, FoldCheck and Structure services
    # We use patch for async functions inside main.py
    with patch("app.main.run_blast_search") as mock_blast, \
         patch("app.main.get_uniprot_annotation") as mock_uniprot, \
         patch("app.main.check_foldability") as mock_fold, \
         patch("app.main.get_3d_structure") as mock_structure:
        
        from app.models import BlastResult, ProteinAnnotation, FoldCheck, Structure3D
        
        mock_blast.return_value = BlastResult(
            top_hit_name="Cellular tumor antigen p53",
            top_hit_organism="Homo sapiens",
            identity_percent=99.7,
            e_value=0.0,
            coverage_percent=100.0,
            uniprot_id="P04637",
            protein_status="known"
        )
        mock_uniprot.return_value = ProteinAnnotation(function="Tumor suppressor")
        mock_fold.return_value = FoldCheck(foldable=True, method="api")
        mock_structure.return_value = Structure3D(source="AlphaFold DB", status="success")

        from httpx import ASGITransport
        async with AsyncClient(transport=ASGITransport(app=app), base_url="http://test") as ac:
            response = await ac.post("/api/v1/analyze", json=input_payload)
        
        assert response.status_code == 200
        data = response.json()
        assert data["gene_id"] == "TP53"
        assert data["blast"]["protein_status"] == "known"
        assert data["annotation"]["function"] == "Tumor suppressor"
