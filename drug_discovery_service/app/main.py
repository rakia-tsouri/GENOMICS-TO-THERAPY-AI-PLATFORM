from fastapi import FastAPI, HTTPException, BackgroundTasks
from .models import PredictRequest, PredictResponse, DrugCandidate, Pocket
from .services.pocket_service import PocketService
from .services.drug_service import DrugService
from .services.graph_service import GraphService
from .services.gnn_service import gnn_service
from .services.report_service import ReportService
from .utils.cache import cache
import os
import uuid

app = FastAPI(title="Drug Discovery Microservice", version="1.0.0")

@app.get("/")
async def root():
    return {"message": "Drug Discovery Microservice is running"}

@app.post("/api/v1/predict", response_model=PredictResponse)
async def predict_drugs(request: PredictRequest):
    # Check cache
    cache_key = f"predict:{request.gene_id}:{request.structure_3d.pdb_file_path}"
    cached_result = cache.get(cache_key)
    if cached_result:
        return PredictResponse(**cached_result)

    # 1. Pocket Detection (Step A)
    uniprot_sites = []
    if request.annotation and request.annotation.binding_sites:
        uniprot_sites = [p for site in request.annotation.binding_sites for p in site.positions]
    pocket_data = PocketService.detect_pockets(request.structure_3d.pdb_file_path, uniprot_sites)
    
    pocket = Pocket(
        residues=pocket_data.get("residues", []),
        center=pocket_data.get("center", [0.0, 0.0, 0.0]),
        volume=pocket_data.get("volume", 0.0),
        score=pocket_data.get("score", 0.0)
    )

    # 2. Fetch Drug Candidates (Step B)
    uniprot_id = request.annotation.uniprot_id if request.annotation else None
    molecules = await DrugService.get_candidates(uniprot_id)
    
    # 3. Protein Graph Conversion (Step C)
    prot_graph = GraphService.protein_to_graph(request.structure_3d.pdb_file_path, pocket.residues)
    if prot_graph is None:
        raise HTTPException(status_code=500, detail="Failed to generate protein graph")

    # 4. GNN Inference and Ranking (Step D & E)
    candidates = []
    for mol in molecules:
        drug_graph = GraphService.smiles_to_graph(mol["smiles"])
        if drug_graph is None: continue
        
        binding_score = gnn_service.predict_affinity(prot_graph, drug_graph)
        
        admet = ReportService.calculate_admet(mol["smiles"])
        toxicity = ReportService.check_toxicity(mol["smiles"])
        
        candidates.append(DrugCandidate(
            chembl_id=mol["chembl_id"],
            smiles=mol["smiles"],
            binding_score=binding_score,
            pchembl_value=mol.get("pchembl_value"),
            lipinski_pass=True, # Already filtered in DrugService
            admet=admet,
            toxicity_risk=toxicity
        ))

    # Rank top 10
    top_candidates = ReportService.rank_candidates(candidates)[:10]

    # 5. Visualization (Step E)
    viz_filename = f"viz_{uuid.uuid4()}.html"
    viz_path = os.path.join(os.getenv("PDB_VOLUME", "/data/structures"), "visualizations", viz_filename)
    os.makedirs(os.path.dirname(viz_path), exist_ok=True)
    
    if top_candidates:
        ReportService.generate_3d_viz(
            request.structure_3d.pdb_file_path, 
            top_candidates[0].smiles, 
            viz_path
        )

    response = PredictResponse(
        gene_id=request.gene_id,
        pocket=pocket,
        top_candidates=top_candidates,
        visualization_html_path=viz_path
    )

    # Cache result
    cache.set(cache_key, response.dict())

    return response
