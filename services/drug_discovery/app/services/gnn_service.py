import torch
import os
from ml.model import DrugProteinGNN
from torch_geometric.data import Batch

class GNNService:
    def __init__(self):
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.model_path = os.getenv("MODEL_PATH", "models_saved/gnn_v1.pt")
        self.model = DrugProteinGNN()
        
        if os.path.exists(self.model_path):
            try:
                self.model.load_state_dict(torch.load(self.model_path, map_location=self.device))
                print(f"Loaded GNN model from {self.model_path}")
            except Exception as e:
                print(f"Failed to load model: {e}. Using untrained weights.")
        else:
            print(f"Model file {self.model_path} not found. Using untrained weights.")
            
        self.model.to(self.device)
        self.model.eval()

    def predict_affinity(self, prot_graph, drug_graph) -> float:
        """Step D: Predict binding affinity pIC50."""
        # Wrap in Batch for single inference
        prot_batch = Batch.from_data_list([prot_graph]).to(self.device)
        drug_batch = Batch.from_data_list([drug_graph]).to(self.device)
        
        with torch.no_grad():
            score = self.model(prot_batch, drug_batch)
            
        return float(score.item())

gnn_service = GNNService()
