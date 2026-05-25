import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool

class DrugProteinGNN(nn.Module):
    def __init__(self, protein_feat_dim=25, drug_feat_dim=5, hidden=128):
        super().__init__()
        # Encodeur protéine
        self.prot_conv1 = GCNConv(protein_feat_dim, hidden)
        self.prot_conv2 = GCNConv(hidden, hidden)
        
        # Encodeur drogue
        self.drug_conv1 = GCNConv(drug_feat_dim, hidden)
        self.drug_conv2 = GCNConv(hidden, hidden)
        
        # Couche de fusion
        self.fusion = nn.Sequential(
            nn.Linear(hidden * 2, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, 64),
            nn.ReLU(),
            nn.Linear(64, 1)   # score d'affinité (pIC50)
        )

    def forward(self, prot_data, drug_data):
        # Encodage protéine
        xp = F.relu(self.prot_conv1(prot_data.x, prot_data.edge_index))
        xp = self.prot_conv2(xp, prot_data.edge_index)
        xp = global_mean_pool(xp, prot_data.batch)
        
        # Encodage drogue
        xd = F.relu(self.drug_conv1(drug_data.x, drug_data.edge_index))
        xd = self.drug_conv2(xd, drug_data.edge_index)
        xd = global_mean_pool(xd, drug_data.batch)
        
        # Fusion et prédiction
        combined = torch.cat([xp, xd], dim=1)
        return self.fusion(combined)
