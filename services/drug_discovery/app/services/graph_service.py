from ml.dataset import smiles_to_graph, pdb_to_graph
from torch_geometric.data import Data
from typing import List, Optional

class GraphService:
    @staticmethod
    def protein_to_graph(pdb_path: str, pocket_residues: List[int]) -> Data:
        """Step C: Convert protein structure to GNN graph."""
        return pdb_to_graph(pdb_path, pocket_residues)

    @staticmethod
    def smiles_to_graph(smiles: str) -> Data:
        """Step C: Convert molecule SMILES to GNN graph."""
        return smiles_to_graph(smiles)
