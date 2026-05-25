import os
import torch
import pandas as pd
from torch_geometric.data import Dataset, Data
from Bio.PDB import PDBParser
from rdkit import Chem
import numpy as np
from tqdm import tqdm

class ProteinLigandDataset(Dataset):
    def __init__(self, root, csv_file, protein_dir, transform=None, pre_transform=None, pre_filter=None):
        self.csv_file = csv_file
        self.protein_dir = protein_dir
        self.df = pd.read_csv(csv_file)
        super().__init__(root, transform, pre_transform, pre_filter)
        
        # Filtrer pour ne garder que les fichiers qui existent vraiment
        self.valid_indices = []
        # On définit processed_dir ici car super().__init__ le crée
        p_dir = os.path.join(root, 'processed')
        for i in range(len(self.df)):
            if os.path.exists(os.path.join(p_dir, f'data_{i}.pt')):
                self.valid_indices.append(i)
        print(f"Dataset filtered: {len(self.valid_indices)} valid samples out of {len(self.df)}")

    @property
    def raw_file_names(self):
        return [os.path.basename(self.csv_file)]

    @property
    def processed_file_names(self):
        # We will save processed data as data_0.pt, data_1.pt, etc.
        return [f'data_{i}.pt' for i in range(len(self.df))]

    def download(self):
        # Data is already local
        pass

    def process(self):
        protein_cache = {} # Cache pour ne pas recalculer la même protéine 100 fois
        
        for idx, row in tqdm(self.df.iterrows(), total=len(self.df), desc="Processing graphs"):
            processed_path = os.path.join(self.processed_dir, f'data_{idx}.pt')
            if os.path.exists(processed_path): continue # Skip si déjà fait
            
            pdb_id = row['pdb_id']
            smiles = row['smiles']
            affinity = row['affinity']
            
            # 1. Convert Drug SMILES to Graph
            drug_data = smiles_to_graph(smiles)
            
            # 2. Convert Protein PDB to Graph (avec Cache)
            if pdb_id in protein_cache:
                prot_data = protein_cache[pdb_id]
            else:
                pdb_path = os.path.join(self.protein_dir, pdb_id, f"{pdb_id}_protein.pdb")
                prot_data = pdb_to_graph(pdb_path)
                protein_cache[pdb_id] = prot_data
            
            if drug_data and prot_data:
                data = {
                    'prot_x': prot_data.x,
                    'prot_edge_index': prot_data.edge_index,
                    'drug_x': drug_data.x,
                    'drug_edge_index': drug_data.edge_index,
                    'y': torch.tensor([affinity], dtype=torch.float)
                }
                torch.save(data, processed_path)

    def len(self):
        return len(self.valid_indices)

    def get(self, idx):
        real_idx = self.valid_indices[idx]
        data = torch.load(os.path.join(self.processed_dir, f'data_{real_idx}.pt'))
        
        # Create separate Data objects for protein and drug
        prot_data = Data(x=data['prot_x'], edge_index=data['prot_edge_index'])
        drug_data = Data(x=data['drug_x'], edge_index=data['drug_edge_index'])
        y = data['y']
        
        return prot_data, drug_data, y

def smiles_to_graph(smiles: str) -> Data:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return None
    
    node_features = []
    for atom in mol.GetAtoms():
        features = [
            atom.GetAtomicNum(),
            atom.GetDegree(),
            atom.GetFormalCharge(),
            int(atom.GetIsAromatic()),
            atom.GetTotalNumHs(),
        ]
        node_features.append(features)
        
    edge_index = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        edge_index += [[i,j],[j,i]]
        
    if not edge_index: # Case of single atom molecules
        edge_index = [[0, 0]]
        
    return Data(
        x=torch.tensor(node_features, dtype=torch.float),
        edge_index=torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    )

def pdb_to_graph(pdb_path: str) -> Data:
    if not os.path.exists(pdb_path):
        return None
        
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("protein", pdb_path)
    except:
        return None

    nodes, coords = [], []
    for residue in structure.get_residues():
        # Clean residues (ignore water etc)
        if residue.get_id()[0] != ' ': continue
        
        aa = residue.get_resname()
        if 'CA' not in residue: continue
        ca = residue['CA'].get_vector()
        coords.append([ca[0], ca[1], ca[2]])
        nodes.append(aa_to_features(aa))

    if not coords: return None

    # Edges: KNN or distance based
    edge_index = []
    coords_np = np.array(coords)
    # Simple distance based edges (8A)
    for i in range(len(coords)):
        for j in range(i+1, len(coords)):
            dist = np.linalg.norm(coords_np[i] - coords_np[j])
            if dist < 8.0:
                edge_index += [[i,j],[j,i]]

    if not edge_index:
        edge_index = [[0, 0]]

    return Data(
        x=torch.tensor(nodes, dtype=torch.float),
        edge_index=torch.tensor(edge_index, dtype=torch.long).t().contiguous()
    )

def aa_to_features(aa):
    amino_acids = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 
                   'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    features = [0.0] * 25
    if aa in amino_acids:
        features[amino_acids.index(aa)] = 1.0
    
    hydrophobic = ['ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
    charged = ['ASP', 'GLU', 'LYS', 'ARG', 'HIS']
    
    if aa in hydrophobic: features[20] = 1.0
    if aa in charged: features[21] = 1.0
    
    return features
