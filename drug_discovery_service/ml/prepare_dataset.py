import os
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

def sdf_to_smiles(sdf_path):
    try:
        suppl = Chem.SDMolSupplier(sdf_path)
        mol = next(suppl)
        if mol:
            return Chem.MolToSmiles(mol)
    except Exception as e:
        print(f"Error processing {sdf_path}: {e}")
    return None

def prepare_pdbbind_dataset(index_file, refined_set_dir, output_csv):
    print(f"Reading index file: {index_file}")
    data = []
    
    with open(index_file, 'r') as f:
        lines = f.readlines()
        
    for line in tqdm(lines, desc="Processing PDBbind"):
        if line.startswith('#'):
            continue
            
        parts = line.split()
        if len(parts) < 5:
            continue
            
        pdb_id = parts[0]
        affinity = float(parts[3]) # -logKd/Ki
        
        # Path to the ligand SDF file
        ligand_path = os.path.join(refined_set_dir, pdb_id, f"{pdb_id}_ligand.sdf")
        
        if os.path.exists(ligand_path):
            smiles = sdf_to_smiles(ligand_path)
            if smiles:
                data.append({
                    'pdb_id': pdb_id,
                    'smiles': smiles,
                    'affinity': affinity
                })
        else:
            print(f"Warning: Ligand file not found for {pdb_id}")
            
    df = pd.DataFrame(data)
    df.to_csv(output_csv, index=False)
    print(f"Dataset saved to {output_csv}. Total samples: {len(df)}")

if __name__ == "__main__":
    # Define paths relative to the script location or project root
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")
    
    index_file = os.path.join(data_dir, "PDBbind_v2020_plain_text_index", "index", "INDEX_refined_data.2020")
    refined_set_dir = os.path.join(data_dir, "PDBbind_v2020_refined", "refined-set")
    output_csv = os.path.join(data_dir, "train_data.csv")
    
    if os.path.exists(index_file) and os.path.exists(refined_set_dir):
        prepare_pdbbind_dataset(index_file, refined_set_dir, output_csv)
    else:
        print("Error: Required directories or files not found.")
        print(f"Index file exists: {os.path.exists(index_file)}")
        print(f"Refined set dir exists: {os.path.exists(refined_set_dir)}")
