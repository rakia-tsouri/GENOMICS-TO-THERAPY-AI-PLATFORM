import os
import pandas as pd
import numpy as np
import csv
from tqdm import tqdm

def extract():
    # Chemins
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    pdb_list_path = os.path.join(base_dir, "data", "existing_pdbs.txt")
    bindingdb_path = os.path.join(base_dir, "data", "BindingDB_All_202605_tsv", "BindingDB_All.tsv")
    output_path = os.path.join(base_dir, "data", "bindingdb_extracted.csv")

    if not os.path.exists(bindingdb_path):
        print(f"Error: {bindingdb_path} not found.")
        return

    # 1. Charger la liste des PDB que nous avons déjà
    print("Loading existing PDB IDs...")
    with open(pdb_list_path, 'r') as f:
        # Stocker en minuscules pour la comparaison
        existing_pdbs = set(line.strip().lower() for line in f if line.strip())
    
    print(f"Searching BindingDB for matches with {len(existing_pdbs)} proteins...")

    new_data = []
    
    # 2. Lire BindingDB ligne par ligne
    try:
        with open(bindingdb_path, 'r', encoding='utf-8', errors='ignore') as f:
            reader = csv.reader(f, delimiter='\t')
            header = next(reader)
            
            # On vérifie les index dynamiquement au cas où
            try:
                idx_smiles = header.index('Ligand SMILES')
                idx_ic50 = header.index('IC50 (nM)')
                idx_ki = header.index('Ki (nM)')
                # Colonnes contenant des PDB IDs
                pdb_cols = [
                    header.index('PDB ID(s) for Ligand-Target Complex'),
                    header.index('PDB ID(s) of Target Chain 1')
                ]
            except ValueError:
                idx_smiles, idx_ki, idx_ic50 = 1, 8, 9
                pdb_cols = [30, 41]

            found = 0
            for row in tqdm(reader, desc="Scanning BindingDB (8GB)", unit=" lines"):
                if len(row) <= max(idx_smiles, idx_ic50, idx_ki, max(pdb_cols)): continue
                
                # Chercher un match dans les colonnes PDB identifiées
                match = None
                for col_idx in pdb_cols:
                    pdb_str = row[col_idx].lower().replace(' ', '')
                    if not pdb_str: continue
                    
                    row_pdbs = pdb_str.split(',')
                    for p in row_pdbs:
                        if p in existing_pdbs:
                            match = p
                            break
                    if match: break
                
                if match:
                    smiles = row[idx_smiles]
                    val_str = row[idx_ic50] if row[idx_ic50] else row[idx_ki]
                    
                    if smiles and val_str:
                        try:
                            # Nettoyer la valeur numérique
                            val_str = val_str.replace('>', '').replace('<', '').replace(' ', '').replace('nM', '')
                            val_nM = float(val_str)
                            if 0.0001 < val_nM < 100000000:
                                pIC50 = -np.log10(val_nM * 1e-9)
                                new_data.append([match, smiles, pIC50])
                                found += 1
                        except:
                            continue
                
                if found >= 100000: 
                    print(f"\nReached {found} samples limit. Stopping.")
                    break

        # 3. Sauvegarder
        if new_data:
            df_new = pd.DataFrame(new_data, columns=['pdb_id', 'smiles', 'affinity'])
            df_new = df_new.drop_duplicates(subset=['pdb_id', 'smiles'])
            df_new.to_csv(output_path, index=False)
            print(f"\nSUCCESS: Extracted {len(df_new)} unique interactions to {output_path}")
        else:
            print("\nNo matches found. Checking first 100 rows for debug...")
            # Petit bloc de debug pour voir ce qui cloche
            with open(bindingdb_path, 'r', encoding='utf-8', errors='ignore') as f:
                reader = csv.reader(f, delimiter='\t')
                next(reader)
                for i in range(10):
                    row = next(reader)
                    print(f"Row {i} - Col 30: {row[30]}, Col 41: {row[41]}")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    extract()
