import pandas as pd
import os

def merge():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    base_dir = os.path.dirname(script_dir)
    data_dir = os.path.join(base_dir, "data")
    
    path1 = os.path.join(data_dir, "train_data.csv")
    path2 = os.path.join(data_dir, "bindingdb_extracted.csv")
    
    print("Merging files...")
    df1 = pd.read_csv(path1)
    df2 = pd.read_csv(path2)
    
    combined = pd.concat([df1, df2], ignore_index=True)
    
    # Supprimer les doublons
    before = len(combined)
    combined = combined.drop_duplicates(subset=['pdb_id', 'smiles'])
    after = len(combined)
    
    output_path = os.path.join(data_dir, "train_combined.csv")
    combined.to_csv(output_path, index=False)
    
    print(f"Original samples: {len(df1)}")
    print(f"New samples: {len(df2)}")
    print(f"Duplicates removed: {before - after}")
    print(f"FINAL COMBINED DATASET: {after} samples saved to {output_path}")

if __name__ == "__main__":
    merge()
