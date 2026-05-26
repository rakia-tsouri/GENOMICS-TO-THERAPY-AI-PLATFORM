import torch
from torch_geometric.loader import DataLoader
from model import DrugProteinGNN
from dataset import ProteinLigandDataset
from evaluate import calculate_metrics
import os
import numpy as np

def run_evaluation():
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")
    model_path = os.path.join(base_dir, "models_saved", "gnn_v1.pt")
    
    # Configuration des chemins pour le Dataset
    csv_file = os.path.join(data_dir, "train_data.csv")
    protein_dir = os.path.join(data_dir, "PDBbind_v2020_refined", "refined-set")

    # 1. Charger le Dataset
    print("Loading dataset...")
    dataset = ProteinLigandDataset(root=data_dir, csv_file=csv_file, protein_dir=protein_dir)
    
    # Split Train/Test identique à train.py
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    _, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size], 
                                                  generator=torch.Generator().manual_seed(42))
    
    loader = DataLoader(test_dataset, batch_size=32, shuffle=False, follow_batch=['prot_x', 'drug_x'])

    # 2. Charger le Modèle
    model = DrugProteinGNN().to(device)
    if not os.path.exists(model_path):
        print(f"Error: Model not found at {model_path}")
        return

    model.load_state_dict(torch.load(model_path, map_location=device, weights_only=True))
    model.eval()

    all_preds = []
    all_labels = []

    print(f"Evaluating model on {len(test_dataset)} samples...")
    with torch.no_grad():
        for prot_batch, drug_batch, labels in loader:
            prot_batch = prot_batch.to(device)
            drug_batch = drug_batch.to(device)
            
            preds = model(prot_batch, drug_batch)
            # S'assurer que preds et labels sont concaténables (flatten)
            all_preds.append(preds.cpu().numpy().reshape(-1))
            all_labels.append(labels.cpu().numpy().reshape(-1))

    y_pred = np.concatenate(all_preds)
    y_true = np.concatenate(all_labels)

    # 3. Calculer les métriques
    metrics = calculate_metrics(y_true, y_pred)
    
    print("\n" + "="*30)
    print("   RESULTATS DE L'EVALUATION")
    print("="*30)
    print(f"Pearson R (Corrélation): {metrics['pearson_r']:.4f}")
    print(f"RMSE (Erreur moyenne):   {metrics['rmse']:.4f}")
    print(f"AUC-ROC (Seuil > 6):     {metrics['auc_roc']:.4f}")
    print("="*30)

if __name__ == "__main__":
    run_evaluation()
