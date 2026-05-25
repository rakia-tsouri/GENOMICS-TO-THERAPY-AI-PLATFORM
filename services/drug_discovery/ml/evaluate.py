import os
import torch
import torch.nn as nn
from torch_geometric.loader import DataLoader
from model import DrugProteinGNN
from dataset import ProteinLigandDataset
import numpy as np
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score

def evaluate():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")
    csv_file = os.path.join(data_dir, "train_combined.csv")
    protein_dir = os.path.join(data_dir, "PDBbind_v2020_refined", "refined-set")
    model_save_path = os.path.join(base_dir, "models_saved", "gnn_v1.pt")
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    if not os.path.exists(model_save_path):
        print(f"Error: Model not found at {model_save_path}")
        return

    # Initialisation du Dataset
    print("Loading dataset...")
    dataset = ProteinLigandDataset(root=data_dir, csv_file=csv_file, protein_dir=protein_dir)
    
    # Split Train/Test (80/20) with a fixed generator for reproducibility
    generator = torch.Generator().manual_seed(42)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    _, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size], generator=generator)

    test_loader = DataLoader(test_dataset, batch_size=128, shuffle=False, follow_batch=['prot_x', 'drug_x'])

    # Modèle
    model = DrugProteinGNN().to(device)
    model.load_state_dict(torch.load(model_save_path, map_location=device))
    model.eval()

    print("Evaluating model on test dataset...")
    all_preds = []
    all_labels = []

    with torch.no_grad():
        for prot_batch, drug_batch, labels in test_loader:
            prot_batch = prot_batch.to(device)
            drug_batch = drug_batch.to(device)
            labels = labels.to(device)
            
            pred = model(prot_batch, drug_batch)
            all_preds.extend(pred.cpu().numpy().flatten())
            all_labels.extend(labels.cpu().numpy().flatten())

    all_preds = np.array(all_preds)
    all_labels = np.array(all_labels)

    mse = mean_squared_error(all_labels, all_preds)
    rmse = np.sqrt(mse)
    mae = mean_absolute_error(all_labels, all_preds)
    r2 = r2_score(all_labels, all_preds)
    corr, _ = pearsonr(all_labels, all_preds)

    print("\n" + "="*40)
    print("           EVALUATION RESULTS")
    print("="*40)
    print(f"Test Samples: {len(all_labels)}")
    print(f"Mean Squared Error (MSE):     {mse:.4f}")
    print(f"Root Mean Squared Error (RMSE): {rmse:.4f}")
    print(f"Mean Absolute Error (MAE):     {mae:.4f}")
    print(f"R-squared (R2) Score:         {r2:.4f}")
    print(f"Pearson Correlation (r):      {corr:.4f}")
    print("="*40)

    print("\nExample Predictions vs True Values:")
    print("True Affinity | Predicted Affinity | Difference")
    print("-" * 47)
    for i in range(min(10, len(all_labels))):
        diff = all_preds[i] - all_labels[i]
        print(f"{all_labels[i]:13.4f} | {all_preds[i]:18.4f} | {diff:+10.4f}")

if __name__ == "__main__":
    evaluate()
