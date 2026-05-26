import os
import torch
import torch.nn as nn
from torch_geometric.loader import DataLoader
from model import DrugProteinGNN
from dataset import ProteinLigandDataset

def train(num_epochs=50, lr=1e-3, batch_size=128):
    # Configuration des chemins
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    data_dir = os.path.join(base_dir, "data")
    csv_file = os.path.join(data_dir, "train_combined.csv")
    protein_dir = os.path.join(data_dir, "PDBbind_v2020_refined", "refined-set")
    model_save_path = os.path.join(base_dir, "models_saved", "gnn_v1.pt")
    checkpoint_path = os.path.join(base_dir, "models_saved", "checkpoint.pt")
    
    # Création du dossier de sauvegarde s'il n'existe pas
    os.makedirs(os.path.dirname(model_save_path), exist_ok=True)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Using device: {device}")

    # Initialisation du Dataset
    print("Loading dataset (this might take a while the first time)...")
    dataset = ProteinLigandDataset(root=data_dir, csv_file=csv_file, protein_dir=protein_dir)
    
    # Split Train/Test (80/20)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_dataset, test_dataset = torch.utils.data.random_split(dataset, [train_size, test_size])

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True, follow_batch=['prot_x', 'drug_x'])
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False, follow_batch=['prot_x', 'drug_x'])

    # Modèle
    model = DrugProteinGNN().to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = nn.MSELoss()

    # Chargement du checkpoint si disponible
    start_epoch = 0
    if os.path.exists(checkpoint_path):
        print(f"Loading checkpoint from {checkpoint_path}...")
        checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
        model.load_state_dict(checkpoint['model_state_dict'])
        optimizer.load_state_dict(checkpoint['optimizer_state_dict'])
        start_epoch = checkpoint['epoch'] + 1
        print(f"Resuming from epoch {start_epoch}")
    elif os.path.exists(model_save_path):
        print(f"Loading weights from {model_save_path} (no checkpoint found)...")
        model.load_state_dict(torch.load(model_save_path, map_location=device, weights_only=True))
        # Si on charge seulement les poids, on suppose qu'on a fini 5 époques comme dit par l'utilisateur
        start_epoch = 5 
        print(f"Starting from epoch {start_epoch + 1}")

    print(f"Starting training on {len(train_dataset)} samples...")
    for epoch in range(start_epoch, num_epochs):
        model.train()
        total_loss = 0
        for prot_batch, drug_batch, labels in train_loader:
            prot_batch = prot_batch.to(device)
            drug_batch = drug_batch.to(device)
            labels = labels.to(device)
            
            optimizer.zero_grad()
            pred = model(prot_batch, drug_batch)
            loss = criterion(pred, labels)
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        
        avg_loss = total_loss / len(train_loader)
        print(f"Epoch {epoch+1}/{num_epochs} - Loss: {avg_loss:.4f}")

        # Sauvegarde du checkpoint après chaque époque
        torch.save({
            'epoch': epoch,
            'model_state_dict': model.state_dict(),
            'optimizer_state_dict': optimizer.state_dict(),
            'loss': avg_loss,
        }, checkpoint_path)

        # Validation rapide
        if (epoch + 1) % 5 == 0:
            model.eval()
            val_loss = 0
            with torch.no_grad():
                for prot_batch, drug_batch, labels in test_loader:
                    prot_batch = prot_batch.to(device)
                    drug_batch = drug_batch.to(device)
                    labels = labels.to(device)
                    pred = model(prot_batch, drug_batch)
                    val_loss += criterion(pred, labels).item()
            print(f" >> Validation Loss: {val_loss / len(test_loader):.4f}")

    # Sauvegarde finale du modèle
    torch.save(model.state_dict(), model_save_path)
    print(f"Training completed. Final model saved to {model_save_path}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Train Drug Discovery GNN model")
    parser.add_argument("--epochs", type=int, default=50, help="Number of training epochs")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--batch-size", type=int, default=128, help="Batch size")
    args = parser.parse_args()
    
    train(num_epochs=args.epochs, lr=args.lr, batch_size=args.batch_size)
