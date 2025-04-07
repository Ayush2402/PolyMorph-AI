import torch
import pandas as pd
from pathlib import Path
from tqdm import tqdm
import numpy as np
from sklearn.model_selection import train_test_split

from ..models.gnn_model import DegradabilityGNN, MoleculeToGraph

def load_data(data_path: str = "data/train.csv"):
    df = pd.read_csv(data_path)
    return df["smiles"].values, df["degradability_score"].values

def prepare_data(smiles_list, scores, converter):
    data = []
    for smiles, score in zip(smiles_list, scores):
        graph = converter.convert(smiles)
        if graph is not None:
            graph["y"] = torch.tensor([score], dtype=torch.float)
            data.append(graph)
    return data

def train_model(model, train_data, val_data, num_epochs=50, lr=0.001):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model = model.to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=lr)
    criterion = torch.nn.MSELoss()
    
    best_val_loss = float("inf")
    best_model_state = None
    
    for epoch in range(num_epochs):
        # Training
        model.train()
        total_loss = 0
        for batch in train_data:
            optimizer.zero_grad()
            nodes = batch["nodes"].to(device)
            edges = batch["edges"].to(device)
            batch_idx = batch["batch"].to(device)
            y = batch["y"].to(device)
            
            pred = model(nodes, edges, batch_idx)
            loss = criterion(pred, y)
            
            loss.backward()
            optimizer.step()
            total_loss += loss.item()
        
        # Validation
        model.eval()
        val_loss = 0
        with torch.no_grad():
            for batch in val_data:
                nodes = batch["nodes"].to(device)
                edges = batch["edges"].to(device)
                batch_idx = batch["batch"].to(device)
                y = batch["y"].to(device)
                
                pred = model(nodes, edges, batch_idx)
                val_loss += criterion(pred, y).item()
        
        print(f"Epoch {epoch+1}/{num_epochs}")
        print(f"Training Loss: {total_loss/len(train_data):.4f}")
        print(f"Validation Loss: {val_loss/len(val_data):.4f}")
        
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_model_state = model.state_dict()
    
    return best_model_state

def main():
    # Create output directory
    model_dir = Path("data/models")
    model_dir.mkdir(exist_ok=True)
    
    # Load and prepare data
    smiles, scores = load_data()
    converter = MoleculeToGraph()
    data = prepare_data(smiles, scores, converter)
    
    # Split data
    train_data, val_data = train_test_split(data, test_size=0.2, random_state=42)
    
    # Initialize and train model
    model = DegradabilityGNN(num_node_features=6)
    best_model_state = train_model(model, train_data, val_data)
    
    # Save model
    torch.save(best_model_state, model_dir / "gnn_model.pt")
    print("Model saved successfully!")

if __name__ == "__main__":
    main() 