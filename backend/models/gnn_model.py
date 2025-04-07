import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool
from rdkit import Chem
import numpy as np

class DegradabilityGNN(torch.nn.Module):
    def __init__(self, num_node_features):
        super(DegradabilityGNN, self).__init__()
        self.conv1 = GCNConv(num_node_features, 64)
        self.conv2 = GCNConv(64, 32)
        self.conv3 = GCNConv(32, 16)
        self.fc = torch.nn.Linear(16, 1)

    def forward(self, x, edge_index, batch):
        x = F.relu(self.conv1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))
        x = F.relu(self.conv3(x, edge_index))
        x = global_mean_pool(x, batch)
        x = self.fc(x)
        return torch.sigmoid(x)

class MoleculeToGraph:
    def __init__(self):
        self.atom_features = {
            'C': [1, 0, 0, 0, 0, 0],
            'N': [0, 1, 0, 0, 0, 0],
            'O': [0, 0, 1, 0, 0, 0],
            'S': [0, 0, 0, 1, 0, 0],
            'P': [0, 0, 0, 0, 1, 0],
            'Other': [0, 0, 0, 0, 0, 1]
        }

    def convert(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None

        # Node features
        nodes = []
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            features = self.atom_features.get(symbol, self.atom_features['Other'])
            nodes.append(features)

        # Edge indices
        edges = []
        for bond in mol.GetBonds():
            i = bond.GetBeginAtomIdx()
            j = bond.GetEndAtomIdx()
            # Add both directions for undirected graph
            edges.append([i, j])
            edges.append([j, i])

        return {
            'nodes': torch.tensor(nodes, dtype=torch.float),
            'edges': torch.tensor(edges, dtype=torch.long).t().contiguous(),
            'batch': torch.zeros(len(nodes), dtype=torch.long)
        }

class DegradabilityPredictor:
    def __init__(self, model_path=None):
        self.converter = MoleculeToGraph()
        self.model = DegradabilityGNN(num_node_features=6)
        if model_path and torch.cuda.is_available():
            self.model.load_state_dict(torch.load(model_path))
        elif model_path:
            self.model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
        self.model.eval()

    def predict(self, smiles):
        graph = self.converter.convert(smiles)
        if graph is None:
            return 0.0

        with torch.no_grad():
            pred = self.model(graph['nodes'], graph['edges'], graph['batch'])
            return float(pred[0]) 