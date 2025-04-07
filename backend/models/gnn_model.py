import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool
from rdkit import Chem
import numpy as np
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

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

    def predict_degradability(self, smiles: str) -> float:
        """
        Predict the degradability score using the GNN model.
        """
        graph = self.converter.convert(smiles)
        if graph is None:
            return 0.0

        with torch.no_grad():
            pred = self.model(graph['nodes'], graph['edges'], graph['batch'])
            return float(pred[0])

    def predict_synthesizability(self, smiles: str) -> float:
        """
        Predict the synthesizability of a molecule based on its structural complexity.
        Returns a score between 0 and 1, where 1 means easily synthesizable.
        """
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 0.0

            # Calculate molecular descriptors
            mw = Descriptors.ExactMolWt(mol)
            rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            rings = rdMolDescriptors.CalcNumRings(mol)
            stereo_centers = len(Chem.FindMolChiralCenters(mol))
            
            # Calculate sp3 carbon fraction manually
            num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
            num_sp3_carbons = sum(1 for atom in mol.GetAtoms() 
                                if atom.GetSymbol() == 'C' and atom.GetHybridization() == Chem.HybridizationType.SP3)
            sp3_carbons = num_sp3_carbons / max(num_carbons, 1)  # Avoid division by zero
            
            # Normalize and combine factors
            mw_score = np.clip(1.0 - (mw / 500.0), 0, 1)  # Penalize high molecular weight
            complexity_score = np.clip(1.0 - (rings * 0.1 + stereo_centers * 0.15), 0, 1)  # Penalize complexity
            flexibility_score = np.clip(1.0 - (rotatable_bonds * 0.05), 0, 1)  # Penalize high flexibility
            sp3_score = sp3_carbons  # Higher sp3 fraction is generally easier to synthesize
            
            # Weighted combination of scores
            synthesizability_score = (
                0.3 * mw_score +
                0.3 * complexity_score +
                0.2 * flexibility_score +
                0.2 * sp3_score
            )
            
            return float(np.clip(synthesizability_score, 0, 1))
            
        except Exception as e:
            print(f"Error calculating synthesizability: {e}")
            return 0.0

    def predict(self, smiles: str) -> dict:
        """
        Predict both degradability and synthesizability scores for a molecule.
        """
        degradability = self.predict_degradability(smiles)
        synthesizability = self.predict_synthesizability(smiles)
        
        return {
            "degradability": degradability,
            "synthesizability": synthesizability
        } 