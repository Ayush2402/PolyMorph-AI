import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv, global_mean_pool
from rdkit import Chem
import numpy as np
from rdkit.Chem import Descriptors, rdMolDescriptors

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
    def __init__(self, model_path: str = None):
        self.converter = MoleculeToGraph()
        self.model = DegradabilityGNN(num_node_features=6)
        
        if model_path:
            try:
                if torch.cuda.is_available():
                    self.model.load_state_dict(torch.load(model_path))
                else:
                    self.model.load_state_dict(torch.load(model_path, map_location=torch.device('cpu')))
            except FileNotFoundError:
                print(f"Warning: Model file {model_path} not found. Using default model weights.")
        self.model.eval()
        
        # Define kill switch groups
        self.kill_switch_groups = {
            'ester': '[CX3](=O)[OX2H0][#6]',  # Ester group
            'amide': '[CX3](=O)[NX3H0][#6]',  # Amide group
            'anhydride': '[CX3](=O)[OX2][CX3](=O)',  # Anhydride group
            'carbonate': '[OX2][CX3](=O)[OX2]',  # Carbonate group
            'acetal': '[#6][OX2][#6]([#6])[OX2][#6]',  # Acetal group
            'imine': '[#6][NX2]=[#6]',  # Imine group
            'disulfide': '[SX2][SX2]',  # Disulfide bond
            'thioester': '[CX3](=O)[SX2][#6]',  # Thioester group
            'phosphoester': '[PX4](=O)([OX2])[OX2][#6]',  # Phosphoester group
            'siloxane': '[SiX4]([#6])[OX2][SiX4]([#6])[OX2]',  # Siloxane group
            'thiol': '[SX2H]',  # Thiol group
            'azide': '[NX3]=[NX2]=[NX1]'  # Azide group
        }
        
        # Define degradable functional groups
        self.degradable_groups = {
            'ester': '[CX3](=O)[OX2H0][#6]',  # Ester group
            'amide': '[CX3](=O)[NX3H0][#6]',  # Amide group
            'anhydride': '[CX3](=O)[OX2][CX3](=O)',  # Anhydride group
            'carbonate': '[OX2][CX3](=O)[OX2]',  # Carbonate group
            'acetal': '[#6][OX2][#6]([#6])[OX2][#6]',  # Acetal group
            'imine': '[#6][NX2]=[#6]',  # Imine group
            'disulfide': '[SX2][SX2]',  # Disulfide bond
            'thioester': '[CX3](=O)[SX2][#6]',  # Thioester group
            'phosphoester': '[PX4](=O)([OX2])[OX2][#6]',  # Phosphoester group
            'siloxane': '[SiX4]([#6])[OX2][SiX4]([#6])[OX2]'  # Siloxane group
        }
        
        # Define synthesizability thresholds
        self.synthesizability_thresholds = {
            'molecular_weight': 1000,  # g/mol
            'stereocenters': 3,
            'rings': 3,
            'rotatable_bonds': 10,
            'sp3_fraction': 0.3
        }

    def has_kill_switch(self, smiles: str) -> bool:
        """Check if the polymer has a kill switch group."""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
            
        for group_name, pattern in self.kill_switch_groups.items():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
                return True
        return False

    def predict_degradability(self, smiles: str) -> float :
        """
        Compute intrinsic polymer properties and an overall degradation score.
        
        Properties computed:
            - molecular_weight (Da): Determined from RDKit Descriptors.ExactMolWt.
              Normalized to [5e3, 1e6].
            - crystallinity (%): Approximated inversely by the number of rotatable bonds.
              (Fewer rotatable bonds imply higher crystallinity.)
            - hydrolyzable_bonds (ratio): Fraction of bonds matching degradable groups (ester, amide, acetal).
            - reactive_groups (ratio): Fraction of atoms participating in reactive groups (peroxide, thiol, azide, imine).
            - surface_area (m²/g): Estimated from the Topological Polar Surface Area (TPSA) and molecular weight.
              Normalized to the range [0.1, 10] m²/g.
            - hydrophilicity (%): Fraction of atoms in hydrophilic groups ([OH], [NH2], [C](=O)[O]).
        
        These normalized scores are then combined (via weighted sum) to yield an overall score between 0 and 1.
        
        Args:
            smiles: SMILES string of the polymer
            
        Returns:
            Dictionary with computed properties and overall score.
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        # 1. Molecular Weight
        mw = Descriptors.ExactMolWt(mol)
        # Normalize to 0-1 given range [5e3, 1e6]
        mw_norm = (mw - 5000) / (1e6 - 5000)
        mw_norm = max(0.0, min(1.0, mw_norm))
        
        # 2. Crystallinity
        # Approximate crystallinity (a higher number of rotatable bonds generally reduces crystallinity)
        rot_bonds = Descriptors.NumRotatableBonds(mol)
        # Using a simple formula: crystallinity = 1/(1 + rot_bonds)
        cryst = 1.0 / (1.0 + rot_bonds)
        # For reporting as a percentage:
        cryst_percent = round(cryst * 100, 2)
        
        # 3. Hydrolyzable Bonds Ratio
        hydrolyzable_keys = ["ester", "amide", "acetal"]
        hydrolyzable_count = 0
        for key in hydrolyzable_keys:
            pattern = Chem.MolFromSmarts(self.kill_switch_groups[key])
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                hydrolyzable_count += len(matches)
        total_bonds = mol.GetNumBonds()
        hydrolyzable_ratio = hydrolyzable_count / total_bonds if total_bonds > 0 else 0.0
        
        # 4. Reactive Groups Ratio
        # Consider these groups as reactive: thiol, azide, imine
        reactive_keys = ["thiol", "azide", "imine"]
        reactive_count = 0
        for key in reactive_keys:
            pattern = Chem.MolFromSmarts(self.kill_switch_groups[key])
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                reactive_count += len(matches)
        total_atoms = mol.GetNumAtoms()
        reactive_ratio = reactive_count / total_atoms if total_atoms > 0 else 0.0
        
        # 5. Surface Area (m²/g)
        try:
            from rdkit.Chem import rdMolDescriptors
            tpsa = rdMolDescriptors.CalcTPSA(mol)  # in Å²
            # Estimate surface area: Convert A² to m² (1 A² = 1e-20 m²) and normalize by molecular weight (g/mol ≈ Da)
            surface_area = (tpsa * 1e-20) / (mw / 1e3)  # m²/g
            # Normalize surface area given desired range [0.1, 10] m²/g
            sa_norm = (surface_area - 0.1) / (10 - 0.1)
            sa_norm = max(0.0, min(1.0, sa_norm))
        except Exception:
            surface_area = 0.0
            sa_norm = 0.0
        
        # 6. Hydrophilicity (%)
        # Count hydrophilic groups: hydroxyl ([OH]), amine ([NH2]), and carboxyl ([C](=O)[O])
        hydrophilic_smarts = ["[OH]", "[NH2]", "[C](=O)[O]"]
        hydro_count = 0
        for smarts in hydrophilic_smarts:
            pattern = Chem.MolFromSmarts(smarts)
            if pattern:
                matches = mol.GetSubstructMatches(pattern)
                hydro_count += len(matches)
        hydrophilicity_ratio = hydro_count / total_atoms if total_atoms > 0 else 0.0
        hydrophilicity_percent = hydrophilicity_ratio * 100
        
        # Overall Score Calculation:
        # Define weights for each property. Notice that for molecular weight and crystallinity lower values
        # (i.e., lower molecular weight and lower crystallinity) may favor degradation; hence their contributions are taken as (1 - normalized_value).
        weights = {
            "molecular_weight": 0.15,
            "crystallinity": 0.15,
            "hydrolyzable_bonds": 0.20,
            "reactive_groups": 0.15,
            "surface_area": 0.20,
            "hydrophilicity": 0.15
        }
        overall_score = (
            weights["molecular_weight"] * (1 - mw_norm) +
            weights["crystallinity"] * (1 - cryst) +
            weights["hydrolyzable_bonds"] * hydrolyzable_ratio +
            weights["reactive_groups"] * reactive_ratio +
            weights["surface_area"] * sa_norm +
            weights["hydrophilicity"] * hydrophilicity_ratio
        )
        overall_score = max(0.0, min(1.0, overall_score))
        
        return overall_score

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
    
