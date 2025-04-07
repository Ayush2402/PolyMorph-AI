import random
from typing import List, Optional, Dict, Any
import selfies as sf
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
import numpy as np

class PolymerGenerator:
    def __init__(self):
        self.kill_switch_groups = {
            "ester": "C(=O)O",
            "amide": "C(=O)N",
            "azide": "N=N=N",
            "peroxide": "OO",
            "thiol": "S",
            "photolabile": "C(=O)C"
        }
        
        self.domain_templates = {
            "agriculture": {
                "groups": ["C(=O)O", "C(=O)N"],
                "triggers": ["pH < 5", "temperature > 30", "microbial_presence"],
                "base_structures": ["CCO", "CC(=O)O", "CC(=O)N"]
            },
            "water": {
                "groups": ["S", "N=N=N"],
                "triggers": ["pH > 8", "UV_light", "oxidation"],
                "base_structures": ["CS", "CN=N=N", "CCO"]
            },
            "urban": {
                "groups": ["C(=O)O", "OO"],
                "triggers": ["temperature > 40", "mechanical_stress", "moisture"],
                "base_structures": ["CCO", "COO", "CC(=O)O"]
            }
        }

    def generate_polymer(self, domain: str, trigger_condition: Optional[str] = None) -> str:
        """
        Generate a polymer structure based on the specified domain and trigger condition.
        
        Args:
            domain: The domain for which to generate the polymer (agriculture, water, urban)
            trigger_condition: Optional trigger condition for degradation
            
        Returns:
            SMILES string of the generated polymer
        """
        if domain not in self.domain_templates:
            raise ValueError(f"Invalid domain: {domain}. Must be one of {list(self.domain_templates.keys())}")
        
        template = self.domain_templates[domain]
        base = random.choice(template["base_structures"])
        
        # Add random number of kill switch groups
        num_groups = random.randint(1, 3)
        groups = random.sample(template["groups"], min(num_groups, len(template["groups"])))
        
        # Create polymer structure
        structure = base
        for group in groups:
            if random.random() > 0.5:
                structure += group
        
        # Ensure valid SMILES
        mol = Chem.MolFromSmiles(structure)
        if mol is None:
            return self.generate_polymer(domain, trigger_condition)
        
        # Add hydrogens and clean up
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        return Chem.MolToSmiles(mol)

    def check_kill_switch(self, smiles: str) -> bool:
        """
        Check if a given SMILES string contains any kill switch groups.
        
        Args:
            smiles: SMILES string to check
            
        Returns:
            Boolean indicating if kill switch groups are present
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        
        for group_name, group_smarts in self.kill_switch_groups.items():
            pattern = Chem.MolFromSmarts(group_smarts)
            if mol.HasSubstructMatch(pattern):
                return True
        
        return False

    def predict_degradability(self, smiles: str) -> float:
        """
        Predict the degradability score of a polymer.
        
        Args:
            smiles: SMILES string of the polymer
            
        Returns:
            Degradability score between 0 and 1
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
        
        # Count degradable groups
        degradable_count = 0
        for group_smarts in self.kill_switch_groups.values():
            pattern = Chem.MolFromSmarts(group_smarts)
            matches = mol.GetSubstructMatches(pattern)
            degradable_count += len(matches)
        
        # Calculate base score from number of degradable groups
        base_score = min(0.8, degradable_count * 0.2)
        
        # Adjust score based on molecular weight (smaller molecules degrade faster)
        mw = Descriptors.ExactMolWt(mol)
        mw_factor = max(0, 1 - (mw / 1000))  # Normalize to 1000 g/mol
        mw_score = mw_factor * 0.2
        
        # Final score
        final_score = min(1.0, base_score + mw_score)
        
        return final_score

    def get_molecule_properties(self, smiles: str) -> Dict[str, Any]:
        """
        Get detailed properties of a molecule for visualization.
        
        Args:
            smiles: SMILES string of the molecule
            
        Returns:
            Dictionary containing molecule properties
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {}
        
        # Get atom information
        atoms = []
        for atom in mol.GetAtoms():
            atoms.append({
                "symbol": atom.GetSymbol(),
                "atomic_num": atom.GetAtomicNum(),
                "formal_charge": atom.GetFormalCharge(),
                "is_aromatic": atom.GetIsAromatic(),
                "hybridization": str(atom.GetHybridization())
            })
        
        # Get bond information
        bonds = []
        for bond in mol.GetBonds():
            bonds.append({
                "bond_type": str(bond.GetBondType()),
                "is_conjugated": bond.GetIsConjugated(),
                "is_aromatic": bond.GetIsAromatic(),
                "begin_atom_idx": bond.GetBeginAtomIdx(),
                "end_atom_idx": bond.GetEndAtomIdx()
            })
        
        # Get ring information using GetRingInfo
        ring_info = Chem.GetSymmSSSR(mol)
        num_rings = len(ring_info) if ring_info else 0
        
        return {
            "atoms": atoms,
            "bonds": bonds,
            "num_rings": num_rings,
            "molecular_weight": Descriptors.ExactMolWt(mol),
            "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "num_h_donors": Descriptors.NumHDonors(mol),
            "num_h_acceptors": Descriptors.NumHAcceptors(mol)
        } 