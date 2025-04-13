import random
from typing import List, Optional, Dict, Any
import selfies as sf
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
import numpy as np

class PolymerGenerator:
    def __init__(self):
        self.kill_switch_groups = {
            "ester": "[CX3](=O)[OX2H0]", 
            "amide": "[CX3](=O)[NX3]", 
            "acetal": "[OX2][CX4]([OX2])", 
            "hydrazone": "[CX3]=[NX2][NX3]",
            "imine": "[CX3]=[NX2]",
            "thiol": "[SX2H]",
            "azide": "[NX2]=[NX2]=[NX1]", 
            "anhydride": "[CX3](=O)[OX2][CX3]=O", 
            "carbonate": "[OX2][CX3](=O)[OX2]", 
            "disulfide": "[SX2D2][SX2D2]", 
            "peroxide": "[OX2D2][OX2D2]", 
            "thioketal": "[SX2][CX4]([SX2])",  
            "boronic_ester": "[BX3]([OX2])[OX2]",  
            "siloxane": "[OX2][Si][OX2]"  
        }
        
        self.construction_groups = {
            "ester": "C(=O)O",
            "amide": "C(=O)N",
            "acetal": "OC(O)C",
            "hydrazone": "C=NN",
            "imine": "C=N",
            "thiol": "SH",
            "azide": "N=N=N",
            "anhydride": "C(=O)OC(=O)",
            "carbonate": "OC(=O)O",
            "disulfide": "SS",
            "peroxide": "OO",
            "thioketal": "SC(S)C",
            "boronic_ester": "B(O)O",
            "siloxane": "OSiO"
        }
        
        self.domain_templates = {
            "agriculture": {
                "groups": [
                    self.construction_groups["ester"],
                    self.construction_groups["amide"],
                    self.construction_groups["acetal"],
                    self.construction_groups["hydrazone"],
                    self.construction_groups["imine"]
                ],
                "triggers": ["pH < 5", "temperature > 30", "microbial_presence"],
                "base_structures": [
                    "CCO",      # Simple alcohol
                    "CC(=O)O",  # Carboxylic acid
                    "CC(=O)N",  # Simple amide
                    "CCOC",     # Ether
                    "CC=O"      # Aldehyde
                ]
            },
            "water": {
                "groups": [
                    self.construction_groups["thiol"],
                    self.construction_groups["azide"],
                    self.construction_groups["anhydride"],
                    self.construction_groups["carbonate"],
                    self.construction_groups["disulfide"]
                ],
                "triggers": ["pH > 8", "UV_light", "oxidation"],
                "base_structures": [
                    "CS",       # Sulfide
                    "CN=N=N",   # Azide derivative
                    "CCO",      # Alcohol
                    "CCOC",     # Ether
                    "CSC"       # Thioether
                ]
            },
            "urban": {
                "groups": [
                    self.construction_groups["ester"],
                    self.construction_groups["peroxide"],
                    self.construction_groups["thioketal"],
                    self.construction_groups["boronic_ester"],
                    self.construction_groups["siloxane"]
                ],
                "triggers": ["temperature > 40", "mechanical_stress", "moisture"],
                "base_structures": [
                    "CCO",      # Alcohol
                    "COO",      # Peroxide
                    "CC(=O)O",  # Carboxylic acid
                    "CCOC",     # Ether
                    "CC(C)O"    # Secondary alcohol
                ]
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
        
        # Create polymer structure with proper connections
        structure = base
        for group in groups:
            if random.random() > 0.5:
                # Handle different connection cases
                if structure.endswith('O') and group.startswith('C'):
                    # Connect oxygen to carbon
                    structure = structure + group
                elif structure.endswith('N') and group.startswith('C'):
                    # Connect nitrogen to carbon
                    structure = structure + group
                elif structure.endswith('C') and (group.startswith('O') or group.startswith('N') or group.startswith('S')):
                    # Connect carbon to heteroatom
                    structure = structure + group
                elif structure.endswith('C') and group.startswith('C'):
                    # Connect carbon to carbon, remove one carbon to avoid duplicates
                    structure = structure[:-1] + group
                else:
                    # Default case: add a carbon linker
                    structure = structure + 'C' + group
        
        # Ensure valid SMILES
        mol = Chem.MolFromSmiles(structure)
        if mol is None:
            return self.generate_polymer(domain, trigger_condition)
        
        try:
            # Add hydrogens and clean up
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)  # Add seed for reproducibility
            AllChem.MMFFOptimizeMolecule(mol)
            
            return Chem.MolToSmiles(mol)
        except Exception as e:
            # If embedding fails, try again with a new base
            return self.generate_polymer(domain, trigger_condition)

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
            if pattern is not None and mol.HasSubstructMatch(pattern):
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
            if pattern is not None:
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
        
    def generate_nomenclature(self, smiles: str, domain: str) -> Dict[str, Any]:
        """
        Generate systematic nomenclature for a polymer based on its structure.
        
        Args:
            smiles: SMILES string of the polymer
            domain: The domain for which the polymer was generated
            
        Returns:
            Dictionary containing nomenclature data
        """
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                "name": "Unknown Polymer",
                "formula": "Unknown",
                "functional_groups": [],
                "structure_type": "Unknown"
            }
        
        # Identify functional groups
        functional_groups = []
        for group_name, group_smarts in self.kill_switch_groups.items():
            pattern = Chem.MolFromSmarts(group_smarts)
            if mol.HasSubstructMatch(pattern):
                functional_groups.append(group_name)
        
        # Check for rings to determine structure type
        ring_info = Chem.GetSymmSSSR(mol)
        num_rings = len(ring_info) if ring_info else 0
        structure_type = "Cyclic" if num_rings > 0 else "Linear"
        
        # Generate molecular formula
        atom_counts = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
        
        # Order by C, H, and then alphabetically
        formula_parts = []
        if 'C' in atom_counts:
            formula_parts.append(f"C{atom_counts['C'] if atom_counts['C'] > 1 else ''}")
            del atom_counts['C']
        if 'H' in atom_counts:
            formula_parts.append(f"H{atom_counts['H'] if atom_counts['H'] > 1 else ''}")
            del atom_counts['H']
        
        # Add remaining elements alphabetically
        for symbol in sorted(atom_counts.keys()):
            formula_parts.append(f"{symbol}{atom_counts[symbol] if atom_counts[symbol] > 1 else ''}")
        
        formula = ''.join(formula_parts)
        
        # Generate systematic name
        domain_prefix = {
            'agriculture': 'Agro',
            'water': 'Hydro',
            'urban': 'Urbo'
        }.get(domain, '')
        
        # Generate base name from functional groups
        if not functional_groups:
            base_name = "polymer"
        elif len(functional_groups) == 1:
            base_name = f"poly{functional_groups[0]}"
        else:
            base_name = 'poly' + '-co-'.join(functional_groups)
        
        systematic_name = f"{domain_prefix}{structure_type.lower()}-{base_name}"
        
        return {
            "name": systematic_name,
            "formula": formula,
            "functional_groups": functional_groups,
            "structure_type": structure_type
        } 