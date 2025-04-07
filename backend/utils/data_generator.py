import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import json
from pathlib import Path
import random

class PolymerDataGenerator:
    def __init__(self):
        # SMARTS patterns for matching
        self.kill_switch_groups = {
            # Common degradable groups with valid SMARTS patterns
            "ester": "[CX3](=O)[OX2H0]",  # C(=O)O
            "amide": "[CX3](=O)[NX3]",    # C(=O)N
            "acetal": "[OX2][CX4]([OX2])",  # OC(O)C
            "hydrazone": "[CX3]=[NX2][NX3]",  # C=NN
            "imine": "[CX3]=[NX2]",  # C=N
            "thiol": "[SX2H]",  # SH
            "azide": "[NX2]=[NX2]=[NX1]",  # N=N=N
            "anhydride": "[CX3](=O)[OX2][CX3]=O",  # C(=O)OC(=O)
            "carbonate": "[OX2][CX3](=O)[OX2]",  # OC(=O)O
            "disulfide": "[SX2D2][SX2D2]",  # SS
            "peroxide": "[OX2D2][OX2D2]",  # OO
            "thioketal": "[SX2][CX4]([SX2])",  # SC(S)C
            "boronic_ester": "[BX3]([OX2])[OX2]",  # B(O)O
            "siloxane": "[OX2][Si][OX2]"  # OSiO
        }
        
        # Simplified SMILES patterns for construction
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

    def generate_polymer(self, domain: str, trigger: str = None) -> dict:
        template = self.domain_templates[domain]
        base = random.choice(template["base_structures"])
        
        # Add random number of kill switch groups
        num_groups = random.randint(1, 3)
        groups = random.sample(template["groups"], min(num_groups, len(template["groups"])))
        
        # Create polymer structure with proper connections
        structure = base
        for group in groups:
            if random.random() > 0.5:
                # Clean up any potential SMARTS patterns that might have slipped through
                group = ''.join(c for c in group if c not in '[]')
                
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
            # If invalid, try again with a new base
            return self.generate_polymer(domain, trigger)
        
        try:
            # Add hydrogens and clean up
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)  # Add seed for reproducibility
            AllChem.MMFFOptimizeMolecule(mol)
            
            smiles = Chem.MolToSmiles(mol)
            
            # Calculate degradability score
            score = self._calculate_degradability_score(smiles, domain, trigger)
            
            return {
                "smiles": smiles,
                "domain": domain,
                "trigger_condition": trigger or random.choice(template["triggers"]),
                "degradability_score": score,
                "has_kill_switch": True
            }
        except Exception as e:
            # If embedding fails, try again with a new base
            return self.generate_polymer(domain, trigger)

    def _calculate_degradability_score(self, smiles: str, domain: str, trigger: str) -> float:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
            
        # Base score from number of degradable groups
        score = 0.0
        for group_name, group_smarts in self.kill_switch_groups.items():
            pattern = Chem.MolFromSmarts(group_smarts)
            if pattern is not None and mol.HasSubstructMatch(pattern):
                score += 0.2
                
        # Adjust score based on trigger condition
        if trigger:
            if "pH" in trigger:
                score += 0.1
            elif "temperature" in trigger:
                score += 0.15
            elif "UV" in trigger:
                score += 0.2
            elif "microbial" in trigger:
                score += 0.25
                
        return min(1.0, score)

    def generate_dataset(self, num_samples: int = 10000) -> pd.DataFrame:
        data = []
        domains = list(self.domain_templates.keys())
        
        for _ in range(num_samples):
            domain = random.choice(domains)
            trigger = random.choice(self.domain_templates[domain]["triggers"])
            polymer = self.generate_polymer(domain, trigger)
            data.append(polymer)
            
        return pd.DataFrame(data)

def save_dataset(df: pd.DataFrame, output_dir: str = "data"):
    Path(output_dir).mkdir(exist_ok=True)
    
    # Save as CSV
    df.to_csv(f"{output_dir}/polymer_dataset.csv", index=False)
    
    # Save as JSON
    df.to_json(f"{output_dir}/polymer_dataset.json", orient="records", indent=2)
    
    # Save training and test splits
    train_size = int(0.8 * len(df))
    train_df = df[:train_size]
    test_df = df[train_size:]
    
    train_df.to_csv(f"{output_dir}/train.csv", index=False)
    test_df.to_csv(f"{output_dir}/test.csv", index=False)

if __name__ == "__main__":
    generator = PolymerDataGenerator()
    dataset = generator.generate_dataset(num_samples=1000)
    save_dataset(dataset) 