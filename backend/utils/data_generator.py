import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import json
from pathlib import Path
import random
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class PolymerDataGenerator:
    def __init__(self):
        self.kill_switch_groups = {
            # Existing groups
            "ester": "C(=O)O",
            "amide": "C(=O)N",
            "azide": "N=N=N",
            "peroxide": "OO",
            "thiol": "S",
            "photolabile": "C(=O)C",
            
            # New additions
            "anhydride": "C(=O)OC(=O)",
            "carbonate": "OC(=O)O",
            "oxalate": "OC(=O)C(=O)O",
            "acetal": "OCOC",
            "disulfide": "SS",
            "thioketal": "SCCS",
            "hydrazone": "C=NNC",
            "imine": "C=N",
            "boronic_ester": "B(O)O",
            "urethane": "NC(=O)O",
            "urea": "NC(=O)N",
            "siloxane": "OSiO"
        }
        
        self.domain_templates = {
            "agriculture": {
                "groups": ["C(=O)O", "C(=O)N", "acetal", "hydrazone", "imine"],
                "triggers": ["pH > 7", "temperature > 30", "microbial_presence"],
                "base_structures": ["CCO", "CC(=O)O", "CC(=O)N"]
            },
            "water": {
                "groups": ["S", "N=N=N", "anhydride", "carbonate", "disulfide"],
                "triggers": ["pH > 8", "UV_light", "oxidation"],
                "base_structures": ["CS", "CN=N=N", "CCO"]
            },
            "urban": {
                "groups": ["C(=O)O", "OO", "thioketal", "boronic_ester", "siloxane"],
                "triggers": ["temperature > 40", "mechanical_stress", "moisture"],
                "base_structures": ["CCO", "COO", "CC(=O)O"]
            }
        }

    def generate_polymer(self, domain: str, trigger: str = None) -> dict:
        template = self.domain_templates[domain]
        base = random.choice(template["base_structures"])
        
        # Add random number of kill switch groups
        num_groups = min(random.randint(1, 3), len(template["groups"]))
        groups = random.sample(template["groups"], num_groups)
        
        # Create polymer structure
        structure = base
        for group in groups:
            if random.random() > 0.5:
                structure += group
        
        # Ensure valid SMILES
        mol = Chem.MolFromSmiles(structure)
        if mol is None:
            return self.generate_polymer(domain, trigger)
        
        # Add hydrogens and clean up
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
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

    def _calculate_degradability_score(self, smiles: str, domain: str, trigger: str) -> float:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return 0.0
            
        # Base score from number of degradable groups
        score = 0.0
        for group in self.domain_templates[domain]["groups"]:
            pattern = Chem.MolFromSmarts(group)
            if mol.HasSubstructMatch(pattern):
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

    def generate_dataset(self, num_samples: int = 1000) -> pd.DataFrame:
        data = []
        domains = list(self.domain_templates.keys())
        
        for _ in range(num_samples):
            domain = random.choice(domains)
            trigger = random.choice(self.domain_templates[domain]["triggers"])
            polymer = self.generate_polymer(domain, trigger)
            data.append(polymer)
            
        return pd.DataFrame(data)

def save_dataset(df: pd.DataFrame, output_dir: str = "data"):
    try:
        Path(output_dir).mkdir(exist_ok=True)
        logging.info(f"Directory '{output_dir}' created or already exists.")
    except Exception as e:
        logging.error(f"Failed to create directory '{output_dir}': {e}")
        return

    try:
        # Save as CSV
        df.to_csv(f"{output_dir}/polymer_dataset.csv", index=False)
        logging.info("CSV file saved successfully.")
    except Exception as e:
        logging.error(f"Failed to save CSV file: {e}")

    try:
        # Save as JSON
        df.to_json(f"{output_dir}/polymer_dataset.json", orient="records", indent=2)
        logging.info("JSON file saved successfully.")
    except Exception as e:
        logging.error(f"Failed to save JSON file: {e}")

    try:
        # Save training and test splits
        train_size = int(0.8 * len(df))
        train_df = df[:train_size]
        test_df = df[train_size:]

        train_df.to_csv(f"{output_dir}/train.csv", index=False)
        test_df.to_csv(f"{output_dir}/test.csv", index=False)
        logging.info("Training and test CSV files saved successfully.")
    except Exception as e:
        logging.error(f"Failed to save training/test CSV files: {e}")

if __name__ == "__main__":
    generator = PolymerDataGenerator()
    dataset = generator.generate_dataset(num_samples=1000)
    save_dataset(dataset)