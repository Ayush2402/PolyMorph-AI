from fastapi import APIRouter, HTTPException
from typing import List, Optional
import json
from pathlib import Path

from models.polymer_generator import PolymerGenerator
from models.gnn_model import DegradabilityPredictor
from schemas.polymer import PolymerRequest, PolymerResponse, PredictionRequest, NomenclatureData

router = APIRouter()
generator = PolymerGenerator()
predictor = DegradabilityPredictor()

# Initialize history file
history_file = Path("data/history.json")
history_file.parent.mkdir(exist_ok=True)
if not history_file.exists():
    with open(history_file, "w") as f:
        json.dump([], f)

@router.post("/generate", response_model=PolymerResponse)
async def generate_polymer(request: PolymerRequest):
    try:
        # Generate polymer
        smiles = generator.generate_polymer(request.domain, request.trigger_condition)
        
        # Check for kill switch and predict degradability
        has_kill_switch = generator.check_kill_switch(smiles)
        degradability_score = predictor.predict(smiles)
        
        # Get molecule properties for visualization
        molecule_properties = generator.get_molecule_properties(smiles)
        
        # Generate nomenclature data
        nomenclature_data = generator.generate_nomenclature(smiles, request.domain)
        
        # Create response
        response = PolymerResponse(
            smiles=smiles,
            has_kill_switch=has_kill_switch,
            degradability_score=degradability_score,
            visualization=molecule_properties,
            nomenclature=NomenclatureData(
                name=nomenclature_data["name"],
                formula=nomenclature_data["formula"],
                functional_groups=nomenclature_data["functional_groups"],
                structure_type=nomenclature_data["structure_type"]
            )
        )
        
        # Save to history
        with open(history_file, "r+") as f:
            history = json.load(f)
            history.append(response.dict())
            f.seek(0)
            json.dump(history, f)
        
        return response
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.post("/predict")
async def predict_degradability(request: PredictionRequest):
    try:
        score = predictor.predict(request.smiles)
        return {"score": score}
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))

@router.get("/history")
async def get_history():
    try:
        with open(history_file, "r") as f:
            history = json.load(f)
        return history
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e)) 