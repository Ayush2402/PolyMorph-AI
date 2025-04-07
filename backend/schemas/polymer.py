from pydantic import BaseModel, Field
from typing import Optional, List, Dict, Any
from enum import Enum

class Domain(str, Enum):
    agriculture = "agriculture"
    water = "water"
    urban = "urban"

class PolymerRequest(BaseModel):
    domain: Domain
    trigger_condition: Optional[str] = Field(None, description="Environmental trigger condition")

class PolymerResponse(BaseModel):
    smiles: str = Field(..., description="SMILES representation of the polymer")
    has_kill_switch: bool = Field(..., description="Whether the polymer has a kill switch")
    degradability_score: float = Field(..., ge=0, le=1, description="Predicted degradability score")
    visualization: Optional[Dict[str, Any]] = Field(None, description="Molecule properties for visualization")

class PredictionRequest(BaseModel):
    smiles: str = Field(..., description="SMILES representation of the molecule to predict") 