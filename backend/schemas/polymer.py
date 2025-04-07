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

class NomenclatureData(BaseModel):
    name: str = Field(..., description="IUPAC-like name of the polymer")
    formula: str = Field(..., description="Molecular formula")
    functional_groups: List[str] = Field(default_factory=list, description="Identified functional groups")
    structure_type: str = Field(..., description="Structure type (linear, cyclic, etc.)")

class PolymerResponse(BaseModel):
    smiles: str = Field(..., description="SMILES representation of the polymer")
    has_kill_switch: bool = Field(..., description="Whether the polymer has a kill switch")
    degradability_score: float = Field(..., ge=0, le=1, description="Predicted degradability score")
    visualization: Optional[Dict[str, Any]] = Field(None, description="Molecule properties for visualization")
    nomenclature: Optional[NomenclatureData] = Field(None, description="Polymer nomenclature data")

class PredictionRequest(BaseModel):
    smiles: str = Field(..., description="SMILES representation of the molecule to predict") 