"""
Edit outcome simulation models 
"""

from __future__ import annotations
from typing import Annotated
from pydantic import BaseModel, Field
from crispr_designer.config import CasType 



class NHEJSimulationRequest(BaseModel):
  """"Request for NHEJ indel simulation"""
  guide_sequence: Annotated[str, Field(min_length=17, max_length=30)]
  target_sequence: Annotated[
    str, 
    Field(min_length=50, description="Sequence context around cut site") 
  ]
  cas_system: CasType 
  simulation_count: Annotated[ 
    int, 
    Field(default=100, ge=10, le=1000)
  ]
  random_seed: Annotated[int | None, Field(default=None, description="Seed for reproducable outcomes")]


class IndelOutcome(BaseModel):
  """ A single simulated indel outcome. """

  indel_type: Annotated[str, Field(description="'Insertion' Or 'Deletion'")]
  indel_size: Annotated[int, Field(description="Size in bp (positive mean insertion, negative means deletion)")]
  edited_sequence: str 
  frameshift: Annotated[
    bool, 
    Field(description="Whether this indel causes a frameshift"),
  ]
  microhomology_mediated: Annotated[ 
    bool,
    Field(description="Whether deletion shows microhomology signature")
  ]
  frequency: Annotated[
    float, 
    Field(ge=0.0, le=1.0, description="Fraction of simulations with this outcome")
  ]


class NHEJSimulationResponse(BaseModel):
  """ Response from NHEJ simulation endpoint """
  original_sequence: str
  cut_position: int 
  outcomes: list[IndelOutcome]
  frameshift_probability: Annotated[
    float,
    Field(
      ge=0.0,
      le=1.0,
      description="Fraction of outcomes causing frameshift" 
    )
  ]
  most_common_outcome: IndelOutcome 
  indel_size_distribution: dict[int, float] 



class HDRSimulationRequest(BaseModel):
  """ Request For HDR outcome prediction"""
  guide_sequence: Annotated[str, Field(min_length=17,max_length=30)]
  target_sequence: Annotated[str, Field(min_length=50)]

  donor_template: Annotated[str, Field(min_length=10, description="HDR donor template sequence")]
  cas_system: CasType 
  left_homology_arm_length: Annotated[
    int, 
    Field(default=50, ge=20, le=500)
  ]
  right_homology_arm_length: Annotated[
    int,
    Field(default=50,ge=20, le=500)
  ]




class HDRSimulationResponse(BaseModel):
  """Response from HDR Simulation Endpoint"""
  original_sequence: str 
  edited_sequence: str 
  donor_template: str 
  cut_position: int 
  insertion_start: int 
  insertion_end: int 
  left_homology_arm: str 
  right_homology_arm: str 
  edit_description: str 


  
