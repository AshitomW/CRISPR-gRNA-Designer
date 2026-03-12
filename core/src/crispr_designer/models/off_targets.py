"""
OFF Target analysis data models.
"""

from __future__ import annotations
from typing import Annotated
from pydantic import BaseModel, Field


class MismatchPosition(BaseModel):
  """Details of a single mismatch in an offtarget alignment"""
  position: Annotated[
    int, 
    Field(ge=0,description="0 indexed position in guide PAM Proximal = 0")
  ]

  guide_nucleotide: str 
  target_nucleotide: str 
  is_in_seed_region: bool 
  position_weight: Annotated[
    float,
    Field(
      ge=0.0,
      le=1.0,
      description="Activity retention weight at this position" 
    ),
  ]



class OffTargetHit(BaseModel):
  """ Potential off target site with alignment details"""
  sequence: Annotated[
    str, 
    Field(description="Off target sequence at this site") 
  ]
  chromosome: Annotated[
    str | None , 
    Field(default=None, description="Chromosome/contig name if from FASTA") 
  ]
  position: Annotated[
    int, 
    Field(ge=0, description="0 indexed position in reference") 
  ]
  strand: str 
  pam_sequence: str 
  total_mismatches: Annotated[int, Field(ge=0, description="Total number of mismatches")] 
  seed_region_mismatches: Annotated[int, Field(ge=0,description="Mismatches in PAM proxmimal seed region")]
  mismatch_details: list[MismatchPosition] 
  alignment_score: Annotated[
    float, 
    Field(description="Smith Waterman alignemnt score if gapped ")
  ]
  cfd_score: Annotated[
    float, 
    Field(
      ge=0.0,
      le=1.0,
      description="Cutting Frequency determination score activity prediction." 
    )
  ]
  risk_category: Annotated[
    str, 
    Field(description="HIGH, MEDIUM or LOW Based on mismatch pattern") 
  ]



class OffTargetAnalysisResponse(BaseModel):
  """ Response from off target analyssis endpoint"""
  guide_sequence: str 
  reference_length: int 
  off_targets: list[OffTargetHit]
  total_sites_scanned: int 
  high_risk_count: int 
  medium_risk_count: int 
  low_risk_count: int 
  aggregate_off_target_score: Annotated[
    float, 
    Field(
      ge=0.0,
      le=100.0,
      description="Overall specificity score , 100 => no off targets"
    )
  ]
