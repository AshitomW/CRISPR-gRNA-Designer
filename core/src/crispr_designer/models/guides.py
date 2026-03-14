"""
Pydantic models for guide RNA data structures.
All API request/reponse schemas are defined here with validations.
"""


from __future__ import annotations
import re 
from typing import Annotated
from pydantic import BaseModel, ConfigDict, Field, field_validator
from crispr_designer.config import CasType, MAX_SEQUENCE_LENGTH



DNA_PATTERN = re.compile(r"[ACGTacgt]+$")



def validate_dna_sequence(sequence: str) -> str:
  sequence = sequence.upper().strip()
  if not sequence:
    raise ValueError("Sequence cannot be empty.")
  if not DNA_PATTERN.match(sequence):
    invalid_charts = set(sequence) - set("ACGT")
    raise ValueError(
      f"Invalid nucleotides in sequence: {invalid_charts}. "
      "Only A,C,G,T are allowed."
    )
  return sequence





class GuideDesignRequest(BaseModel):
  model_config = ConfigDict(str_strip_whitespace=True)
  sequence: Annotated[
    str,
    Field(
      min_length=50,
      max_length=MAX_SEQUENCE_LENGTH,
      description="Target DNA sequence (5'->3',sense strand)",
      examples=["ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"]
    )
  ]
  cas_system: Annotated[
    CasType,
    Field(
      description="CRISPR-Cas system to use for guide design",
      examples=[CasType.SPCAS9], 
    )
  ]
  include_both_strands: Annotated[
    bool,
    Field(
      default=True,
      description="Search for PAM sites on both strands",   
    )
  ]
  max_results: Annotated[ 
    int,
    Field(
      default=50,
      ge=1,
      le=500,
      description="Maximum number of guides to return",
    )
  ]


  @field_validator("sequence")
  @classmethod
  def validate_sequence(cls, value:str) -> str:
    return validate_dna_sequence(value) 
  



class OffTargetRequest(BaseModel):
  model_config = ConfigDict(str_strip_whitespace=True)
  guide_sequence: Annotated[
    str,
    Field(
      min_length=17,
      max_length=30,
      description="Guide RNS Spacer sequence (5'->3', DNA Representation)"
    ),
  ]


  reference_sequence: Annotated[
    str, 
    Field(
      min_length=100,
      description="Reference genome or sequence to scan (FASTA or raw)" 
    )
  ]


  cas_system: CasType
  
  max_mismatches: Annotated[
    int,
    Field(
      default=4,
      ge=0,
      le=6,
      description="Maxmimum mismatches to report as potential off target"
    )
  ]

  max_results: Annotated[
    int,
    Field(default=100, ge=1, le=1000),
  ]


  @field_validator("guide_sequence")
  @classmethod
  def validate_guide(cls, value: str) -> str:
    return validate_dna_sequence(value)
  

# Some Response Models


class GCDistribution(BaseModel):
  ''' GC Content Distribution across guide regions'''
  five_prime_third: Annotated[float, Field(ge=0.0,le=1.0,description="GC fraction in 5' third of guide")]

  middle_third: Annotated[
    float,
    Field(ge=0.0, le=1.0, description="GC Fraction in middle third")
  ]

  three_prime_third: Annotated[
    float,
    Field(ge=0.0,le=1.0,description="GC fraction in 3' third")
  ]

  is_balanced: Annotated[bool,Field(description="Whether distribution is within +- 15% balance")]




class ThermodynamicProfile(BaseModel):
  ''' Thermodynamic analysis result for a guide'''
  melting_temperature_celsius: Annotated[
    float,
    Field(description="Predicted Tm using nearest neighbor model")
  ]

  gc_content: Annotated[
    float,
    Field(ge=0.0, le=1.0, description="Overall GC fraction")
  ]

  gc_distribution: GCDistribution

  has_hairpin: Annotated[
    bool, 
    Field(description="Whether self complementary structure detected")
  ]

  hairpin_details: Annotated[
    list[dict[str, int | str]] | None, 
    Field(default=None, description="Stem positions if hairpin detected")
  ]

class EfficiencyScores(BaseModel):
  ''' Breakdown of efficiency scoring components'''
  gc_score: Annotated[
    float,
    Field(description="Score from GC content from 0 to 1")
  ]
  seed_stability_score: Annotated[
    float,
    Field(description="Score from seed region stability from 0 to 1")
  ]
  off_target_score: Annotated[ 
    float,
    Field(description="Penalty from predicted off targets , 0 to 1 , lower value mens worse")
  ]

  hairpin_score: Annotated[
    float, 
    Field(description="Penalty from self complementarily , 0-1 lower is worse") 
  ]

  composite_score: Annotated[ 
    float,
    Field(description="Weighted composite score (0-100)") 
  ]


class GuideCandidate(BaseModel):
  """ Single guide RNA candidate with full analysis"""
  id: Annotated[
    str,
    Field(description="Unique identifier for this guide") 
  ]
  spacer_sequence: Annotated[str, Field(description="20-25nt spacer sequence (DNA,5'->3')")]

  pam_sequence: Annotated[str, Field(description="PAM Sequence as foudn in target")]

  full_target: Annotated[
    str,
    Field(description="Full target including PAM"),
  ]

  strand: Annotated[
    str,
    Field(description="'+' for sense, '-' for antisense")
  ]

  position: Annotated[
    int, 
    Field(description="Predicted cut site position") 
  ]

  cas_system: CasType
  thermodynamics: ThermodynamicProfile
  efficiency: EfficiencyScores

  rank: Annotated[int, Field(ge=1, description="Rank among all candidates, 1 means it's the best")]






class GuideDesignResponse(BaseModel):
  ''' Response from guide design endpoint.'''
  guides: list[GuideCandidate]
  total_pam_sites_found: Annotated[
    int, 
    Field(description="Total PAM sites detected before filtering")
  ]
  input_sequence_length: int 
  cas_system_used: CasType
  analysis_parameters: dict[str,int|float|bool|str]