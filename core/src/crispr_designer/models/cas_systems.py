'''
CRISPR Cas systems definitions
Each nuclease system is defined with it's PAM sequence, guide length constraints, and cleavage characteristics.
'''

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import ClassVar

from crispr_designer.config import (
  CAS12A_GUIDE_LENGTH_MAX,
  CAS12A_GUIDE_LENGTH_MIN,
  CasType,
  SACAS9_GUIDE_LENGTH_MAX,
  SACAS9_GUIDE_LENGTH_MIN,
  SPCAS9_GUIDE_LENGTH,
)


@dataclass(frozen=True,slots=True)
class CasSystem:
  '''
  Immutable definition of a CRISPR_Cas nuclease system.

  Attributes:
    cas_type: Enum identifier for the nuclease.
    pam_pattern: Regex pattern for PAM recognition (IUPAC codes expanded).
    pam_location: 'downstream' (3' of protospacer) or 'upstream' (5').
    guide_length_min: Minimum spacer sequence length in nucleotides.
    guide_length_max: Maximum spacer sequence length in nucleotides.
    cut_offset_sense: Cleavage position on sense strand relative to PAM.
      Negative = upstream of PAM, positive = downstream
    cut_offset_antisense: Cleavage position on antisense strand
    creates_sticky_ends: Whether the nuclease creates staggered cuts.
  '''
  cas_type: CasType
  pam_pattern: str 
  pam_location: str 
  guide_length_min: int 
  guide_length_max: int
  cut_offset_sense: int 
  cut_offset_antisense: int 
  creates_sticky_ends: bool 


  # IUPAC nucleotide code expansion for regex compilation
  IUPAC_CODES: ClassVar[dict[str,str]] = {
    "N": "[ACGT]",
    "R": "[AG]",
    "Y": "[CT]",
    "S": "[GC]",
    "W": "[AT]",
    "K": "[GT]",
    "M": "[AC]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
  }

  @property
  def pam_regex(self) -> re.Pattern[str]:
    """ Compile PAM Pattern To Regex With IUPAC Expansion"""
    pattern = self.pam_pattern
    for code, expansion in self.IUPAC_CODES.items():
      pattern = pattern.replace(code,expansion)
    return re.compile(pattern)
  
  @property
  def pam_length(self) -> int:
    """ Return the length of the PAM Sequence."""
    return len(self.pam_pattern)
  
  @property
  def default_guide_length(self) -> int:
    """ Return the preferred/defafult guide length for this system. """
    return self.guide_length_min
  
  def validate_guide_length(self,length:int) -> bool:
    """Check if a guide length is valid for this system."""
    return self.guide_length_min <= length <= self.guide_length_max
  





# ================================
# Predefined Cas systems
# ================================


SPCAS9 = CasSystem(
    cas_type=CasType.SPCAS9,
    pam_pattern="NGG",
    pam_location="downstream",
    guide_length_min=SPCAS9_GUIDE_LENGTH,
    guide_length_max=SPCAS9_GUIDE_LENGTH,
    cut_offset_sense=-3, # custs 3 bp upstream of PAM on sense strand
    cut_offset_antisense=-3, # Blunt Cut
    creates_sticky_ends=False
)


SACAS9=CasSystem(
  cas_type=CasType.SACAS9,
  pam_pattern="NNGRRT", # R = A or G
  pam_location="downstream",
  guide_length_min=SACAS9_GUIDE_LENGTH_MIN,
  guide_length_max=SACAS9_GUIDE_LENGTH_MAX,
  cut_offset_sense=-3,
  cut_offset_antisense=-3,
  creates_sticky_ends=False
)


CAS12A = CasSystem(
  cas_type=CasType.CAS12A,
  pam_pattern="TTTV", # V = A, C, or G (not T)
  pam_location="upstream",
  guide_length_min=CAS12A_GUIDE_LENGTH_MIN,
  guide_length_max=CAS12A_GUIDE_LENGTH_MAX,
  cut_offset_sense=18, # Cuts 18 bp downstream of PAM on sense
  cut_offset_antisense=23, # Cuts 23 bp downstream on antisense
  creates_sticky_ends=True, # Creates 5 nt 5' overhangs 
)


# Registry for lookup by type
CAS_SYSTEM_REGISTRY: dict[CasType,CasSystem] = {
  CasType.SPCAS9: SPCAS9,
  CasType.SACAS9: SACAS9,
  CasType.CAS12A: CAS12A,
}




def get_cas_system(cas_type: CasType) -> CasSystem:
  """ Retrieve a CasSystem by its type identifier"""
  if cas_type not in CAS_SYSTEM_REGISTRY:
    registered = ", ".join(ct.value for ct in CAS_SYSTEM_REGISTRY)
    raise KeyError(
      f"Unknown Cas system: {cas_type}. Registered systems: {registered}"
    )
  return CAS_SYSTEM_REGISTRY[cas_type]