'''
Thermodynamic analysis of guide RNA Sequences.
'''




from __future__ import annotations
import logging
import math
from dataclasses import dataclass
from typing import Sequence


from crispr_designer.config import(
  GC_DISTRIBUTION_MAX_IMBALANCE,
  MIN_HAIRPIN_LOOP_LENGTH,
  MIN_HAIRPIN_STEM_LENGTH,
  NN_DELTA_H,
  NN_DELTA_S,
  NN_INITIATION_DELTA_H,
  NN_INITIATION_DELTA_S,
  STANDARD_NA_CONCENTRATION_MM,
  STANDARD_OLIGO_CONCENTRATION_NM
)

from crispr_designer.models.guides import GCDistribution, ThermodynamicProfile
from crispr_designer.utils.sequence import reverse_complement

logger = logging.getLogger(__name__)


GAS_CONSTANT_R: float = 1.987 # GAS Constant in cal/(mol.K)

@dataclass(frozen=True)
class HairpinStem:
  """Detected hairpin stem structure"""

  stem_start: int 
  stem_end: int 
  loop_start: int 
  loop_end: int 
  stem_length: int 
  sequence: str



  def to_dict(self) -> dict[str, int | str]:
    return {
      "stem_start": self.stem_start,
      "stem_end": self.stem_end,
      "loop_start": self.loop_start,
      "loop_end": self.loop_end,
      "stem_length": self.stem_length,
      "sequence": self.sequence
    }
  


class ThermodynamicAnalyzer:
  """
  Performs thermodynamic analysis on guide RNA sequences.
  Calculations use standard conditions unless overridden.
  """


  def __init__(
      self,
      *,
      na_concentration_mm: float = STANDARD_NA_CONCENTRATION_MM,
      oligo_contentation_nm: float = STANDARD_OLIGO_CONCENTRATION_NM,
  ) -> None:
    self._na_mm = na_concentration_mm
    self._oligo_nm = oligo_contentation_nm

  def calculate_melting_temperature(self, sequence: str) -> float:
    sequence = sequence.upper()
    n = len(sequence)

    if n < 2:
      raise ValueError("Sequence must be at least 2 nucleotides")
    
    delta_h = NN_INITIATION_DELTA_H
    delta_s = NN_INITIATION_DELTA_S

    for i  in range(n - 1):
      dinucleotide = sequence[i:i+2]
      complement = reverse_complement(dinucleotide)
      key = f"{dinucleotide}/{complement}"

      if key in NN_DELTA_H:
        delta_h += NN_DELTA_H[key]
        delta_s += NN_DELTA_S[key]
      else:
        reverse_key = f"{complement}/{dinucleotide}"
        if reverse_key in NN_DELTA_H:
          delta_h += NN_DELTA_H[reverse_key]
          delta_s += NN_DELTA_S[reverse_key]
        else:
          delta_h += -8.0
          delta_s += -22.0 
          logger.warning("Unknown diculeotide pair: %s", key)

    # Convert oligo concentration to molar
    oligo_m = self._oligo_nm * 1e-9
    tm_kelvin = (delta_h * 1000) / (delta_s + GAS_CONSTANT_R * math.log(oligo_m))
    na_m = self._na_mm  * 1e-3
    salt_correction = 16.6 * math.log10(na_m)
    tm_celcius = tm_kelvin - 273.15 + salt_correction

    return tm_celcius
  


  def calculate_gc_content(self, sequence: str) -> float:
    sequence = sequence.upper()
    gc_count = sequence.count("G") + sequence.count("C")
    return gc_count / len(sequence) if sequence else 0.0 
    
  


  def calculate_gc_distribution(self,sequence:str) -> GCDistribution:
    sequence = sequence.upper()
    n = len(sequence)
    third = n//3
    regions=[
      sequence[:third],
      sequence[third:2*third],
      sequence[2*third:],
    ]

    gc_fractions = [self.calculate_gc_content(r) for r in regions]
    max_gc = max(gc_fractions)
    min_gc = min(gc_fractions)
    is_balanced = (max_gc - min_gc) <= GC_DISTRIBUTION_MAX_IMBALANCE


    return GCDistribution(
      five_prime_third=round(gc_fractions[0],3),
      middle_third=round(gc_fractions[1],3),
      three_prime_third=round(gc_fractions[2],3),
      is_balanced=is_balanced 
    )
  



  def detect_hairpins(self,sequence:str) -> list[HairpinStem]:
    """
    Detect potential haripin structures via self complementarity
    A hairpin requires:
    - A stem of at least MIN_HAIRPIN_STEM_LENGTH base pairs
    - A loop of at least MIN_HAIRPIN_LOOP_LENGTH nucleotides
    - Non overlapping stem sequences
    """


    sequence = sequence.upper()
    n = len(sequence)
    hairpins: list[HairpinStem] = []
    min_stem = MIN_HAIRPIN_STEM_LENGTH
    min_loop = MIN_HAIRPIN_LOOP_LENGTH

    for stem_len in range(min_stem, n//2 - min_loop//2 + 1):
      for start in range(n-2*stem_len-min_loop+1):
        stem_5prime = sequence[start:start + stem_len]

        for loop_len in range(min_loop, n-start-2*stem_len+1):
          stem_3prime_start = start+stem_len+loop_len
          stem_3prime_end = stem_3prime_start + stem_len 
          if stem_3prime_end > n:
            break 
          
          stem_3prime = sequence[stem_3prime_start:stem_3prime_end]
          rev_comp = reverse_complement(stem_3prime)
          if stem_5prime == rev_comp:
            hairpin = HairpinStem(
              stem_start=start, 
              stem_end=start+stem_len ,
              loop_start=start+stem_len,
              loop_end=stem_3prime_start,
              stem_length=stem_len, 
              sequence=f"{stem_5prime}...{stem_3prime}",
            )
            hairpins.append(hairpin)



    # Remove redundant (nested) hairpins, keep longest
    if hairpins:
      hairpins.sort(key=lambda h: h.stem_length, reverse=True)
      filtered: list[HairpinStem] = []
      for hp in hairpins:
        is_nested = (existing.stem_start <= hp.stem_start and existing.loop_end >= hp.loop_end for existing in filtered)

        if not is_nested:
          filtered.append(hp)
      
      return filtered
    

    return hairpins


  



  def calculate_seed_region_tm(guide_sequence: str, seed_length: int = 8) -> float:
    analyzer = ThermodynamicAnalyzer()
    # seed is PAM-proximal, which is 3' end for downstream PAM systems
    seed = guide_sequence[-seed_length:]
    return analyzer.calculate_melting_temperature(seed)
  



  




  def analyze(self, guide_sequence: str) -> ThermodynamicProfile:
    guide_sequence = guide_sequence.upper()
    tm = self.calculate_melting_temperature(guide_sequence)
    gc_content = self.calculate_gc_content(guide_sequence)
    gc_dist = self.calculate_gc_distribution(guide_sequence)
    hairpins = self.detect_hairpins(guide_sequence)

    hairpin_details = [h.to_dict() for h in hairpins] if hairpins else None 

    return ThermodynamicProfile(
      melting_temperature_celsius=round(tm,1),
      gc_content=round(gc_content,3),
      gc_distribution=gc_dist,
      has_hairpin=len(hairpins) > 0, 
      hairpin_details=hairpin_details
    )




