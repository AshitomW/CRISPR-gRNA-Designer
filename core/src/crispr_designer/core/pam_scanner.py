''' 
PAM Scanner for Pam site detection and guide extraction.
This will scan input sequences for valid PAM sites and extracts spacer sequences according to the CAS System configuration.

'''


from __future__ import annotations
import hashlib
import logging
import re
from dataclasses import dataclass
from typing import Iterator



from crispr_designer.models.cas_systems import CasSystem, get_cas_system
from crispr_designer.config import CasType
from crispr_designer.utils.sequence import reverse_complement


logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class PAMSite:
  """
  A detected pam site with extracted spacer sequence.

  spacer_sequence: guide RNA spacer / protospacer sequence.
  pam_sequence: actual pam sequence found.
  position: 0 indexed start position of the spacer in the input
  pam_position: 0 indexed start position of the pam in the input
  strand: '+' for sense strand and '-' for antisense
  cut_position: predicted cut site position.
  """

  spacer_sequence: str 
  pam_sequence: str
  position: int 
  pam_position: int 
  strand: str 
  cut_position: int 

  @property
  def full_target(self) -> str:
    """ Return the full target sequence including PAM."""
    # For downstream PAM (SpCas9, SaCas9)
    if self.pam_position > self.position:
      return self.spacer_sequence + self.pam_sequence
    # For upstream pam (Cas12a)
    return self.pam_sequence + self.spacer_sequence
  

  def generate_id(self) -> str:
    """ID Generator for this PAM Site"""
    unique_string = f"{self.spacer_sequence}:{self.position}:{self.strand}"
    return hashlib.sha256(unique_string.encode()).hexdigest()[:12]




class PAMScanner: 
  """
  Scans DNA sequence for PAM sites and extracts guide candidates.

  System Specific parameters come from the CasSystem configuration, so this class is basically stateless
  """

  def __init__(self, cas_system: CasSystem | CasType) -> None:
    """
      Initialize scanner with a Cas system configuration
      cas system takes either a CasSystem object or a CasType enumeration.
    """

    if isinstance(cas_system, CasType):
      cas_system = get_cas_system(cas_system)
    self._system = cas_system
    self._pam_regex = cas_system.pam_regex

    logger.debug(
      "Initialized PAM Scanner for %s (PAM: %s, guide length: %d-%d)",
      cas_system.cas_type.value,
      cas_system.pam_pattern,
      cas_system.guide_length_min,
      cas_system.guide_length_max
    )

  @property
  def system(self) -> CasSystem:
    return self._system
  
  def _scan_strand(self, sequence: str, strand: str, guide_length: int) -> Iterator[PAMSite]:
    """
    Scanning a sinel strand for PAM Sites
    Sequence: DNA sequence / already oriented for this strand.
    strand: + or - indicator
    guide_length: length of the spacer to extract

    Will Return Valid PAMSite objects for each valid site found.
    """
    if self._system.pam_location == "downstream":
      # PAM IS 3' of protoscpacer for SpCas9, SaCas9
      yield from self._scan_downstream_pam(sequence, strand, guide_length)
    else:
      # PAM is 5' of protospacer for Cas12a
      yield from self._scan_upstream_pam(sequence, strand, guide_length)



  def _scan_downstream_pam(
      self,
      sequence: str,
      strand: str, 
      guide_length: int, 

  ) -> Iterator[PAMSite]:
    """Handle PAM sites that are downstream 3' of the protospacer."""
    # we need guide length nt before the PAM
    min_start = guide_length
    for match in self._pam_regex.finditer(sequence):
      pam_start = match.start()
      # Checking if there's enough sequence for the spacer.
      if pam_start < guide_length:
        continue 

      spacer_start = pam_start - guide_length
      spacer = sequence[spacer_start:pam_start]
      pam = match.group()

      # calculating cut position , negative offset means upstream of PAM
      cut_position = pam_start + self._system.cut_offset_sense

      yield PAMSite(
        spacer_sequence=spacer, 
        pam_sequence=pam, 
        position=spacer_start, 
        pam_position=pam_start, 
        strand=strand,
        cut_position=cut_position
      )


  def _scan_upstream_pam(
      self,
      sequence: str, 
      strand: str,
      guide_length: int, 
  ):
    """Handle PAM sites that are upstream 5' of the protospacer"""
    for match in self._pam_regex.finditer(sequence):
      pam_start = match.start()
      pam_end = match.end()

      spacer_start = pam_end 
      spacer_end = spacer_start + guide_length
      if spacer_end > len(sequence):
        continue 
      spacer = sequence[spacer_start:spacer_end]
      pam = match.group()

      # For cas12a , cut is downstream of pam
      cut_position = pam_end + self._system.cut_offset_sense

      yield PAMSite(
        spacer_sequence=spacer,
        pam_sequence=pam, 
        position=spacer_start,
        pam_position=pam_start,
        strand=strand,
        cut_position=cut_position
      )


  def find_all_pam_sites(sequence: str, cas_type: CasType,*,include_both_strands: bool=True)-> list[PAMSite]:
    """
    Convenience function to find all PAM sites in a sequence

    sequence: DNA sequence to scan
    cas_type: which cas system to use
    include_both_strands: whether to scan reverse complement.


    Will Return all PAM sites found 
    """

    scanner = PAMScanner(cas_type)
    return scanner.scan(sequence, include_both_strands=include_both_strands)



  def scan(self, sequence: str, * , include_both_strands: bool = True, guide_length: int | None = None) -> list[PAMSite]:

    """
      Scan a pam sequence for PAM sites and extract guide candidates

      sequence: DNA sequence to scan (5'->3', uppercase)
      include_both_strands: if true, this will also can the revertse complement
      guide_length: specific guide length to use: defaults to system minimum.

        In return we get a list of PamSite Objects.

        Value Error , given the guide_length is invalid for the cas system.

    """

    if guide_length is None:
      guide_length = self._system.default_guide_length
    
    if not self._system.validate_guide_length(guide_length):
      raise ValueError(
        f"Guide Length: {guide_length} is invalid for {self._system.cas_type}"
        f"Valid range: {self._system.guide_length_min}-{self._system.guide_length_max}"
      )
  
    sequence = sequence.upper()
    sites: list[PAMSite] = []

    # Scan sense strand 
    sites.extend(self._scan_strand(sequence, "+", guide_length))

    if include_both_strands:
      antisense = reverse_complement(sequence)
      antisense_sites = self._scan_strand(antisense,"-",guide_length)

      for site in antisense_sites:
        full_length = guide_length + self._system.pam_length
        sense_position = len(sequence) - site.position - full_length

        sites.append(PAMSite(
          spacer_sequence=site.spacer_sequence,
          pam_sequence=site.pam_sequence,
          position=sense_position,
          pam_position=len(sequence) - site.pam_position - self._system.pam_length,
          strand="-",
          cut_position=len(sequence) - site.cut_position
        ))

    logger.debug(
      'Found %d PAM sites in %d bp sequence (both strands: %s)',
      len(sites),
      len(sequence),
      include_both_strands
    )

    return sites


    


