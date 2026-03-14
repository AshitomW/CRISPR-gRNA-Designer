"""
Off-target site detection and scoring.

"""
from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field

from crispr_designer.config import (
    MAX_MISMATCHES_FOR_OFF_TARGET,
    MISMATCH_POSITION_WEIGHTS,
    SEED_REGION_LENGTH,
    SW_GAP_EXTEND_PENALTY,
    SW_GAP_OPEN_PENALTY,
    SW_MATCH_SCORE,
    SW_MISMATCH_PENALTY,
)
from crispr_designer.models.cas_systems import CasSystem
from crispr_designer.models.off_targets import MismatchPosition, OffTargetHit
from crispr_designer.utils.alignment import smith_waterman_align
from crispr_designer.utils.sequence import reverse_complement

logger = logging.getLogger(__name__)


def compute_hamming_distance(seq1: str, seq2: str) -> int:

    if len(seq1) != len(seq2):
        raise ValueError(
            f"Sequences must have equal length for Hamming distance. "
            f"Got {len(seq1)} and {len(seq2)}."
        )
    
    return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))


@dataclass
class OffTargetCandidate:
  
    
    sequence: str
    position: int
    strand: str
    pam: str
    chromosome: str | None = None
    mismatches: list[tuple[int, str, str]] = field(default_factory=list)
    
    @property
    def mismatch_count(self) -> int:
       
        return len(self.mismatches)
    
    @property
    def seed_mismatches(self) -> int:
        """Return number of mismatches in seed region."""
        return sum(1 for pos, _, _ in self.mismatches if pos < SEED_REGION_LENGTH)


class OffTargetAnalyzer:

    
    def __init__(
        self,
        cas_system: CasSystem,
        *,
        max_mismatches: int = MAX_MISMATCHES_FOR_OFF_TARGET,
    ) -> None:

        self._system = cas_system
        self._max_mismatches = max_mismatches
        self._position_weights = self._get_position_weights()
    
    def _get_position_weights(self) -> tuple[float, ...]:

        guide_len = self._system.default_guide_length
        base_weights = MISMATCH_POSITION_WEIGHTS
        
        if guide_len == len(base_weights):
            return base_weights
        
        if guide_len < len(base_weights):
            return base_weights[:guide_len]
        
        # Extend with the last weight for longer guides
        extension = (base_weights[-1],) * (guide_len - len(base_weights))
        return base_weights + extension
    
    def analyze(
        self,
        guide_sequence: str,
        reference_sequence: str,
        *,
        chromosome: str | None = None,
        max_results: int = 100,
    ) -> list[OffTargetHit]:

        guide_sequence = guide_sequence.upper()
        reference_sequence = reference_sequence.upper()
        guide_len = len(guide_sequence)
        
        candidates: list[OffTargetCandidate] = []
        
        # Scan forward strand
        forward_candidates = self._scan_strand(
            guide=guide_sequence,
            reference=reference_sequence,
            strand="+",
            chromosome=chromosome,
        )
        candidates.extend(forward_candidates)
        
        # Scan reverse strand
        rev_ref = reverse_complement(reference_sequence)
        reverse_candidates = self._scan_strand(
            guide=guide_sequence,
            reference=rev_ref,
            strand="-",
            chromosome=chromosome,
        )
        
        # Convert reverse strand positions to forward strand coordinates
        for candidate in reverse_candidates:
            candidate.position = len(reference_sequence) - candidate.position - guide_len
        
        candidates.extend(reverse_candidates)
        
        # Convert to OffTargetHit models and score
        hits = [self._candidate_to_hit(c, guide_sequence) for c in candidates]
        
        # Sort by CFD score descending (highest activity risk first)
        hits.sort(key=lambda h: h.cfd_score, reverse=True)
        
        logger.info(
            "Found %d potential off-targets (max mismatches: %d)",
            len(hits),
            self._max_mismatches,
        )
        
        return hits[:max_results]
    
    def _scan_strand(
        self,
        guide: str,
        reference: str,
        strand: str,
        chromosome: str | None,
    ) -> list[OffTargetCandidate]:

        guide_len = len(guide)
        pam_regex = self._system.pam_regex
        
        if self._system.pam_location == "downstream":
            return self._scan_downstream_pam_sites(
                guide=guide,
                reference=reference,
                strand=strand,
                chromosome=chromosome,
                guide_len=guide_len,
                pam_regex=pam_regex,
            )
        
        return self._scan_upstream_pam_sites(
            guide=guide,
            reference=reference,
            strand=strand,
            chromosome=chromosome,
            guide_len=guide_len,
            pam_regex=pam_regex,
        )
    
    def _scan_downstream_pam_sites(
        self,
        guide: str,
        reference: str,
        strand: str,
        chromosome: str | None,
        guide_len: int,
        pam_regex: re.Pattern[str],
    ) -> list[OffTargetCandidate]:

        candidates: list[OffTargetCandidate] = []
        
        for match in pam_regex.finditer(reference):
            pam_start = match.start()
            
            # Check if there's enough sequence upstream for the target
            if pam_start < guide_len:
                continue
            
            target_start = pam_start - guide_len
            target = reference[target_start:pam_start]
            
            # Fast Hamming distance filter
            hamming_dist = compute_hamming_distance(guide, target)
            if hamming_dist > self._max_mismatches:
                continue
            
            # Passed filter — do detailed mismatch analysis
            mismatches = self._find_mismatches(guide, target)
            
            candidates.append(OffTargetCandidate(
                sequence=target,
                position=target_start,
                strand=strand,
                pam=match.group(),
                chromosome=chromosome,
                mismatches=mismatches,
            ))
        
        return candidates
    
    def _scan_upstream_pam_sites(
        self,
        guide: str,
        reference: str,
        strand: str,
        chromosome: str | None,
        guide_len: int,
        pam_regex: re.Pattern[str],
    ) -> list[OffTargetCandidate]:

        candidates: list[OffTargetCandidate] = []
        
        for match in pam_regex.finditer(reference):
            pam_end = match.end()
            target_end = pam_end + guide_len
            
            # Check if there's enough sequence downstream for the target
            if target_end > len(reference):
                continue
            
            target = reference[pam_end:target_end]
            
            # Fast Hamming distance filter
            hamming_dist = compute_hamming_distance(guide, target)
            if hamming_dist > self._max_mismatches:
                continue
            
            # Passed filter — do detailed mismatch analysis
            mismatches = self._find_mismatches(guide, target)
            
            candidates.append(OffTargetCandidate(
                sequence=target,
                position=pam_end,
                strand=strand,
                pam=match.group(),
                chromosome=chromosome,
                mismatches=mismatches,
            ))
        
        return candidates
    
    def _find_mismatches(
        self,
        guide: str,
        target: str,
    ) -> list[tuple[int, str, str]]:

        mismatches: list[tuple[int, str, str]] = []
        
        if self._system.pam_location == "downstream":
            # For downstream PAM, position 0 is the 3' end (PAM-proximal)
            for i in range(len(guide)):
                if guide[i] != target[i]:
                    pam_proximal_pos = len(guide) - 1 - i
                    mismatches.append((pam_proximal_pos, guide[i], target[i]))
        else:
            # For upstream PAM, position 0 is at the 5' end (PAM-proximal)
            for i in range(len(guide)):
                if guide[i] != target[i]:
                    mismatches.append((i, guide[i], target[i]))
        
        return mismatches
    
    def _candidate_to_hit(
        self,
        candidate: OffTargetCandidate,
        guide: str,
    ) -> OffTargetHit:
        """Convert internal candidate to API response model."""
        mismatch_details: list[MismatchPosition] = []
        
        for pos, guide_nt, target_nt in candidate.mismatches:
            weight = (
                self._position_weights[pos]
                if pos < len(self._position_weights)
                else 1.0
            )
            mismatch_details.append(MismatchPosition(
                position=pos,
                guide_nucleotide=guide_nt,
                target_nucleotide=target_nt,
                is_in_seed_region=pos < SEED_REGION_LENGTH,
                position_weight=weight,
            ))
        
        cfd_score = self._calculate_cfd_score(candidate.mismatches)
        risk = self._categorize_risk(
            candidate.mismatch_count,
            candidate.seed_mismatches,
            cfd_score,
        )
        
        return OffTargetHit(
            sequence=candidate.sequence,
            chromosome=candidate.chromosome,
            position=candidate.position,
            strand=candidate.strand,
            pam_sequence=candidate.pam,
            total_mismatches=candidate.mismatch_count,
            seed_region_mismatches=candidate.seed_mismatches,
            mismatch_details=mismatch_details,
            alignment_score=self._calculate_alignment_score(guide, candidate.sequence),
            cfd_score=cfd_score,
            risk_category=risk,
        )
    
    def _calculate_cfd_score(
        self,
        mismatches: list[tuple[int, str, str]],
    ) -> float:

        if not mismatches:
            return 1.0
        
        score = 1.0
        for pos, _, _ in mismatches:
            if pos < len(self._position_weights):
                score *= self._position_weights[pos]
        
        return round(score, 4)
    
    def _calculate_alignment_score(self, guide: str, target: str) -> float:

        _, _, score = smith_waterman_align(
            guide,
            target,
            match_score=SW_MATCH_SCORE,
            mismatch_penalty=SW_MISMATCH_PENALTY,
            gap_open=SW_GAP_OPEN_PENALTY,
            gap_extend=SW_GAP_EXTEND_PENALTY,
        )
        return float(score)
    
    def _categorize_risk(
        self,
        total_mismatches: int,
        seed_mismatches: int,
        cfd_score: float,
    ) -> str:

        # High risk: few mismatches, especially if none in seed
        if total_mismatches <= 1 and seed_mismatches == 0:
            return "HIGH"
        if cfd_score >= 0.5:
            return "HIGH"
        
        # Medium risk: moderate mismatches or some in seed
        if total_mismatches <= 2 or cfd_score >= 0.1:
            return "MEDIUM"
        
        return "LOW"