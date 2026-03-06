'''
Configuration / Constants
'''
from __future__ import annotations

from enum import Enum
from typing import Final 

# Some Theory for myself to learn while building this, i dont know shit about biology

'''
Some Theory:

CRISPR:
  - CRISPR : Clustered Regularly Interspaced Short Palindromic Repeats?
  - In bacteria, CRISPR sequences in genome stores snippets of DNA from past viral infections
  - These snippets act as a memory , that allows bacteria to recognize and defend against the same viruses in the future.


CAS Proteins:

- Cas : CRISPR-associated proteins
- Many cas proteins function as nucleases, meaning they can cut DNA OR RNA.
- THe most commonly used Cas protien in genome editing is Cas9, which cuts double stranded DNA at specification locations



Guide RNA:

- Short RNA Sequence that matches the target DNA in the genome.
- It forms a complex with a Cas nuclease, guiding it to the precise location to cut.

\
'''


class CasType(str, Enum):
  ''' Allowed CRISPR-Cas nuclease systems.'''
  SPCAS9 = "SpCas9"
  SACAS9 = "SaCas9"
  CAS12A = "Cas12a"

  # Guide RNA Length Constratins per system
  # From Some of the papers about SpCas9, SaCas9, Cas12A

'''
    Guide RNA length constraints per system is specific length requirement of the guide RNA (gRNA) for different CRISPR nuclease systems.

    The guide RNA is the part that directs Cas proteins to correct the sequence of dna rna . so length matters for efficiency and specificity.

    The length of grna affects , target specificity: too short means it can bind unintended sequences off-target effects.

    Cutting efficiency: Too long means it may form secondary structures , reducing efficiency. 

'''
SPCAS9_GUIDE_LENGTH: Final[int] = 20
SACAS9_GUIDE_LENGTH_MIN: Final[int] = 21
SACAS9_GUIDE_LENGTH_MAX: Final[int] = 23 
CAS12A_GUIDE_LENGTH_MIN: Final[int] = 23
CAS12A_GUIDE_LENGTH_MAX: Final[int] = 25



'''
  GC : Percentage of nucleotides in a DNA sequence that are either Guanine(G) or Cytosine(C).

  DNA has four bases:

  - A : Adenine
  - T : Thymine
  - G : Guanine
  - C : Cytosine

  So GC Content  measures how much of the sequence is G OR C instead of A or T  

  For CRISPR guides like SpCas9, GC content affects how well the guide binds to the target DNA.   

'''

# ============
# GC Content Thresholds
# ============
# These Optimal GC content range come from various literatures. 

GC_CONTENT_MIN_OPTIMAL: Final[float] = 0.40
GC_CONTENT_MAX_OPTIMAL: Final[float] = 0.70
# Extreme GC Content threshold that severly impact efficiency
GC_CONTENT_MIN_ACCEPTABLE: Final[float] = 0.30
GC_CONTENT_MAX_ACCEPTABLE: Final[float] = 0.80
# Maximum allowed imbalance between GC distribution in thirds of guide
GC_DISTRIBUTION_MAX_IMBALANCE: Final[float] = 0.15



# ===================
# EFFICIENCY SCORING WEIGHTS
# ===================
# These weights come from multiple studies and independent validations.  
# Let the literature do the heavy lifting while I just code.


# GC Balance contributes about 25% to overall efficiency predction.
GC_BALANCE_WEIGHT: Final[float] = 0.25
# Seed Region (PAM-proximal 8nt) stability is the strongest predictor. 
SEED_STABILITY_WEIGHT: Final[float] = 0.30
# Off target potential is a negative contributor
OFF_TARGET_WEIGHT: Final[float] = -0.35
# Self-Complementarity reduces guide availability
HAIRPIN_WEIGHT: Final[float] = -0.10



'''
  PAM : Protospacer Adjacent Motif
  It is a short DNA sequence located next to the target site that a CRISPR nuclease must recognize before it can cut the DNA.

  Kind of like a required permission tag for the cas protein.

  
  CRISPR Targeting needs two things:

  1) Guide RNA match -> Sequence complementary to the target DNA.
  2) Correct PAM nearby -> Short motif the cas protein recognizes

  If the pam is missing, the cas enzymes will not cut even if the guide rna perfectly matches the dna

'''


# ==============================
# OFF-Target Analysis Parameters
# ==============================

# Positions-weighted mismatch penalties for 20 nt guide
# Positions 1-8 (seed region, PAM proximal) are most critical
# Index 0 = position 1 (PAM-proximal) , Index 19 = poition 20 (PAM-distal)
MISMATCH_POSITION_WEIGHTS: Final[tuple[float,...]] = (
  # Seed Region (Positions 1-8: High penalty)
  0.0, # Position 1: Mismatch here abolishes activity
  0.0, # Position 2 
  0.1,
  0.1,
  0.2,
  0.3,
  0.4,
  0.5, # Position 8
  # Transition Region (positions 9 - 12)
  0.6,
  0.7,
  0.7,
  0.8,
  # PAM-distal region (positions 13-20): more tolerant
  0.85,
  0.85,
  0.9,
  0.9,
  0.95,
  0.95,
  1.0,
  1.0
)

# Seed region spans PAM-proximal 8 nucleotides
SEED_REGION_LENGTH: Final[int] = 8
# Maximum mismatches to consider a potential off-target
MAX_MISMATCHES_FOR_OFF_TARGET: Final[int] = 4

# Smith-Waterman alignment parameters
SW_MATCH_SCORE: Final[int] = 2
SW_MISMATCH_PENALTY: Final[int] = -1
SW_GAP_OPEN_PENALTY: Final[int] = -3
SW_GAP_EXTEND_PENALTY: Final[int] = -1



# ========================
# Thermodynamic Parameters
# ========================

# Nearest neighbor thermodynamic parameters
# Delta(H) in kcal/mol, Delta(S) in cal/(mol.K)
# Keys are dinucleotide steps (5'->3'/3'->5')
NN_DELTA_H: Final[dict[str,float]] = {
    "AA/TT": -7.9,
    "AT/TA": -7.2,
    "TA/AT": -7.2,
    "CA/GT": -8.5,
    "GT/CA": -8.4,
    "CT/GA": -7.8,
    "GA/CT": -8.2,
    "CG/GC": -10.6,
    "GC/CG": -9.8,
    "GG/CC": -8.0,
}

NN_DELTA_S: Final[dict[str, float]] = {
    "AA/TT": -22.2,
    "AT/TA": -20.4,
    "TA/AT": -21.3,
    "CA/GT": -22.7,
    "GT/CA": -22.4,
    "CT/GA": -21.0,
    "GA/CT": -22.2,
    "CG/GC": -27.2,
    "GC/CG": -24.4,
    "GG/CC": -19.9,
}


# Initiation Parameters for terminal base pairs
NN_INITIATION_DELTA_H: Final[float] = 0.1
NN_INITIATION_DELTA_S: Final[float] = -2.8

# Standard conditions for Tm calculations
STANDARD_NA_CONCENTRATION_MM: Final[float] = 50.0 # mM Na+
STANDARD_OLIGO_CONCENTRATION_NM: Final[float] = 250.0 # nM oligonucleotide



# Minimum stem length for hairpin detection
MIN_HAIRPIN_STEM_LENGTH: Final[int]  = 4
# Minimum loop length for stalbe hairpin
MIN_HAIRPIN_LOOP_LENGTH: Final[int] = 3



#==============================
# EDITING SIMULATION PARAMETERS
#==============================



# Number of NHEJ outcomes to simulate
# provides statistically meaningful distribution of repair outcomes
NHEJ_SIMULATION_COUNT: Final[int] = 100


# Indel size probability distribution (favors small indels)
# Maps indel size to relative probability weight
INDEL_SIZE_PROBABILITIES: Final[dict[int,float]] = {
  1: 0.35, # +1/-1 bp most common
  2: 0.25,
  3: 0.15,
  4: 0.10,
  5: 0.05,
  6: 0.04,
  7: 0.03,
  8: 0.02,
  9: 0.005,
  10: 0.005,
}


# Insertion vs deletion probability
# Deletions are more common than insertions in NHEJ repair
DELETION_PROBABILITY: Final[float] = 0.65
INSERTION_PROBABILITY: Final[float] = 0.35

# Microhomology-mediated deletion probability
MICROHOMOLOGY_SEARCH_WINDOW: Final[int] = 20



# ===============================
# API CONFIGURATION
# ===============================


API_PREFIX: Final[str] = "/api/v1"
DEFAULT_PAGE_SIZE: Final[int] = 50
MAX_PAGE_SIZE: Final[int] = 200
MAX_SEQUENCE_LENGTH: Final[int] = 100_000 # 100 kb limit per request
MAX_REFERENCE_LENGTH: Final[int] = 10_000_000 # 10 Mb reference genome

