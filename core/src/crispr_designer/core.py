'''
Configuration / Constants
'''
from __future__ import annotations

from enum import Enum
from typing import Final 

# Some Theory for myself to learn before building this, i dont know shit about biology

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


  '''
  SPCAS9_GUIDE_LENGTH: Final[int] = 20
  SACAS9_GUIDE_LENGTH_MIN: Final[int] = 21
  SACAS9_GUIDE_LENGTH_MAX: Final[int] = 23 
  CAS12A_GUIDE_LENGTH_MIN: Final[int] = 23
  CAS12A_GUIDE_LENGTH_MAX: Final[int] = 25



# GC Content Threshold