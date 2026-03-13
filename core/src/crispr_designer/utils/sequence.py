"""
All sequence manipulation based functions
"""
from __future__ import annotations
from typing import Final 


DNA_COMPLEMENT: Final[dict[str,str]] = {
    "A": "T",
    "T": "A",
    "G": "C",
    "C": "G",
    "N": "N",
    # Support lowercase
    "a": "t",
    "t": "a",
    "g": "c",
    "c": "g",
    "n": "n",
}


def complement(sequence: str) -> str:
  """
    This returns the complement of a dna sequence.
    Think of it as a mirror sequence of DNA where each base is replaced
    by its natural partner A with T and C with G according to base pairing rules.
  """
  try:
    return "".join(DNA_COMPLEMENT[nt] for nt in sequence)
  except KeyError as e:
    invalid_char = str(e).strip("'")
    raise ValueError(
      f"Invalid nucleotide '{invalid_char}' in sequence. "
      "Only A, C, G, T, N are allowed."
    ) from None 



def reverse_complement(sequence: str) -> str:
  return complement(sequence)[::-1]

def gc_content(sequence: str) -> float:
  sequence = sequence.upper()
  if not sequence:
    return 0.0
  gc_count = sequence.count("G") + sequence.count("C")
  return gc_count / len(sequence)


def is_valid_dna(sequence: str) -> bool:
  valid_chars = set("ACGTacgt")
  return all(char in valid_chars for char in sequence)



def transcribe(dna_sequence: str) -> str:
  return dna_sequence.upper().replace("T","U")


def sliding_window(sequence: str, window_size: int, step: int = 1) -> list[tuple[int,str]]:
  """Sliding window over a sqeuence"""
  windows: list[tuple[int,str]] =[]
  for i in range(0,len(sequence) - window_size + 1, step):
    windows.append((i,sequence[i:i + window_size]))
  return windows 