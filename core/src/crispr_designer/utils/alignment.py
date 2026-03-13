"""Utitlities for Sequence Alignments
Implements Smith Waterman Local Alignment For Off Target Analysis
"""



from __future__ import annotations
import numpy as np
from numpy.typing import NDArray
from crispr_designer.config import SW_MATCH_SCORE, SW_MISMATCH_PENALTY, SW_GAP_OPEN_PENALTY, SW_GAP_EXTEND_PENALTY


def smith_waterman_align(seq_o: str, seq_t: str, *, match_score = SW_MATCH_SCORE, mismatch_penalty= SW_MISMATCH_PENALTY, gap_open= SW_GAP_OPEN_PENALTY, gap_extend = SW_GAP_EXTEND_PENALTY):
  """ 
    Performing smith waterman local alignment
  """
  seq_o = seq_o.upper()
  seq_t = seq_t.upper()

  m, n = len(seq_o), len(seq_t)

  score_matrix: NDArray[np.int32] = np.zeros((m+1,n+1),dtype=np.int32)
  # 0 -> stop, 1-> diagonal, 2->up, 3->left
  traceback: NDArray[np.int8] = np.zeros((m+1,n+1),dtype=np.int8)
  max_score = 0
  max_pos = (0,0)
  for i in range(1, m+1):
    for j in range(1,n+1):
      if seq_o[i-1] == seq_t[j-1]:
        diag_score = score_matrix[i-1,j-1] + match_score
      else:
        diag_score = score_matrix[i-1,j-1] + mismatch_penalty

      # Gap score, this will be a simplified version , not tracking affine gaps fully
      up_score = score_matrix[i-1,j] + gap_extend
      left_score = score_matrix[i,j-1] + gap_extend

    scores = [0, diag_score, up_score, left_score]
    best_score = max(scores)
    score_matrix[i,j] = best_score
    traceback[i,j] = scores.index(best_score)

    if best_score > max_score:
      max_score = best_score
      max_pos = (i,j)
  aligned_o: list[str] = []
  aligned_t: list[str] = []
  i,j = max_pos


  while traceback[i,j] != 0:
    if traceback[i,j] == 1: # Diagonal
      aligned_o.append(seq_o[i-1])
      aligned_t.append(seq_t[j-1])
      i -= 1
      j -= 1
    elif traceback[i,j] == 2: # Up (gap in sequence two)
      aligned_o.append(seq_o[i-1])
      aligned_t.append("-")
      i-=1
    else: # Left (gap in sequence one)
      aligned_o.append("-")
      aligned_t.append(seq_t[j-1])
      j -= 1


  aligned_o_str = "".join(reversed(aligned_o))
  aligned_t_str = "".join(reversed(aligned_t))

  return aligned_o_str, aligned_t_str, int(max_score)




def calculate_alignment_identity(aligned_o: str, aligned_t: str)-> float:
  """
    Calculating percent indentity from aligned sequences.


    seq_o / sequence 1: first aligned sequence (may contain gaps)
    seq_t / sequence 2: second aligned sequence (may contain gaps)

    Will return percent identiy as float in [0.0, 1.0].
  """

  if len(aligned_o) != len(aligned_t):
    raise ValueError("Aligned sequences must have equal lengths")
  
  if not aligned_o:
    return 0.0 

  matches = sum(
    c1 == c2 and c1 != "-" for c1, c2 in zip(aligned_o, aligned_t)
  )

  # Denominator excludes positions where both have gaps
  aligned_positions = sum(not (c1 == "-" and c2 == "-") for c1, c2 in zip(aligned_o, aligned_t))


  if aligned_positions == 0:
    return 0.0

  return matches / aligned_positions

  