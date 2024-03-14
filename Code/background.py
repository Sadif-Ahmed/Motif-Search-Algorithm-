import numpy as np
from collections import Counter

def calculate_background(sequences):
  background = Counter()
  for sequence in sequences:
    background.update(sequence)
  total_count = sum(background.values())
  return {nucleotide: count / total_count for nucleotide, count in background.items()}

def score_kmer_background(kmer, background):
  score = 0
  for nucleotide in kmer:
    score += np.log2( (1 / background[nucleotide]) )  # Log-likelihood score
  return score

def find_greedy_motif_background(sequences, k, max_iterations=100):
 
  background = calculate_background(sequences)

  best_motif = None
  best_score = float('-inf')

  for _ in range(max_iterations):
    print(_)
    motif_scores = {}
    for sequence in sequences:
      for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        motif_scores[kmer] = score_kmer_background(kmer, background)

    if not motif_scores:  # Check if motif_scores is empty
       break  # Terminate if no k-mers meet the conditions
    # Find the highest scoring k-mer
    best_kmer, score = max(motif_scores.items(), key=lambda item: item[1])
    if score > best_score:
      best_motif = best_kmer
      best_score = score

    # Remove sequences containing the motif
    sequences = [seq for seq in sequences if best_kmer not in seq]

  return best_motif, best_score


