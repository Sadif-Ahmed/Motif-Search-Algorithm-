import numpy as np
from collections import Counter

def calculate_background(sequences):
  """
  Calculates the background nucleotide frequencies from the sequences.

  Args:
      sequences: A list of DNA sequences (strings).

  Returns:
      A dictionary containing background nucleotide frequencies (A, C, G, T).
  """
  background = Counter()
  for sequence in sequences:
    background.update(sequence)
  total_count = sum(background.values())
  return {nucleotide: count / total_count for nucleotide, count in background.items()}

def score_kmer_background(kmer, background):
  """
  Calculates the score of a k-mer based on its frequency compared to background.

  Args:
      kmer: A string representing the k-mer sequence.
      background: A dictionary containing background nucleotide frequencies (A, C, G, T).

  Returns:
      A float representing the score of the kmer.
  """
  score = 0
  for nucleotide in kmer:
    score += np.log2( (1 / background[nucleotide]) )  # Log-likelihood score
  return score

def find_greedy_motif_background(sequences, k, max_iterations=100):
  """
  Implements the greedy motif search algorithm with a background model.

  Args:
      sequences: A list of DNA sequences (strings).
      k: The length of the motif to search for.
      max_iterations: The maximum number of iterations for the algorithm (optional).

  Returns:
      A tuple containing the consensus motif (string) and its score (float).
  """
  background = calculate_background(sequences)

  best_motif = None
  best_score = float('-inf')

  for _ in range(max_iterations):
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


