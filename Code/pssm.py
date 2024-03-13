import numpy as np
from collections import Counter
def calculate_pwm(sequences, k):
  """
  Calculates a basic position weight matrix (PWM) from sequences.

  Args:
      sequences: A list of DNA sequences (strings).
      k: The length of the motif to search for.

  Returns:
      A numpy array representing the PWM (shape: (4, k)).
  """
  pwm = np.zeros((4, k))  # 4 for A, C, G, T
  for sequence in sequences:
    for i in range(len(sequence) - k + 1):
      kmer = sequence[i:i+k]
      for j, nucleotide in enumerate(kmer):
        pwm[get_nucleotide_index(nucleotide), j] += 1
  # Normalize to probabilities
  pwm = pwm / (len(sequences) * k)
  return pwm

def get_nucleotide_index(nucleotide):
  """
  Maps a nucleotide character to its corresponding index in the PWM.

  Args:
      nucleotide: A DNA nucleotide character (A, C, G, or T).

  Returns:
      An integer representing the index (0 for A, 1 for C, etc.).
  """
  return {"A": 0, "C": 1, "G": 2, "T": 3}[nucleotide]

def score_kmer_pwm(kmer, pwm):
  """
  Calculates the score of a k-mer based on the PWM.

  Args:
      kmer: A string representing the k-mer sequence.
      pwm: A numpy array representing the position weight matrix.

  Returns:
      A float representing the score of the kmer.
  """
  score = 1.0  # Initialize with 1 to avoid underflow
  for i, nucleotide in enumerate(kmer):
    score *= pwm[get_nucleotide_index(nucleotide), i]
  return np.log2(score)  # Log-likelihood score

def find_greedy_motif_pwm(sequences, k, max_iterations=100):
 """
 Implements the greedy motif search algorithm with a basic PWM.

 Args:
     sequences: A list of DNA sequences (strings).
     k: The length of the motif to search for.
     max_iterations: The maximum number of iterations for the algorithm (optional).

 Returns:
     A tuple containing the consensus motif (string) and its score (float).
 """

 best_motif = None
 best_score = float('-inf')

 for _ in range(max_iterations):
   pwm = calculate_pwm(sequences, k)

   # Initialize motif with the highest scoring k-mer based on the PWM
   motif_scores = {}
   for sequence in sequences:
     for i in range(len(sequence) - k + 1):
       kmer = sequence[i:i+k]
       motif_scores[kmer] = score_kmer_pwm(kmer, pwm)
   initial_motif, _ = max(motif_scores.items(), key=lambda item: item[1])

   # Iteratively refine the motif
   motif = initial_motif
   while True:
     motif_scores = {}
     for sequence in sequences:
       for i in range(len(sequence) - k + 1):
         kmer = sequence[i:i+k]
         motif_scores[kmer] = score_kmer_pwm(kmer, pwm)

     best_kmer, score = max(motif_scores.items(), key=lambda item: item[1])
     if score > best_score:
       best_motif = best_kmer
       best_score = score
       motif = best_kmer  # Update for the next iteration
     else:
       break  # No improvement, stop refining

 return best_motif, best_score


