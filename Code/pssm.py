import numpy as np
from collections import Counter
def calculate_pwm(sequences, k):
 
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

  return {"A": 0, "C": 1, "G": 2, "T": 3}[nucleotide]

def score_kmer_pwm(kmer, pwm):
 
  score = 1.0  # Initialize with 1 to avoid underflow
  for i, nucleotide in enumerate(kmer):
    score *= pwm[get_nucleotide_index(nucleotide), i]
  return np.log2(score)  # Log-likelihood score

def find_greedy_motif_pwm(sequences, k, max_iterations=100):


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


