from collections import Counter
import random

def generate_all_kmers(sequence, k):
  """
  Generates all possible k-mers of a specific length from a sequence.

  Args:
      sequence: A DNA sequence string.
      k: The length of the k-mer to generate.

  Returns:
      A set containing all possible k-mers of length k in the sequence.
  """
  return set([sequence[i:i+k] for i in range(len(sequence) - k + 1)])

def calculate_kmer_score(kmer, sequences):
  """
  Calculates a score for a k-mer based on its frequency in the sequences.

  Args:
      kmer: A DNA k-mer string.
      sequences: A list of DNA sequences.

  Returns:
      A float representing the score of the kmer (number of occurrences).
  """
  count = sum(kmer in seq for seq in sequences)
  return count

def find_motif_kmer_scanning_refined(sequences, k, min_score=0):
  """
  Finds the potential motif using k-mer scanning with thresholding and refinement.

  Args:
      sequences: A list of DNA sequences.
      k: The length of the motif to search for.
      min_score: The minimum score threshold for k-mers (optional, defaults to 0).

  Returns:
      A string representing the potential motif and its score (as a tuple).
  """
  kmer_scores = {}
  for sequence in sequences:
    for kmer in generate_all_kmers(sequence, k):
      score = calculate_kmer_score(kmer, sequences)
      if score >= min_score:  # Apply threshold
        kmer_scores[kmer] = score

  if not kmer_scores:
    return None, None  # No k-mers met the threshold

  # Refine the top-scoring k-mer
  best_kmer, best_score = max(kmer_scores.items(), key=lambda item: item[1])
  for _ in range(10):  # Limit refinement iterations (optional)
    new_kmer = best_kmer
    for i in range(k):
      neighbor_kmers = {kmer[:i] + c + kmer[i+1:] for c in 'ACGT' if c != best_kmer[i]}
      neighbor_scores = {kmer: calculate_kmer_score(kmer, sequences) for kmer in neighbor_kmers}
      best_neighbor, neighbor_score = max(neighbor_scores.items(), key=lambda item: item[1])
      if neighbor_score > best_score:
        best_kmer = best_neighbor
        best_score = neighbor_score
        break  # Stop replacing if no improvement

  return best_kmer, best_score

