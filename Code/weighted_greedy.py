from collections import Counter

def score_kmer(kmer, weights):
  """
  Calculates the weighted score of a k-mer based on position weights.

  Args:
      kmer: A string representing the k-mer sequence.
      weights: A dictionary where keys are nucleotide characters (A, C, G, T) 
              and values are corresponding weights for each position in the kmer.

  Returns:
      A float representing the weighted score of the kmer.
  """
  score = 0
  for i, nucleotide in enumerate(kmer):
    score += weights.get(nucleotide, 0) * (i + 1)  # Weight by position
  return score

def find_greedy_motif(sequences, k, weights, max_iterations=100):
  """
  Implements the greedy motif search algorithm with weighted scoring.

  Args:
      sequences: A list of DNA sequences (strings).
      k: The length of the motif to search for.
      weights: A dictionary where keys are nucleotide characters (A, C, G, T) 
              and values are corresponding weights for each position in the kmer.
      max_iterations: The maximum number of iterations for the algorithm (optional).

  Returns:
      A tuple containing the consensus motif (string) and its score (float).
  """
  best_motif = None
  best_score = float('-inf')

  for _ in range(max_iterations):
    motif_counts = Counter()
    for sequence in sequences:
      for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        motif_counts[kmer] += score_kmer(kmer, weights)

    # Find the k-mer with the highest weighted score
    best_kmer, score = motif_counts.most_common(1)[0]
    if score > best_score:
      best_motif = best_kmer
      best_score = score

    # Remove sequences that don't contain the motif
    sequences = [seq for seq in sequences if best_kmer in seq]

  return best_motif, best_score

# Example usage
sequences = ["ACGTGGCT", "TTAGATCC", "ACTGGTCA", "CCAGTTAC"]
k = 5
weights = {"A": 1, "C": 2, "G": 0.5, "T": 3}  # Example weights (modify as needed)

motif, score = find_greedy_motif(sequences, k, weights)

print(f"Motif: {motif}, Score: {score}")