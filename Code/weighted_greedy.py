from collections import Counter

def score_kmer(kmer, weights):
 
  score = 0
  for i, nucleotide in enumerate(kmer):
    score += weights.get(nucleotide, 0) * (i + 1)  # Weight by position
  return score

def find_greedy_motif_weighted(sequences, k, weights, max_iterations=100):
 
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

