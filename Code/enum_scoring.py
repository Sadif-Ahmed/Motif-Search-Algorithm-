from collections import Counter

def generate_all_kmers(sequence, k):

  return set([sequence[i:i+k] for i in range(len(sequence) - k + 1)])

def calculate_kmer_score(kmer, sequences):
 
  count = sum(kmer in seq for seq in sequences)
  return count

def find_motif_kmer_scanning(sequences, k):
 
  kmer_scores = {}
  for sequence in sequences:
    for kmer in generate_all_kmers(sequence, k):
      kmer_scores[kmer] = calculate_kmer_score(kmer, sequences)

  best_kmer, best_score = max(kmer_scores.items(), key=lambda item: item[1])
  return best_kmer, best_score


