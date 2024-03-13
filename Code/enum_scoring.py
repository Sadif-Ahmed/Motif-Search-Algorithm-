from collections import Counter

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

def find_motif_kmer_scanning(sequences, k):
  """
  Finds the potential motif using k-mer scanning based on frequency.

  Args:
      sequences: A list of DNA sequences.
      k: The length of the motif to search for.

  Returns:
      A string representing the potential motif and its score (as a tuple).
  """
  kmer_scores = {}
  for sequence in sequences:
    for kmer in generate_all_kmers(sequence, k):
      kmer_scores[kmer] = calculate_kmer_score(kmer, sequences)

  best_kmer, best_score = max(kmer_scores.items(), key=lambda item: item[1])
  return best_kmer, best_score


