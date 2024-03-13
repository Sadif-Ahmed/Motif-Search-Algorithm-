import time
import numpy as np;
from weighted_greedy import find_greedy_motif_weighted;
from randomized_greedy import find_greedy_motif_randomized;
from pssm import find_greedy_motif_pwm;
from multiple_greedy import find_greedy_motifs;
from background import find_greedy_motif_background;
def read_lines_from_file(filename):
  try:
    with open(filename, 'r') as file:
      lines = file.readlines()  # Read all lines
      return lines
  except FileNotFoundError:
    print(f"Error: File '{filename}' not found.")
    return []  # Return empty list on file not found

def weighted_greedy_motif(sequences,k):
    weights = {"A": 1, "C": 2, "G": 0.5, "T": 3}  # Example weights (modify as needed)
    motif, score = find_greedy_motif_weighted(sequences, k, weights)
    return motif,score

def randomised_greedy_motif(sequences,k):
    motif, score = find_greedy_motif_randomized(sequences, k)
    return motif,score

def multiple_greedy_motif(sequences,k):
    num_motifs = 2
    motifs = find_greedy_motifs(sequences, k, num_motifs)
    return motifs

def pssm_motif(sequences,k):
    motif, score = find_greedy_motif_pwm(sequences, k)
    return motif,score

def background_motif(sequences,k):
    motif, score = find_greedy_motif_background(sequences, k)
    return motif,score

if __name__ == "__main__":
    filename = "data/hm03.txt"

    sequences = read_lines_from_file(filename)
    seq_count=1

    # if sequences:
    #     print("Found Sequences:  "+str(len(sequences)))
    #     for sequence in sequences:
    #         print(str(seq_count)+"  :  "+sequence.strip())  # Remove trailing newline character
    #         seq_count+=1
    # else:
    #     print("No sequences read from the file.")
    
    start_time = time.time()  # Get the current time in seconds

    # Your code to be timed here
    motif,score = weighted_greedy_motif(sequences=sequences,k=20)
    end_time = time.time()
    elapsed_time = end_time - start_time

    print(f"Elapsed time: {elapsed_time:.2f} seconds")
    print("Motif: " + motif)
    print("Score: " + str(score))
