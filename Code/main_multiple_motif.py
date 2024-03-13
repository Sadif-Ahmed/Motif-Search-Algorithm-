import csv
import time
import numpy as np;
from weighted_greedy import find_greedy_motif_weighted;
from randomized_greedy import find_greedy_motif_randomized;
from pssm import find_greedy_motif_pwm;
from multiple_greedy import find_greedy_motifs;
from background import find_greedy_motif_background;
from enum_scoring import find_motif_kmer_scanning;
from thresh_refine import find_motif_kmer_scanning_refined;
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
    num_motifs = 5
    motifs = find_greedy_motifs(sequences, k, num_motifs)
    return motifs

def pssm_motif(sequences,k):
    motif, score = find_greedy_motif_pwm(sequences, k)
    return motif,score

def background_motif(sequences,k):
    motif, score = find_greedy_motif_background(sequences, k)
    return motif,score

def enumeration_counting_motif(sequences,k):
    motif, score = find_motif_kmer_scanning(sequences, k)
    return motif,score

def enumeration_counting_refined_motif(sequences,k):
    min_score = 2
    motif, score = find_motif_kmer_scanning_refined(sequences, k, min_score)
    if motif:
        return motif,score
    else:
        print("No motif found above the threshold.")
    

if __name__ == "__main__":
    filename = "data/hm03.txt"

    sequences = read_lines_from_file(filename)
    data = [
  ["Input file", "K", "Motif","Score","Time"],]
    for i in range(8,16):
        K=i
        start_time = time.time()  # Get the current time in seconds
        motifs =multiple_greedy_motif(sequences=sequences,k=K)
        print(motifs)
        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f"Elapsed time: {elapsed_time:.2f} seconds")
        for j in range(len(motifs)):
            data.append(["hm03.txt", K, motifs[j][0].strip(),motifs[j][1],elapsed_time/len(motifs)])

with open('results/multiple_greedy/test1.csv', 'w', newline='') as csvfile:
    # Create a csv writer object
    writer = csv.writer(csvfile)

    # Write the header row
    writer.writerow(data[0])

    # Write the data rows
    for row in data[1:]:
        writer.writerow(row)

    print("CSV file created successfully!")
