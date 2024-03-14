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
from background import score_kmer_background;
from background import calculate_background;
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
    num_motifs = 1
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
    
    return motif,score
    

if __name__ == "__main__":
    filename = "data/yst08r.txt"

    sequences = read_lines_from_file(filename)

    for i in range(len(sequences)):
        sequences[i] = sequences[i].strip()

    algo_tools=['Weighted Greedy','Randomised Greedy','Multiple Greedy','Enumeration & Scoring','PSSM','Threshold Refining','Background Model','MEME(cd)','MEME(de)','STREME(cd)','STREME(de)']
    motifs1 = ['AATCAGCACTCTGT','CCCGCCGCCCGCCC','GCGCCGGCCGCGC','CAGCTTAGTGCCTG','AAAAAAAAAAAAAA','AGCTTAGTGCCTGA','GCGCCGGCCGCGC','TTTTCCTCTTCATC','GGTATTATTAAA?C','AAAAAAAAGTGAAA','AAAACAAAGCGAAG']
    motifs2 = ['AAAAAAAAAAAAAA','AGGTGGCGGAGGGG','ATAAAGCTAGAGA','AAGGAAGAAAAAAA','TTTATTTTATTTTT','AAGGAAGAAAAAAA','ATAAAGCTAGAGA','G?ATT?A?ACTGTT','AAAGAAGAAGTGAA','GAAAAGAAAAAAAA','GAAAAGAAAAAAAA']
    motifs3 = ['AAAAAAAAAAAAAA','GGCGCCGCGCCGCC','GTAAGAATAAACG','AAGGAAGAAAAAAA','TTTTTATTTTATTT','AAGGAAGAAAAAAA','GTAAGAATAAACG','AAGGGTT?T?TACT','AATAAACACATACA','GGAAAAAAAAAAAD','AAAAAAAAAAATAG']
    # Your code to be timed here
    data = [
  ["Algorithm/Tool", "K", "Motif","Score"],]
    for i in range(len(algo_tools)):
       motifs3[i]=motifs3[i].replace("?","A")
       motifs3[i]=motifs3[i].replace("D","A")
       score = score_kmer_background(motifs3[i],calculate_background(sequences=sequences))
       data.append([algo_tools[i],14,motifs3[i],score]);    

with open('results/comparison/dataset3.csv', 'w', newline='') as csvfile:
    # Create a csv writer object
    writer = csv.writer(csvfile)

    # Write the header row
    writer.writerow(data[0])

    # Write the data rows
    for row in data[1:]:
        writer.writerow(row)

    print("CSV file created successfully!")
