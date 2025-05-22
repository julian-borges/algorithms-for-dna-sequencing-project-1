# Bioinformatics - Genomic Data Science
# Algorithms for DNA Sequencing
# The Johns Hopkins University
# Student: Julian Borges
# Programming Homework 1

import matplotlib.pyplot as plt
from collections import Counter

# -----------------------------
# Utility: Genome Reader
# -----------------------------
def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

# -----------------------------
# Utility: Reverse Complement
# -----------------------------
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

# -----------------------------
# Utility: Naive Exact Match
# -----------------------------
def naive_match(pattern, text):
    positions = []
    for i in range(len(text) - len(pattern) + 1):
        match = True
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                match = False
                break
        if match:
            positions.append(i)
    return positions

# -----------------------------
# Utility: Approximate Match ≤ 2 mismatches
# -----------------------------
def naive_2mm(pattern, text):
    positions = []
    for i in range(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            positions.append(i)
    return positions

# -----------------------------
# Utility: Count Exact Matches
# -----------------------------
def count_occurrences(pattern, sequence):
    count = 0
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i + len(pattern)] == pattern:
            count += 1
    return count

# -----------------------------
# Question 1 - Count AGGT and ACCT
# -----------------------------
genome_path = '/Users/FxMED/PycharmProjects/Introduction to Python/Introduction/Comments/lambda_virus.fa'
genome = read_genome(genome_path)
pattern1 = 'AGGT'
revcomp1 = reverse_complement(pattern1)
print("Q1:")
print(f"Forward '{pattern1}': {count_occurrences(pattern1, genome)}")
print(f"Reverse '{revcomp1}': {count_occurrences(revcomp1, genome)}\n")

# -----------------------------
# Question 2 - Count TTAA and reverse
# -----------------------------
pattern2 = 'TTAA'
revcomp2 = reverse_complement(pattern2)
print("Q2:")
print(f"Forward '{pattern2}': {count_occurrences(pattern2, genome)}")
print(f"Reverse '{revcomp2}': {count_occurrences(revcomp2, genome)}\n")

# -----------------------------
# Question 3 - Leftmost match for ACTAAGT or reverse
# -----------------------------
pattern3 = 'ACTAAGT'
rev3 = reverse_complement(pattern3)
positions3 = naive_match(pattern3, genome)
rev_positions3 = naive_match(rev3, genome)
min3 = min(
    positions3[0] if positions3 else float('inf'),
    rev_positions3[0] if rev_positions3 else float('inf')
)
print("Q3:")
print(f"Leftmost occurrence: {min3}\n")

# -----------------------------
# Question 4 - Leftmost match for AGTCGA or reverse
# -----------------------------
pattern4 = 'AGTCGA'
rev4 = reverse_complement(pattern4)
pos4 = naive_match(pattern4, genome)
rev_pos4 = naive_match(rev4, genome)
min4 = min(
    pos4[0] if pos4 else float('inf'),
    rev_pos4[0] if rev_pos4 else float('inf')
)
print("Q4:")
print(f"Leftmost occurrence: {min4}\n")

# -----------------------------
# Question 5 - Approximate matches ≤2 for TTCAAGCC
# -----------------------------
pattern5 = 'TTCAAGCC'
matches5 = naive_2mm(pattern5, genome)
print("Q5:")
print(f"Matches with ≤2 mismatches: {len(matches5)}\n")

# -----------------------------
# Question 6 - First match for AGGAGGTT with ≤2 mismatches
# -----------------------------
pattern6 = 'AGGAGGTT'
matches6 = naive_2mm(pattern6, genome)
print("Q6:")
print(f"Leftmost match with ≤2 mismatches: {matches6[0] if matches6 else 'None'}\n")

# -----------------------------
# Question 7 - Analyze FASTQ quality scores
# -----------------------------
def phred33_to_q(char):
    return ord(char) - 33

def read_fastq(filename):
    qualities = []
    with open(filename, 'r') as fh:
        while True:
            fh.readline()
            seq = fh.readline()
            fh.readline()
            qual = fh.readline().strip()
            if not qual:
                break
            qualities.append(qual)
    return qualities

def quality_by_position(qualities):
    max_len = max(len(q) for q in qualities)
    pos_scores = [Counter() for _ in range(max_len)]
    for qual in qualities:
        for i, char in enumerate(qual):
            pos_scores[i][phred33_to_q(char)] += 1
    avg_scores = []
    for counter in pos_scores:
        total = sum(score * count for score, count in counter.items())
        count = sum(counter.values())
        avg = total / count if count > 0 else 0
        avg_scores.append(avg)
    return avg_scores

fastq_path = '/Users/FxMED/PycharmProjects/Introduction to Python/Introduction/Comments/ERR037900_1.first1000.fastq'
qualities = read_fastq(fastq_path)
avg_scores = quality_by_position(qualities)

plt.figure(figsize=(12, 5))
plt.plot(range(len(avg_scores)), avg_scores)
plt.title("Average Quality Score per Sequencing Cycle")
plt.xlabel("Cycle (Base Position)")
plt.ylabel("Average Phred Score")
plt.grid(True)
plt.tight_layout()
plt.show()

worst_cycle = avg_scores.index(min(avg_scores))
print("Q7:")
print(f"Lowest average quality cycle: {worst_cycle}")