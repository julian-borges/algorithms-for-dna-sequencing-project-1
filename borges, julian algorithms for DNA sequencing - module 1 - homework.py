'''
# Bioinformatics - Genomic Data Science
# Algorithms for DNA Sequencing
# The Johns Hopkins University
# Author: Julian Borges
# Programming Homework 1
'''
# -----------------------------
# Question 1
# -----------------------------
'''How many times does AGGT or its reverse complement (ACCT)
# occur in the lambda virus genome?'''
# -----------------------------

# -----------------------------
# Step 1: Read genome from FASTA
# -----------------------------
def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

# -----------------------------
# Local path to your genome
# -----------------------------
genome_path = '/Users/FxMED/PycharmProjects/Introduction to Python/Introduction/Comments/lambda_virus.fa'
genome = read_genome(genome_path)

# -----------------------------
# Step 2: Reverse complement
# -----------------------------
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# -----------------------------
# Step 3: Count matches
# -----------------------------
def count_occurrences(pattern, sequence):
    count = 0
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i + len(pattern)] == pattern:
            count += 1
    return count

# -----------------------------
# Step 4: Run count for 'AGGT' and 'ACCT'
# -----------------------------
pattern = 'AGGT'
revcomp = reverse_complement(pattern)

count_forward = count_occurrences(pattern, genome)
count_reverse = count_occurrences(revcomp, genome)

# -----------------------------
# Step 5: Results
# -----------------------------
print(f"Forward pattern '{pattern}': {count_forward} occurrence(s)")
print(f"Reverse complement '{revcomp}': {count_reverse} occurrence(s)")
print(f"Total matches (forward or reverse): {count_forward + count_reverse}")

# -----------------------------
# Question 2
# -----------------------------
'''How many times does TTAA,or its reverse complement occur in the lambda virus genome?'''
# -----------------------------

# -----------------------------
# Step 1: Read genome from FASTA
# -----------------------------
def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

# Your genome file
genome_path = '/Users/FxMED/PycharmProjects/Introduction to Python/Introduction/Comments/lambda_virus.fa'
genome = read_genome(genome_path)

# -----------------------------
# Step 2: Reverse complement function
# -----------------------------
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# -----------------------------
# Step 3: Count pattern occurrences
# -----------------------------
def count_occurrences(pattern, sequence):
    count = 0
    for i in range(len(sequence) - len(pattern) + 1):
        if sequence[i:i + len(pattern)] == pattern:
            count += 1
    return count

# -----------------------------
# Step 4: Count TTAA and its reverse
# -----------------------------
motif = 'TTAA'
revcomp = reverse_complement(motif)

count_forward = count_occurrences(motif, genome)
count_reverse = count_occurrences(revcomp, genome)

# -----------------------------
# Step 5: Output results
# -----------------------------
print(f"Motif '{motif}' occurs {count_forward} time(s)")
print(f"Reverse complement '{revcomp}' occurs {count_reverse} time(s)")
print(f"Total (forward or reverse): {count_forward + count_reverse}")

# -----------------------------
# Question 3
# -----------------------------
'''What is the offset of the leftmost occurrence of ACTAAGT or its reverse
# complement in the Lambda virus genome? E.g. if the leftmost occurrence of
# ACTAAGT is at offset 40 (0-based) and the leftmost occurrence of its reverse
# complement ACTTAGT is at offset 29, then report 29.
'''
# -----------------------------

# -----------------------------
# Step 1: Read genome from FASTA
# -----------------------------

def read_genome(filename):
    genome = ''
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

# -----------------------------
# Step 2: Reverse complement function
# -----------------------------
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# -----------------------------
# Step 3: Naive match finder - returns list of all match positions
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
# Loading genome from local file
# -----------------------------
genome_path = '/Users/FxMED/PycharmProjects/Introduction to Python/Introduction/Comments/lambda_virus.fa'
genome = read_genome(genome_path)

# -----------------------------
# Define the pattern and compute reverse complement
# -----------------------------
pattern = 'ACTAAGT'
revcomp = reverse_complement(pattern)

# -----------------------------
# Find all positions
# -----------------------------
positions_pattern = naive_match(pattern, genome)
positions_revcomp = naive_match(revcomp, genome)

# -----------------------------
# Find leftmost match
# -----------------------------
min_pos = min(
    positions_pattern[0] if positions_pattern else float('inf'),
    positions_revcomp[0] if positions_revcomp else float('inf')
)

# -----------------------------
# Print result
# -----------------------------
if min_pos != float('inf'):
    print(f"Leftmost occurrence of '{pattern}' or its reverse complement '{revcomp}' is at offset {min_pos}")
else:
    print("Pattern and its reverse complement were not found in the genome.")

# -----------------------------
# Question 4
# -----------------------------
'''What is the offset of the leftmost occurrence of AGTCGA or its reverse
complement in the Lambda virus genome?'''
# -----------------------------

# -----------------------------
# Load the genome from FASTA
# -----------------------------
def read_genome(filepath):
    genome = ''
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

# -----------------------------
# Reverse complement function
# -----------------------------
def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

# -----------------------------
# Naive pattern matcher
# -----------------------------
def naive_match(pattern, text):
    for i in range(len(text) - len(pattern) + 1):
        match = True
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                match = False
                break
        if match:
            return i  # Return only first match
    return None

# -----------------------------
# Apply to genome + pattern
# -----------------------------
genome_path = '/Users/FxMED/PycharmProjects/Introduction to Python/Introduction/Comments/lambda_virus.fa'
genome = read_genome(genome_path)

pattern = 'AGTCGA'
rev_pattern = reverse_complement(pattern)

pos_pattern = naive_match(pattern, genome)
pos_revcomp = naive_match(rev_pattern, genome)

# Compare positions
if pos_pattern is not None and pos_revcomp is not None:
    result = min(pos_pattern, pos_revcomp)
elif pos_pattern is not None:
    result = pos_pattern
elif pos_revcomp is not None:
    result = pos_revcomp
else:
    result = "Pattern not found in genome."

# Print result
print(f"Leftmost occurrence of '{pattern}' or its reverse complement '{rev_pattern}' is at offset: {result}")

# -----------------------------
# Question 4
# -----------------------------
'''As we will discuss, sometimes we would like to find approximate matches for Pin T. 
That is, we want to find occurrences with one or more differences.
For Questions 5 and 6, make a new version of the naive function called naive _2mm that 
allows up to 2 mismatches per occurrence. Unlike for the previous questions, 
do not consider the reverse complement here. We're looking for approximate matches for 
Pitself, not its reverse complement.

For example, ACTTTA occurs twice in ACTTACTTGATAAAGT, once at offset 0 with 2 mismatches, 
and once at offset 4 with 1 mismatch. So naive_2mm ('ACITTA', 'ACITACITGATAAAGI') should 
return the list [O, 4].
Hint: See the attached notebook for a few examples i used to test the naive_2mm function.
How many times does TTCAAGCC occur in the Lambda virus genome when allowing up to 2 
mismatches?'''
# -----------------------------

# -----------------------------
# Step 1: Load the genome
def read_genome(filepath):
    genome = ''
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

# Step 2: Define naive_2mm function
def naive_2mm(p, t):
    occurrences = []
    for i in range(len(t) - len(p) + 1):
        mismatches = 0
        for j in range(len(p)):
            if t[i + j] != p[j]:
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)
    return occurrences

# Step 3: Load genome and search
genome_path = '/Users/FxMED/PycharmProjects/Introduction to Python/Introduction/Comments/lambda_virus.fa'
genome = read_genome(genome_path)

pattern = 'TTCAAGCC'
approx_matches = naive_2mm(pattern, genome)

# Step 4: Print result
print(f"Occurrences of '{pattern}' with ≤2 mismatches: {len(approx_matches)}")

# -----------------------------
# Question 5
# -----------------------------
'''What is the offset of the leftmost occurrence of AGGAGGIT in the Lambda virus 
genome when allowing up to 2 mismatches?'''
# -----------------------------
# ------------------------------------------------------------
# Step 1: Read genome from local FASTA
# ------------------------------------------------------------
def read_genome(filepath):
    genome = ''
    with open(filepath, 'r') as f:
        for line in f:
            if not line.startswith('>'):
                genome += line.strip()
    return genome

# ------------------------------------------------------------
# Step 2: Define naive_2mm matcher
# ------------------------------------------------------------
def naive_2mm(pattern, text):
    occurrences = []
    for i in range(len(text) - len(pattern) + 1):
        mismatches = 0
        for j in range(len(pattern)):
            if text[i + j] != pattern[j]:
                mismatches += 1
                if mismatches > 2:
                    break
        if mismatches <= 2:
            occurrences.append(i)
    return occurrences

# ------------------------------------------------------------
# Step 3: Run search
# ------------------------------------------------------------
genome_path = '/Users/FxMED/PycharmProjects/Introduction to Python/Introduction/Comments/lambda_virus.fa'
genome = read_genome(genome_path)

pattern = 'AGGAGGTT'  # The pattern to search
matches = naive_2mm(pattern, genome)

# ------------------------------------------------------------
# Step 4: Report result
# ------------------------------------------------------------
if matches:
    print(f"Leftmost occurrence of '{pattern}' with ≤2 mismatches is at offset: {matches[0]}")
else:
    print(f"No match found for '{pattern}' with ≤2 mismatches.")


# -----------------------------
# Question 7
# -----------------------------
'''Finally, download and parse the provided FAST file containing real DNA sequencing reads 
derived from a human:
https://d28rh4a8wqOiu5.cloudfront.net/ads1/data/ERR0379001.first1000.fasto E3
Note that the file has many reads in it and you should examine all of them together when 
answering this question. The reads are taken from this study:

Ajay, S. S., Parker, S. C., Abaan, H. O., Fajardo, K. V. F., & Margulies, E. H. (2011).
 
Accurate and comprehensive sequencing of personal genomes. Genome research, 21(9), 1498-1505.

This dataset has something wrong with it; one of the sequencing cycles is poor quality.
Report which sequencing cycle has the problem. Remember that a sequencing cycle corresponds 
to a particular offset in all the reads. For example, if the leftmost read position seems 
to have a problem consistently across reads, report o. If the fourth position from the 
left has the problem, report 3. Do whatever analysis you think is needed to identify the 
bad cycle.'''
# -----------------------------
import matplotlib.pyplot as plt
from collections import Counter

# -----------------------------
# Step 1: Phred33 Quality Decoder
# -----------------------------
def phred33_to_q(char):
    return ord(char) - 33

# -----------------------------
# Step 2: Read .fastq File
# -----------------------------
def read_fastq(filename):
    qualities = []
    with open(filename, 'r') as fh:
        while True:
            fh.readline()           # name
            seq = fh.readline()    # sequence
            fh.readline()          # plus line
            qual = fh.readline().strip()  # quality line
            if not qual:
                break
            qualities.append(qual)
    return qualities

# -----------------------------
# Step 3: Analyze Per-Position Quality
# -----------------------------
def quality_by_position(qualities):
    max_len = max(len(q) for q in qualities)
    pos_scores = [Counter() for _ in range(max_len)]

    for qual in qualities:
        for i, char in enumerate(qual):
            pos_scores[i][phred33_to_q(char)] += 1

    avg_scores = []
    for counter in pos_scores:
        total_score = sum(score * count for score, count in counter.items())
        total_count = sum(counter.values())
        avg = total_score / total_count if total_count > 0 else 0
        avg_scores.append(avg)

    return avg_scores

# -----------------------------
# Step 4: Run Analysis
# -----------------------------
fastq_path = '/Users/FxMED/PycharmProjects/Introduction to Python/Introduction/Comments/ERR037900_1.first1000.fastq'
qualities = read_fastq(fastq_path)
avg_qual = quality_by_position(qualities)

# -----------------------------
# Step 5: Plot Results
# -----------------------------
plt.figure(figsize=(12, 5))
plt.plot(range(len(avg_qual)), avg_qual)
plt.title("Average Quality Score per Sequencing Cycle")
plt.xlabel("Cycle (Base Position)")
plt.ylabel("Average Phred Score")
plt.grid(True)
plt.tight_layout()
plt.show()

# -----------------------------
# Step 6: Report Worst Cycle
# -----------------------------
worst = avg_qual.index(min(avg_qual))
print(f"❗ Lowest average quality at cycle: {worst} (0-based index)")
