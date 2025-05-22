# algorithms-for-dna-sequencing-project-1
Bioinformatics - Genomic Data Science Specialization - Module 1 
Project for the Algorithms for DNA Sequencing course by Johns Hopkins University. 
Includes FASTA/FASTQ parsing, naive exact and approximate matching, and read quality analysis.

# DNA Sequencing Algorithms – JHU Genomic Data Science Specialization

**Author:** Julian Borges MD
**Course:** Algorithms for DNA Sequencing  
**Institution:** Johns Hopkins University  

## Project Description

This project implements core algorithms and analyses from **Module 1** of the JHU Genomic Data Science specialization. It covers:

- Reading and parsing **FASTA** and **FASTQ** files
- Performing **naive exact string matching**
- Implementing a version of the matcher allowing **≤2 mismatches**
- Finding the **reverse complements** of DNA sequences
- Analyzing **sequencing quality scores** per read cycle
- Answering biological questions about specific patterns in the **lambda virus genome**

---

## Project Structure

```bash
algorithms-dna-sequencing-project-1/
├── dna_sequencing_analysis.py        # Main script with all answers
├── dna_sequencing_analysis.ipynb     # Optional: Notebook version
├── lambda_virus.fa                   # Reference genome (FASTA)
├── ERR037900_1.first1000.fastq       # Human reads dataset (FASTQ)
├── README.md                         # This file
└── requirements.txt                  # Python dependencies
