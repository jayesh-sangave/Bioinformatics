# 🧬 Computational Genomics and Advanced Bioinformatics Algorithms

This repository documents the implementation of core and advanced biological algorithms developed for the rebelScience Computational Genomics curriculum (CCA3, CCA4, and CCA5). The project demonstrates proficiency in Python programming for bioinformatics, focusing on scalable sequence manipulation, statistical analysis, pattern recognition, and protein translation.

All code adheres to strict PEP 8 standards, includes comprehensive docstrings, and is accompanied by unit tests and detailed performance analysis documents.

---

## Technical Highlights and Core Competencies

| Technical Area | Demonstrated Capability | Key Assignments |
| :--- | :--- | :--- |
| **Scalability** | Implementation of **Python Generators** and **Parallel Processing** to handle gigabyte-scale FASTA files without memory overload. | CCA5 |
| **Efficiency** | Algorithms optimized for $O(L)$ linear time complexity, avoiding quadratic complexity ($O(L^2)$) through optimized Pythonic data structures. | CCA3, CCA4 |
| **Core Algorithms** | Implementation of **Dynamic Programming** models for sequence comparison (LCS) and robust motif finding. | CCA4, CCA5 |
| **Tool Development** | Creation of a professional Command-Line Interface (CLI) using `argparse` with comprehensive error handling and logging. | CCA5 |

---

## 1. Assignment 1: DNA Fundamentals and Basic Tools (CCA3)

**Focus:** Core data structure design, sequence validation, and fundamental genetic operations.

### Theoretical Foundation

This foundational assignment establishes the computational representation of DNA sequences.

* **Sequence Validation:** Ensures data integrity by strictly accepting only standard nucleotides (A, T, G, C), preventing errors in downstream processing.
* **DNA Transcription:** The biological process of converting a DNA template strand to its corresponding functional RNA sequence by substituting **Thymine (T)** with **Uracil (U)**.
* **Reverse Complement:** A crucial step in molecular biology for studying the opposite strand, requiring efficient implementation that handles both reversal and base complementation in a single $O(L)$ pass.
* **Efficiency:** All manipulation algorithms are optimized for $O(L)$ linear time complexity (where $L$ is sequence length), leveraging `collections.Counter` and list comprehensions for high performance.

### File Structure (CCA3)

| File Name | Content Description | Assignment Questions Covered | Key Features | 
| :----- | :----- | :----- | :----- | 
| `dna_fundamentals.py` | Implementation of the `DNA` class, counting methods, and core algorithms. | Q1, Q2, Q3, Q4, Q5 | Validation, GC Content, Codon Splitting, Transcription, Reverse Complement. | 
| `testing_suite.py` | Comprehensive unit tests for all CCA3 functionality. | Q7 (Testing) | Edge cases: empty sequences, invalid input, long sequences. | 

---

## 2. Assignment 2: Advanced Genomic Analysis and Pattern Recognition (CCA4)

**Focus:** Sequence statistics, distance metrics, motif discovery, and structural analysis of protein-coding regions.

### Theoretical Foundation

This module addresses analytical and pattern-recognition challenges critical for annotation:

* **GC Skew Analysis:** Tracking the cumulative difference between G and C bases. The minimum skew point is often a strong predictor for the **Origin of Replication (oriC)** in bacterial genomes.
* **Hamming Distance & Distance Matrix:** Implemented as a measure of substitution dissimilarity, strictly defined for equal-length sequences. The distance matrix enables phylogenetic comparisons and clustering.
* **Fuzzy Matching:** Pattern recognition that allows for a defined number of mismatches (substitutions), essential for finding regulatory **motifs** that are highly similar but not perfectly conserved.
* **ORF (Open Reading Frame) Detection:** Identifying potential protein-coding regions by scanning a DNA sequence across **all six reading frames** (three forward, three reverse) to accurately locate and extract potential protein sequences between START (ATG) and STOP codons.

### File Structure (CCA4)

| File Name | Content Description | Assignment Questions Covered | Key Features | 
| :----- | :----- | :----- | :----- | 
| `cca4_algorithms.py` | GC analysis, Hamming matrix, fuzzy matching, and robust ORF translation logic. | Q1, Q2, Q3, Q4, Q5, Q6 | GC Skew, Dinucleotide Freq., Distance Matrix, Degenerate Pattern Matching. | 
| `genetic_code.py` | Dictionary constant for the Standard Genetic Code. | Q5, Q6 | Clean separation of data (codon: amino acid mapping). | 
| `cca4_testing_suite.py` | Unit tests for distance metrics and ORF extraction. | Q7 (Testing) | Tests for unequal Hamming sequences and longest ORF identification. | 

---

## 3. Assignment 3: Genomic Databases and Advanced Applications (CCA5)

**Focus:** File I/O for large data, complex algorithmic solutions, and production-ready tool packaging.

### Theoretical Foundation

This final phase focuses on producing production-ready tools capable of handling large genomic data efficiently:

* **FASTA Parsing (Generators):** Solves memory constraints by using **Python Generators** (`yield`) to read large FASTA files line-by-line, ensuring memory usage is constant regardless of file size.
* **LCS (Longest Common Subsequence):** Implemented using the classic **Dynamic Programming** approach to find the longest common string between two sequences, a core concept for sequence alignment and homology.
* **Scalability & Parallelism:** Achieved through the `concurrent.futures.ThreadPoolExecutor`, which divides time-consuming linear tasks across multiple CPU cores, dramatically reducing the wall-clock execution time for large-scale analysis.
* **Tool Packaging (CLI):** Utilizes the `argparse` module to transform the bioinformatics script into a professional, user-friendly command-line tool with integrated help and error reporting.

### File Structure (CCA5)

| File Name | Content Description | Assignment Questions Covered | Key Features | 
| :----- | :----- | :----- | :----- | 
| `cca5_file_handlers.py` | FASTA reading/writing (generator-based) and `SequenceDatabase` class. | Q1, Q2 | **Memory-efficient parsing**, data integrity checks. | 
| `cca5_advanced_algs.py` | LCS, consensus sequence generation, and conceptual suffix structure problems. | Q3, Q4 | Dynamic programming, phylogenetic readiness. | 
| `cca5_tools_cli.py` | Command-line interface logic, primary `main()` function, and parallel processing scaffolding. | Q5, Q6 | CLI setup, robust error handling, **Parallel execution**. | 
| `cca5_documentation.py` | Embedded Markdown documentation for CCA5 complexity and usage. | Q5, Q6 | Stores performance notes internally for CLI usage. | 
| `cca5_testing_suite.py` | Unit tests for file parsing, LCS, and consensus generation. | Q7 (Testing) | Focus on speed and correctness of I/O operations. | 

---

## 🚀 Usage and Testing

### Setup

```bash
# 1. Clone the repository:
git clone [Your Repository URL]
cd [Your Repository Name]

# 2. Install dependencies:
# NumPy is required for distance matrices (CCA4).
pip install numpy

# Note: Other packages (like concurrent.futures) are part of the Python standard library.
```

### Running Tests

To verify correctness and functionality for any assignment, run the associated test suite:

```bash
# Run CCA3 Tests
python dna_fundamentals.py 

# Run CCA4 Tests
python cca4_testing_suite.py

# Run CCA5 Tests (for core algorithms and file handlers)
python cca5_testing_suite.py
```

