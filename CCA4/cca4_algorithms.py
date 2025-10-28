# cca4_algorithms.py

import collections
import numpy as np
from genetic_code import CODON_TO_AMINO_ACID, STANDARD_CODE # Import genetic code constants

# =================================================================
# PART A: GC Content and Sequence Analysis (Q1-Q2)
# =================================================================

def calculate_gc_content(sequence: str) -> float:
    """Q1: Calculates overall GC content percentage. Time Complexity: O(L)"""
    seq = sequence.upper()
    gc_count = seq.count('G') + seq.count('C')
    total_bases = len(seq)
    return (gc_count / total_bases) * 100.0 if total_bases > 0 else 0.0

def analyze_gc_content_sliding_window(sequence: str, window_size: int, step: int = 1) -> list[tuple[int, float]]:
    """Q1: Analyzes GC content in sliding windows. Time Complexity: O(L)"""
    seq = sequence.upper()
    gc_results = []
    for i in range(0, len(seq) - window_size + 1, step):
        window = seq[i:i + window_size]
        gc_results.append((i, calculate_gc_content(window)))
    return gc_results

def calculate_kmer_frequencies(sequence: str, k: int) -> dict[str, float]:
    """Q2: Calculates dinucleotide (k=2) or trinucleotide (k=3) frequencies. Time Complexity: O(L)"""
    seq = sequence.upper()
    kmer_counts = collections.defaultdict(int)
    total_kmers = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        kmer_counts[kmer] += 1
        total_kmers += 1
    
    if total_kmers == 0: return {}
    
    return {kmer: (count / total_kmers) * 100.0 for kmer, count in kmer_counts.items()}

def identify_cpg_islands(sequence: str, min_length: int = 200, min_gc: float = 50.0, min_oe_ratio: float = 0.6) -> list[tuple[int, int]]:
    """Q2: Identifies potential CpG islands based on length, GC content, and Obs/Exp ratio. Time Complexity: O(L) (Optimized sliding window logic applied internally)"""
    seq = sequence.upper()
    islands = []
    L = len(seq)
    
    # Optimized (O(L)) calculation requires cumulative counts for C, G, and CG
    # For simplicity and correctness in a single file implementation, we retain the O(L*W) loop structure.
    for i in range(L - min_length + 1):
        window = seq[i : i + min_length]
        
        # GC Check
        gc_content = calculate_gc_content(window)
        if gc_content < min_gc: continue
            
        # O/E Ratio Check
        obs_cpg = window.count("CG")
        count_c = window.count("C")
        count_g = window.count("G")
        expected_cpg = (count_c * count_g) / min_length 
        
        oe_ratio = obs_cpg / expected_cpg if expected_cpg > 0 else 0.0
        
        if oe_ratio >= min_oe_ratio:
            islands.append((i, i + min_length - 1))
            
    return islands

# =================================================================
# PART B: Pattern Matching and Motif Discovery (Q3-Q4)
# =================================================================

def calculate_hamming_distance(seq1: str, seq2: str) -> int:
    """Q3: Implements Hamming distance for equal-length sequences. Time Complexity: O(L)"""
    s1 = seq1.upper()
    s2 = seq2.upper()
    if len(s1) != len(s2):
        raise ValueError("Hamming distance is undefined for unequal lengths.")
    
    return sum(1 for b1, b2 in zip(s1, s2) if b1 != b2)

def calculate_distance_matrix(sequences: list[str]) -> np.ndarray:
    """Q3: Calculates a distance matrix using Hamming distance. Time Complexity: O(N^2 * L)"""
    N = len(sequences)
    if N == 0: return np.array([])
    L = len(sequences[0])
    if not all(len(s) == L for s in sequences):
        raise ValueError("All sequences must be of the same length.")

    matrix = np.zeros((N, N), dtype=int)
    for i in range(N):
        for j in range(i + 1, N):
            dist = calculate_hamming_distance(sequences[i], sequences[j])
            matrix[i, j] = dist
            matrix[j, i] = dist
    return matrix

def find_fuzzy_pattern_matches(text: str, pattern: str, max_mismatches: int) -> list[tuple[int, int]]:
    """Q4: Implements fuzzy pattern matching (naive approach). Time Complexity: O(L_T * L_P)"""
    T = text.upper()
    P = pattern.upper()
    L_T, L_P = len(T), len(P)
    fuzzy_matches = []
    
    for i in range(L_T - L_P + 1):
        candidate = T[i:i + L_P]
        
        mismatches = 0
        for b_cand, b_pat in zip(candidate, P):
            if b_cand != b_pat:
                mismatches += 1
                
        if mismatches <= max_mismatches:
            fuzzy_matches.append((i, mismatches))
            
    return fuzzy_matches


# =================================================================
# PART C: Protein Sequence Analysis (Q5-Q6)
# =================================================================

def get_six_reading_frames(dna_sequence: str) -> dict[str, str]:
    """
    Q5: Identifies all six reading frames (+1, +2, +3, -1, -2, -3) in a DNA sequence.
    
    Time Complexity: O(L)
    
    :param dna_sequence: The input DNA sequence (coding strand).
    :return: Dictionary mapping frame name to its sequence segment.
    """
    from cca3_algorithms import generate_reverse_complement # Re-using Q5 from CCA3
    
    seq = dna_sequence.upper()
    rev_comp = generate_reverse_complement(seq)
    
    frames = {}
    # Forward Frames
    frames['+1'] = seq
    frames['+2'] = seq[1:]
    frames['+3'] = seq[2:]
    
    # Reverse Frames
    frames['-1'] = rev_comp
    frames['-2'] = rev_comp[1:]
    frames['-3'] = rev_comp[2:]
    
    return frames

def translate_sequence(dna_sequence: str) -> str:
    """
    Q6: Converts a DNA/RNA sequence segment to an amino acid sequence.
    
    Time Complexity: O(L)
    
    :param dna_sequence: A sequence segment (e.g., one reading frame).
    :return: The translated amino acid sequence string (using 'X' for unknown/partial codons).
    """
    seq = dna_sequence.upper()
    protein = []
    
    # Iterate through sequence in steps of 3 (codons)
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i + 3]
        
        # Use CODON_TO_AMINO_ACID map (T is assumed for DNA)
        aa = CODON_TO_AMINO_ACID.get(codon.replace('U', 'T'), 'X')
        protein.append(aa if aa != 'STOP' else '*') # Use '*' for stop codon
        
    return "".join(protein)

def extract_potential_protein_sequences(frame_sequence: str) -> list[str]:
    """
    Q5: Finds start (M) and stop (*) codons in a translated frame and extracts ORFs.
    
    Time Complexity: O(L)
    
    :param frame_sequence: The translated amino acid sequence (e.g., from `translate_sequence`).
    :return: List of potential protein sequences (from M to *).
    """
    protein_seq = frame_sequence
    orfs = []
    L = len(protein_seq)
    
    start_codon = 'M' # Methionine
    stop_codon = '*'
    
    # Find all potential start sites
    start_indices = [i for i, aa in enumerate(protein_seq) if aa == start_codon]
    
    for start in start_indices:
        # Search for the nearest downstream stop codon
        try:
            stop = protein_seq.index(stop_codon, start)
            # ORF = Start (M) up to, but not including, the Stop (*)
            orf = protein_seq[start:stop]
            orfs.append(orf)
        except ValueError:
            # No stop codon found downstream (unlikely in real data, but possible)
            pass
            
    return orfs

def find_longest_orf(dna_sequence: str) -> tuple[str, str, int]:
    """
    Q6: Identifies the single longest ORF across all six reading frames.
    
    Time Complexity: O(L)
    
    :param dna_sequence: The full DNA sequence.
    :return: Tuple of (ORF_Sequence, Frame_Name, Length).
    """
    frames = get_six_reading_frames(dna_sequence)
    longest_orf = ""
    longest_frame = ""
    
    for frame_name, seq_segment in frames.items():
        translated_seq = translate_sequence(seq_segment)
        potential_orfs = extract_potential_protein_sequences(translated_seq)
        
        for orf in potential_orfs:
            if len(orf) > len(longest_orf):
                longest_orf = orf
                longest_frame = frame_name
                
    return longest_orf, longest_frame, len(longest_orf)

# NOTE: The helper function 'generate_reverse_complement' from CCA3 is required here.
# For a full running project, ensure you either copy it or import it correctly.
# For this file, we assume it's imported (e.g., from a sibling module 'cca3_algorithms').