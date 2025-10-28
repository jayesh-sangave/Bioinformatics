import collections
import random
import timeit

# --- Global Data for Q5 (Reverse Complement) ---
IUPAC_COMPLEMENT_MAP = {
    'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G',
    'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W',
    'K': 'M', 'M': 'K', 'B': 'V', 'V': 'B',
    'D': 'H', 'H': 'D', 'N': 'N', '-': '-',
}

# =================================================================
# PART A: DNA Representation and Basic Operations (Q1-Q3)
# =================================================================

class DNA:
    """
    Question 1: DNA Data Structures (15 marks)
    
    Implements a DNA class using a string internally, enforcing validation 
    and providing methods for length, counting, and basic statistics.
    """
    VALID_NUCLEOTIDES = {'A', 'T', 'G', 'C'}

    def __init__(self, sequence: str):
        """
        Initializes the DNA object and validates the sequence.
        
        :param sequence: The raw DNA sequence string.
        :raises ValueError: If the sequence contains invalid nucleotides.
        """
        # PEP 8 Note: Internal storage is prefixed with an underscore
        self._sequence = sequence.upper()
        self._validate_sequence()

    def _validate_sequence(self):
        """
        Validates that the sequence contains only A, T, G, C.
        
        Time Complexity: O(L) - where L is the sequence length.
        """
        invalid_chars = set(self._sequence) - self.VALID_NUCLEOTIDES
        if invalid_chars:
            raise ValueError(
                f"Invalid DNA sequence. Found non-ATGC nucleotides: {', '.join(sorted(invalid_chars))}"
            )

    def __str__(self):
        """Return the sequence string (Python standard method)."""
        return self._sequence

    def __len__(self):
        """Implements method for sequence length (Python standard method)."""
        return len(self._sequence)
    
    def nucleotide_count(self) -> dict:
        """
        Implements method for nucleotide counting.
        
        Time Complexity: O(L) | Space Complexity: O(1) [Fixed size dict]
        :return: A dictionary with counts for A, T, G, C.
        """
        # Utilizes collections.Counter for efficient O(L) counting
        return collections.Counter(self._sequence)

    def get_basic_statistics(self) -> dict:
        """
        Returns basic statistics: length, individual counts, and GC content.
        
        :return: Dictionary containing statistical metrics.
        """
        counts = self.nucleotide_count()
        length = len(self)
        gc_count = counts.get('G', 0) + counts.get('C', 0)
        
        stats = {
            "length": length,
            "A": counts.get('A', 0), "T": counts.get('T', 0),
            "G": counts.get('G', 0), "C": counts.get('C', 0),
            "GC_content": (gc_count / length) if length > 0 else 0.0
        }
        return stats


def count_nucleotides(sequence: str) -> dict:
    """
    Question 2: Nucleotide Counting and Analysis (12 marks)
    
    Counts individual nucleotides in a DNA sequence.
    Time Complexity: O(L)
    :param sequence: The DNA sequence string.
    :return: A dictionary with counts for A, T, G, C.
    """
    return collections.Counter(sequence.upper())

def calculate_nucleotide_frequencies(sequence: str) -> dict:
    """
    Calculates nucleotide frequencies as percentages.
    Time Complexity: O(L)
    :param sequence: The DNA sequence string.
    :return: A dictionary with percentages for A, T, G, C.
    """
    counts = count_nucleotides(sequence)
    total_length = len(sequence)
    frequencies = {}
    
    if total_length == 0: return {'A': 0.0, 'T': 0.0, 'G': 0.0, 'C': 0.0}

    for nuc in ['A', 'T', 'G', 'C']:
        frequencies[nuc] = (counts.get(nuc, 0) / total_length) * 100.0

    return frequencies

def generate_analysis_report(sequence: str) -> str:
    """
    Generates a comprehensive nucleotide analysis report.
    :param sequence: The DNA sequence string.
    :return: A formatted string report.
    """
    freqs = calculate_nucleotide_frequencies(sequence)
    gc_content = freqs.get('G', 0.0) + freqs.get('C', 0.0)
    
    report = f"--- Nucleotide Analysis Report ---\n"
    report += f"Sequence Length: {len(sequence)}\n"
    report += f"GC Content: {gc_content:.2f}%\n"
    report += f"\nBreakdown (Percentage):\n"
    for nuc, freq in freqs.items():
        report += f"  {nuc}: {freq:.2f}%\n"
    return report

def compare_nucleotide_composition(seq1: str, seq2: str) -> dict:
    """
    Compares nucleotide composition (frequencies) between two sequences (Seq1 - Seq2).
    :return: A dictionary showing the difference in frequency (Seq1 - Seq2).
    """
    freq1 = calculate_nucleotide_frequencies(seq1)
    freq2 = calculate_nucleotide_frequencies(seq2)
    
    comparison = {}
    for nuc in ['A', 'T', 'G', 'C']:
        diff = freq1.get(nuc, 0.0) - freq2.get(nuc, 0.0)
        comparison[nuc] = f"{diff:.2f}%"
    return comparison


def remove_non_nucleotide_chars(sequence: str) -> str:
    """
    Question 3: String Manipulation for Genomics (13 marks)
    
    Removes non-nucleotide characters (A, T, G, C) from a sequence.
    Time Complexity: O(L) | Uses list comprehension for efficiency.
    :param sequence: The raw sequence string.
    :return: The cleaned DNA sequence string.
    """
    valid_chars = DNA.VALID_NUCLEOTIDES
    cleaned_sequence = [char for char in sequence.upper() if char in valid_chars]
    return "".join(cleaned_sequence)

def split_into_codons(sequence: str) -> list[str]:
    """
    Splits a long, cleaned sequence into codons (groups of 3).
    Time Complexity: O(L) | Space Complexity: O(L)
    :return: A list of 3-base codon strings.
    """
    cleaned_seq = remove_non_nucleotide_chars(sequence)
    codons = [cleaned_seq[i:i + 3] for i in range(0, len(cleaned_seq), 3)]
    return codons

def merge_dna_fragments(fragments: list[str]) -> str:
    """
    Merges multiple DNA fragments into a single sequence.
    Time Complexity: O(L_total) | Uses efficient str.join().
    :param fragments: A list of DNA sequence strings.
    :return: The single merged sequence string.
    """
    return "".join(fragments)


# =================================================================
# PART B: Essential DNA Algorithms (Q4-Q5)
# =================================================================

def transcribe_dna_to_rna(dna_sequence: str, strand_type: str = 'coding') -> str:
    """
    Question 4: DNA Transcription (15 marks)
    
    Converts a DNA sequence to an RNA sequence (T -> U substitution).
    Handles 'coding' (same sequence, T->U) and 'template' (complementary, T->U) strands.
    Time Complexity: O(L)
    """
    # Error checking for invalid sequence (removes non-ATGC)
    dna_sequence = remove_non_nucleotide_chars(dna_sequence)
    if not dna_sequence: return ""

    if strand_type == 'coding':
        # T -> U
        rna_sequence = dna_sequence.replace('T', 'U')
    elif strand_type == 'template':
        # First, complement the template strand (using standard A/T/G/C mapping)
        complement_map = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
        template_strand_comp = "".join([complement_map.get(base, '') for base in dna_sequence])
        # Then, T -> U
        rna_sequence = template_strand_comp.replace('T', 'U')
    else:
        # Error checking for unknown strand type
        raise ValueError("Invalid strand_type. Must be 'coding' or 'template'.")

    return rna_sequence

def batch_transcription(dna_sequences: list[str], strand_type: str = 'coding') -> list[str]:
    """
    Performs batch processing for multiple sequences.
    Time Complexity: O(N * L_avg)
    """
    rna_sequences = []
    for seq in dna_sequences:
        try:
            rna_sequences.append(transcribe_dna_to_rna(seq, strand_type))
        except ValueError:
            rna_sequences.append(None) # Mark as failed
    return rna_sequences


def generate_reverse_complement(
    dna_sequence: str, 
    orientation: str = "5'-3'",
    optimized: bool = True
) -> str:
    """
    Question 5: Reverse Complement Generation (20 marks)
    
    Generates the reverse complement, handling standard and degenerate IUPAC codes.
    The 'optimized' flag demonstrates the performance analysis requirement (Q6).
    
    Time Complexity: O(L) (one pass for complement, one pass for reversal/join)
    Space Complexity: O(L) (to store intermediate character list)
    
    :param dna_sequence: The input DNA sequence string.
    :param orientation: "5'-3'" (default) or "3'-5'".
    :param optimized: Use efficient list comprehension and join/reversed.
    :return: The reverse complement sequence string.
    """
    if not dna_sequence: return ""
    
    # 1. Complement Generation (Robust and handles degenerate bases)
    # Optimized: Use list comprehension for speed
    if optimized:
        complement_chars = [IUPAC_COMPLEMENT_MAP.get(base, base) for base in dna_sequence.upper()]
    else:
        # Non-optimized version (for Q6 profiling)
        complement_chars = []
        for base in dna_sequence.upper():
             complement_chars.append(IUPAC_COMPLEMENT_MAP.get(base, base))

    # 2. Handle Orientation
    if orientation == "5'-3'":
        # Reverse the complement
        reverse_complement = "".join(reversed(complement_chars))
        return reverse_complement
    elif orientation == "3'-5'":
        # This is the complement sequence without reversing it
        return "".join(complement_chars)
    else:
        raise ValueError("Invalid orientation. Must be '5'-3'' or '3'-5''.")

# --- Example Usage ---
if __name__ == '__main__':
    try:
        dna_obj = DNA("ATGCc_ATXGC")
        print(f"Q1 DNA Length: {len(dna_obj)}")
        print(f"Q2 Frequencies: {calculate_nucleotide_frequencies(str(dna_obj))}")
        print(f"Q3 Codons: {split_into_codons(str(dna_obj))}")
        print(f"Q4 RNA: {transcribe_dna_to_rna(str(dna_obj))}")
        print(f"Q5 Reverse Complement: {generate_reverse_complement(str(dna_obj))}")
    except ValueError as e:
        print(f"Caught expected validation error: {e}")