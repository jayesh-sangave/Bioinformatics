import unittest
import random
# Import all functions and classes from the main file
from dna_fundamentals import (
    DNA, generate_reverse_complement, transcribe_dna_to_rna, 
    split_into_codons, calculate_nucleotide_frequencies
)

# --- Test Data Generator (Q7) ---
def generate_random_dna(length: int, include_invalid: bool = False) -> str:
    """
    Generates a random DNA sequence of a given length.
    """
    chars = ['A', 'T', 'G', 'C']
    if include_invalid:
        chars.extend(['X', 'Z', '-'])
    return "".join(random.choices(chars, k=length))

class TestDNAFunctions(unittest.TestCase):
    """
    Question 7: Comprehensive Testing Suite (10 marks)
    
    Provides unit tests for all major DNA manipulation functions and edge cases.
    """
    
    # === Q1/Q3: DNA Class and Basic Ops Tests ===
    def test_dna_validation(self):
        """Test sequence validation and error handling."""
        # Valid sequence test
        dna_valid = DNA("ATGCAGTA")
        self.assertEqual(len(dna_valid), 8)
        
        # Invalid sequence (Error Handling - Edge Case)
        with self.assertRaisesRegex(ValueError, "Found non-ATGC"):
            DNA("ATXGC-")

    def test_basic_statistics(self):
        """Test nucleotide counting and statistics."""
        dna_test = DNA("AATTGGCC")
        stats = dna_test.get_basic_statistics()
        self.assertEqual(stats['A'], 2)
        self.assertEqual(stats['GC_content'], 0.5) # 4/8 = 0.5

    # === Q4: Transcription Tests ===
    def test_dna_transcription(self):
        """Test DNA to RNA conversion."""
        # Coding strand (T -> U)
        self.assertEqual(transcribe_dna_to_rna("ATGCCGT"), "AUGCCGU")
        
        # Template strand (Complement, then T -> U)
        # Template: ACGT -> Complement: TGCA -> RNA: UGCA
        self.assertEqual(transcribe_dna_to_rna("ACGT", strand_type='template'), "UGCA")
        
    def test_transcription_empty_sequence(self):
        """Edge case: Empty sequence."""
        self.assertEqual(transcribe_dna_to_rna(""), "")

    # === Q5: Reverse Complement Tests ===
    def test_reverse_complement_standard(self):
        """Test standard reverse complement functionality (5'-3')."""
        self.assertEqual(generate_reverse_complement("ATGC"), "GCAT")
        self.assertEqual(generate_reverse_complement("T"), "A") # Single nucleotide

    def test_reverse_complement_degenerate(self):
        """Test handling of degenerate IUPAC codes."""
        # R->Y, Y->R, K->M, W->W, N->N. Sequence: RYKWN -> Complement: YRMWM -> Reverse: MWMYR
        self.assertEqual(generate_reverse_complement("RYKWN"), "NMWYR")

    def test_reverse_complement_orientation(self):
        """Test 3'-5' orientation (complement only)."""
        # Sequence: ATGC -> Complement: TACG -> Reverse Complement: GCAT
        self.assertEqual(generate_reverse_complement("ATGC", orientation="3'-5'"), "TACG")

    # === Q3/Q7: Codon Splitting and Long Sequence Test ===
    def test_codon_splitting(self):
        """Test splitting sequences into groups of 3."""
        # Full codons and partial end
        self.assertEqual(split_into_codons("ATGCATGC"), ['ATG', 'CAT', 'GC'])
        
        # Robustness test (invalid chars should be removed before splitting)
        self.assertEqual(split_into_codons("ATG-C-T"), ['ATG', 'CT'])

    def test_long_sequence_performance(self):
        """Test functionality and performance with a very long random sequence."""
        long_seq = generate_random_dna(10000)
        # Check that the reverse complement is of the correct length
        self.assertEqual(len(generate_reverse_complement(long_seq)), 10000)

if __name__ == '__main__':
    unittest.main()