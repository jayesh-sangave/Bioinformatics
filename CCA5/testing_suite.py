import unittest
import os
import io
import tempfile
from cca5_file_handlers import parse_fasta_generator, write_fasta, SequenceDatabase
from cca5_advanced_algs import (
    find_longest_common_subsequence, 
    generate_consensus_sequence,
    find_repeats_by_suffix_array_concept
)

# --- Test Data Setup ---
TEST_FASTA_CONTENT = """
>Seq1 Example header
ATGCAGTA
ATGCAGTA
>Seq2 Different
TACCATGC
GTAG
>Empty Seq
>Invalid Seq
ATXC
"""

class TestCCA5FileHandlers(unittest.TestCase):
    """Q1, Q2: Tests for file parsing, writing, and database functions."""

    def setUp(self):
        # Create a temporary FASTA file for testing
        self.temp_file = tempfile.NamedTemporaryFile(mode='w+', delete=False, suffix='.fasta')
        self.temp_file.write(TEST_FASTA_CONTENT)
        self.temp_file.close()
        self.file_path = self.temp_file.name

    def tearDown(self):
        # Clean up the temporary file
        os.unlink(self.file_path)

    def test_parse_fasta_generator_correctness(self):
        """Test Q1: Correct header/sequence extraction and generator function."""
        sequences = list(parse_fasta_generator(self.file_path))
        
        # Check total sequences parsed
        self.assertEqual(len(sequences), 4) # Seq1, Seq2, Empty, Invalid
        
        # Check sequence length and content (Seq1 is wrapped)
        self.assertEqual(sequences[0][0], "Seq1 Example header")
        self.assertEqual(sequences[0][1], "ATGCAGTAATGCAGTA")
        
        # Check sequence 3 (Empty Seq) - should be empty sequence string
        self.assertEqual(sequences[2][0], "Empty Seq")
        self.assertEqual(sequences[2][1], "")

    def test_sequence_database_load_and_validate(self):
        """Test Q2: Database loading and basic validation."""
        db = SequenceDatabase()
        count = db.load_from_fasta(self.file_path)
        
        self.assertEqual(count, 4)
        
        # Validation should find 1 invalid sequence (ATXC)
        invalid_count = db.validate_data()
        self.assertEqual(invalid_count, 1)

class TestCCA5AdvancedAlgs(unittest.TestCase):
    """Q3, Q4: Tests for advanced sequence problems."""

    def test_longest_common_subsequence(self):
        """Test Q3: LCS using Dynamic Programming."""
        self.assertEqual(find_longest_common_subsequence("GATTACA", "GATCTAC"), "GATTAC")
        self.assertEqual(find_longest_common_subsequence("AABBC", "BBCCA"), "BBCC")
        self.assertEqual(find_longest_common_subsequence("AGGTAB", "GXTXAYB"), "GATAB")

    def test_generate_consensus_sequence(self):
        """Test Q3: Consensus sequence and profile matrix generation."""
        sequences = ["AATTC", "AGATC", "TATTC", "AAGTC"]
        consensus, profile = generate_consensus_sequence(sequences)
        
        # Pos 0: A (3), T (1) -> A
        # Pos 1: A (2), G (2) -> A or G (depends on tie-breaking, max() is safe) -> A
        # Pos 2: T (3), G (1) -> T
        # Pos 3: T (2), C (2) -> T
        # Pos 4: C (3), T (1) -> C
        self.assertEqual(consensus, "AATTC") 
        self.assertEqual(len(profile), 5)
        self.assertEqual(profile[0]['A'], 0.75) # 3/4 A at pos 0

    def test_find_repeats_by_suffix_array_concept(self):
        """Test Q4: Conceptual suffix array for repeat finding."""
        # Simple repeat: ATGCATGC
        repeats = find_repeats_by_suffix_array_concept("ATGCATGC")
        
        # Expected to find ATGC repeat (length 4) or longer common prefixes
        self.assertTrue(any(r == 'ATGC' for r, l in repeats))
        
        # Repeat longer than 5: BANANABANANA
        repeats_banana = find_repeats_by_suffix_array_concept("BANANABANANA")
        self.assertTrue(any(r == 'BANANA' for r, l in repeats_banana))


if __name__ == '__main__':
    unittest.main()
