# cca4_testing_suite.py

import unittest
import numpy as np
import random
from cca4_algorithms import (
    calculate_gc_content, analyze_gc_content_sliding_window, 
    calculate_hamming_distance, calculate_distance_matrix, 
    find_fuzzy_pattern_matches, find_longest_orf, translate_sequence
)

# Test Data Helper
def generate_random_dna(length: int) -> str:
    return "".join(random.choices(['A', 'T', 'G', 'C'], k=length))

class TestCCA4Algorithms(unittest.TestCase):
    
    # === Q1/Q2 GC and Composition Tests ===
    def test_gc_content(self):
        self.assertAlmostEqual(calculate_gc_content("ATGC"), 50.0)
        self.assertAlmostEqual(calculate_gc_content("GGCGCGC"), 85.71, places=2)
        self.assertEqual(calculate_gc_content(""), 0.0) # Edge case: empty sequence

    def test_sliding_window(self):
        seq = "GGCAATGCC"
        results = analyze_gc_content_sliding_window(seq, 3, step=1)
        # (0, 66.67), (1, 33.33), (2, 33.33), (3, 33.33), (4, 66.67), (5, 100.0), (6, 66.67)
        self.assertEqual(len(results), 7)
        self.assertAlmostEqual(results[5][1], 100.0)

    # === Q3 Hamming Distance Tests ===
    def test_hamming_distance_equal_length(self):
        self.assertEqual(calculate_hamming_distance("GAGCCTTA", "GATCTTAC"), 3)
        self.assertEqual(calculate_hamming_distance("AAAA", "AAAA"), 0)

    def test_hamming_distance_unequal_length(self):
        # Edge case: Different lengths
        with self.assertRaises(ValueError):
            calculate_hamming_distance("ATGC", "ATG")
            
    def test_distance_matrix(self):
        seqs = ["ATGC", "ATTC", "GATA"]
        matrix = calculate_distance_matrix(seqs)
        # D(ATGC, ATTC) = 1; D(ATGC, GATA) = 3; D(ATTC, GATA) = 4
        expected = np.array([[0, 1, 3], [1, 0, 4], [3, 4, 0]])
        self.assertTrue(np.array_equal(matrix, expected))

    # === Q4 Motif Finding Tests ===
    def test_fuzzy_matching(self):
        text = "GATTACA"
        pattern = "GATA"
        # GATA (0 mismatches), TACA (1 mismatch)
        matches = find_fuzzy_pattern_matches(text, pattern, max_mismatches=1)
        self.assertEqual(len(matches), 2)
        self.assertIn((0, 0), matches)
        self.assertIn((3, 1), matches)

    # === Q5/Q6 ORF and Translation Tests ===
    def test_translate_sequence(self):
        # Start (ATG), Amino Acid (TTC=F), Stop (TGA=*)
        self.assertEqual(translate_sequence("ATGTTCGATGA"), "M F D *")
        
    def test_find_longest_orf(self):
        # Sequence: ATG TTT TGA GGG ATG AAA AAA TAG
        # Frame +1 ORFs: M, MKK
        seq = "ATgTTTtGagGGATgAAAaAaTAG"
        longest_orf, frame, length = find_longest_orf(seq)
        
        # Expected longest ORF is 'MKK' (length 3)
        self.assertEqual(longest_orf, "MKK") 
        self.assertEqual(length, 3)
        self.assertEqual(frame, '+1')
        
if __name__ == '__main__':
    unittest.main()