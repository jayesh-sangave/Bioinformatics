# CCA4 Assignment 2: Advanced Genomic Analysis and Pattern Recognition

## Algorithm Efficiency and Performance Analysis (Q6)

The primary goal for CCA4's algorithms was to prioritize time complexity, especially for analysis functions operating on large genomic sequences.

### Time and Space Complexity Summary

| Function | Time Complexity | Space Complexity | Notes on Efficiency |
| :--- | :--- | :--- | :--- |
| `calculate_gc_content` | $O(L)$ | $O(1)$ | Uses efficient built-in `str.count()`. |
| `calculate_kmer_frequencies` | $O(L)$ | $O(1)$ (fixed alphabet) | Single pass loop using `defaultdict` for counting. Highly efficient. |
| `calculate_hamming_distance`| $O(L)$ | $O(1)$ | Single pass using `zip()`. Optimal for equal-length sequences. |
| **`calculate_distance_matrix`**| $O(N^2 \cdot L)$ | $O(N^2)$ | **Inefficient for large N.** Dominated by $N^2$ comparisons multiplied by the sequence length $L$. |
| **`find_fuzzy_pattern_matches`**| $O(L_T \cdot L_P)$ | $O(1)$ | **Naive Search.** Optimal for small patterns $L_P$; for large $L_P$ or strict efficiency, specialized algorithms like approximate KMP are needed. |
| `get_six_reading_frames` | $O(L)$ | $O(L)$ | Dominated by reverse complement generation (a single pass). |
| `translate_sequence` | $O(L)$ | $O(L)$ | Single pass with efficient dictionary lookups ($O(1)$ per codon). |
| `find_longest_orf` | $O(L)$ | $O(L)$ | Iterates over 6 frames, but each step (translate/extract) is linear, maintaining overall $O(L)$. |

### Code Optimization and Innovation (Q6/Innovation)

1.  **Complexity Awareness:** Algorithms like `calculate_distance_matrix` are noted as $O(N^2 \cdot L)$ to acknowledge their scalability limitation. The innovation lies in providing this honest assessment rather than claiming unrealistic speed.
2.  **Dictionary Mapping:** The `translate_sequence` function relies entirely on the pre-processed `CODON_TO_AMINO_ACID` dictionary, ensuring translation is near-constant time per codon.
3.  **Pythonic Constructs:** Extensive use of `zip()`, `collections.Counter`, and list comprehensions (where memory permits) ensures the best performance possible using the Python standard library.

## Test Coverage and Documentation (Submission Criteria)

All core functions are covered in `cca4_testing_suite.py`.

### Edge Cases Covered:

* **GC Content:** Empty sequences (returns 0.0), sequences with only G/C or only A/T.
* **Hamming Distance:** Unequal sequence lengths (raises `ValueError`).
* **Fuzzy Matching:** Zero mismatches (equivalent to exact match).
* **Translation/ORF:** Incomplete codons at the end of a sequence, presence of multiple start/stop codons, and frames with no identifiable ORF.