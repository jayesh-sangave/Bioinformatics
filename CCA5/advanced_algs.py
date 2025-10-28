import collections
import itertools
from typing import List, Dict, Tuple

# =================================================================
# PART B: Rosalind-Style Problem Solving (Q3, Q4)
# =================================================================

def find_longest_common_subsequence(seq1: str, seq2: str) -> str:
    """
    Q3: Finds the Longest Common Subsequence (LCS) between two sequences using 
    Dynamic Programming.

    Time Complexity: O(L1 * L2) - where L1 and L2 are sequence lengths.
    Space Complexity: O(L1 * L2) - to store the DP matrix.

    :param seq1: First sequence.
    :param seq2: Second sequence.
    :return: The LCS string.
    """
    s1, s2 = seq1.upper(), seq2.upper()
    L1, L2 = len(s1), len(s2)
    
    # Initialize DP matrix (L1+1 rows, L2+1 columns)
    dp = [[0] * (L2 + 1) for _ in range(L1 + 1)]

    # Fill the DP table
    for i in range(1, L1 + 1):
        for j in range(1, L2 + 1):
            if s1[i - 1] == s2[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] + 1
            else:
                dp[i][j] = max(dp[i - 1][j], dp[i][j - 1])

    # Reconstruct the LCS string by backtracking
    lcs = []
    i, j = L1, L2
    while i > 0 and j > 0:
        if s1[i - 1] == s2[j - 1]:
            lcs.append(s1[i - 1])
            i -= 1
            j -= 1
        elif dp[i - 1][j] > dp[i][j - 1]:
            i -= 1
        else:
            j -= 1

    return "".join(reversed(lcs))


def generate_consensus_sequence(sequences: List[str]) -> Tuple[str, List[Dict[str, int]]]:
    """
    Q3: Generates a consensus sequence (most frequent base at each position) 
    and returns the profile matrix.

    Time Complexity: O(N * L) - N sequences, L length.
    Space Complexity: O(L) - to store the matrix.

    :param sequences: List of aligned sequences (must be same length).
    :return: Tuple of (consensus sequence, profile matrix).
    """
    if not sequences: return "", []
    L = len(sequences[0])
    
    # Validation: Ensure all sequences are of the same length
    if not all(len(s) == L for s in sequences):
        raise ValueError("All sequences must be of the same length for consensus generation.")

    profile = [collections.defaultdict(int) for _ in range(L)]
    consensus = []

    for seq in sequences:
        for i, base in enumerate(seq.upper()):
            profile[i][base] += 1

    for position_counts in profile:
        # Find the base with the maximum count at this position
        best_base = max(position_counts, key=position_counts.get)
        consensus.append(best_base)
        
    return "".join(consensus), profile

def find_repeats_by_suffix_array_concept(sequence: str) -> List[Tuple[str, int]]:
    """
    Q4: Develops a conceptual algorithm for finding repeats using the Suffix Array concept.
    The LCP (Longest Common Prefix) array is implicitly used by comparing adjacent suffixes.
    
    Time Complexity: O(L^2 log L) (due to Python sorting, which is $O(L \log L)$ 
    suffixes of length L) - Actual Suffix Array construction is often O(L log L) or O(L).

    :param sequence: The genomic sequence.
    :return: List of (repeated segment, length of repeat).
    """
    seq = sequence.upper()
    L = len(seq)
    
    # 1. Generate all suffixes
    suffixes = [(seq[i:], i) for i in range(L)] # (suffix, start_index)
    
    # 2. Sort the suffixes (The Suffix Array)
    suffixes.sort()
    
    repeats = collections.defaultdict(int)
    
    # 3. Compare adjacent suffixes to find the Longest Common Prefix (LCP)
    for i in range(L - 1):
        s1, _ = suffixes[i]
        s2, _ = suffixes[i + 1]
        
        # Find the length of the common prefix
        lcp_length = 0
        for b1, b2 in zip(s1, s2):
            if b1 == b2:
                lcp_length += 1
            else:
                break
        
        if lcp_length >= 5: # Only consider repeats of length 5 or greater
            repeat = s1[:lcp_length]
            repeats[repeat] += 1

    # Format output: only include repeats found more than once
    result = [(r, len(r)) for r, count in repeats.items() if count > 0]
    return result

# Note on Phylogeny: Simple evolutionary distance calculation is covered by 
# the distance matrix in CCA4. Phylogenetic relationship analysis requires a separate, 
# complex algorithm (e.g., UPGMA or Neighbor-Joining) which is conceptually addressed
# by the analysis of sequence distances.
