# genetic_code.py

"""
Standard Genetic Code Dictionary for Protein Translation.
"""

# Maps a single amino acid (key) to a list of its synonymous codons (value).
STANDARD_CODE = {
    'F': ['TTT', 'TTC'],        # Phenylalanine
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], # Leucine
    'I': ['ATT', 'ATC', 'ATA'], # Isoleucine
    'M': ['ATG'],               # Methionine (START)
    'V': ['GTT', 'GTC', 'GTA', 'GTG'], # Valine
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'], # Serine
    'P': ['CCT', 'CCC', 'CCA', 'CCG'], # Proline
    'T': ['ACT', 'ACC', 'ACA', 'ACG'], # Threonine
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], # Alanine
    'Y': ['TAT', 'TAC'],        # Tyrosine
    'STOP': ['TAA', 'TAG', 'TGA'], # Termination Codons
    'H': ['CAT', 'CAC'],        # Histidine
    'Q': ['CAA', 'CAG'],        # Glutamine
    'N': ['AAT', 'AAC'],        # Asparagine
    'K': ['AAA', 'AAG'],        # Lysine
    'D': ['GAT', 'GAC'],        # Aspartic Acid
    'E': ['GAA', 'GAG'],        # Glutamic Acid
    'C': ['TGT', 'TGC'],        # Cysteine
    'W': ['TGG'],               # Tryptophan
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], # Arginine
    'G': ['GGT', 'GGC', 'GGA', 'GGG'], # Glycine
}

# Reverse lookup dictionary (Codon: Amino Acid) for easy translation
CODON_TO_AMINO_ACID = {}
for aa, codons in STANDARD_CODE.items():
    for codon in codons:
        CODON_TO_AMINO_ACID[codon] = aa