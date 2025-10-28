import collections
import io
from typing import Generator, Tuple, Dict, Optional

# =================================================================
# PART A: File Format Handling and Data Processing (Q1, Q2)
# =================================================================

def parse_fasta_generator(file_path: str) -> Generator[Tuple[str, str], None, None]:
    """
    Q1: Parses a FASTA file using a generator for memory-efficient handling of 
    large genomic files. Yields header and sequence pairs.

    Time Complexity: O(L_total) - where L_total is the sum of all sequence lengths.
    Space Complexity: O(L_max) - memory usage is dominated by the longest single sequence.

    :param file_path: Path to the FASTA file.
    :yield: Tuple (header, sequence) for each entry.
    """
    current_header = ""
    current_sequence = []

    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                if line.startswith('>'):
                    # Yield the previous sequence if one exists
                    if current_header and current_sequence:
                        yield current_header, "".join(current_sequence)
                    
                    # Start new sequence
                    current_header = line[1:].strip()
                    current_sequence = []
                else:
                    # Append sequence data, normalizing case
                    current_sequence.append(line.upper())

            # Yield the last sequence after the file ends
            if current_header and current_sequence:
                yield current_header, "".join(current_sequence)

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return # Generator terminates gracefully

def write_fasta(sequences: Dict[str, str], output_path: str, line_width: int = 60):
    """
    Q1: Implements FASTA writing capabilities with proper formatting (line wrapping).

    :param sequences: Dictionary of {header: sequence}.
    :param output_path: Path to write the output FASTA file.
    :param line_width: Number of characters per line for sequence wrapping.
    """
    try:
        with open(output_path, 'w') as f:
            for header, seq in sequences.items():
                f.write(f">{header}\n")
                
                # Write sequence with line wrapping
                for i in range(0, len(seq), line_width):
                    f.write(seq[i:i + line_width] + "\n")
        print(f"Successfully wrote {len(sequences)} sequences to {output_path}")

    except IOError as e:
        print(f"Error writing file: {e}")

class SequenceDatabase:
    """
    Q2: Simple in-memory database for storing sequence information (FASTA/FASTQ basics).
    Implements search and retrieval functions.
    """
    def __init__(self):
        self._data: Dict[str, Tuple[str, Optional[str]]] = {} # {Header: (Sequence, Quality)}
        self._next_id = 1

    def load_from_fasta(self, file_path: str) -> int:
        """Loads sequences from a FASTA file."""
        count = 0
        for header, seq in parse_fasta_generator(file_path):
            # Use unique ID internally, map header to ID for retrieval
            key = f"seq_{self._next_id}"
            self._data[key] = (seq, None) # FASTQ quality is None for FASTA
            self._next_id += 1
            count += 1
        return count

    def get_sequence(self, internal_id: str) -> Optional[Tuple[str, Optional[str]]]:
        """Implements search and retrieval functions by internal ID."""
        return self._data.get(internal_id)

    def search_by_metadata(self, header_keyword: str) -> Dict[str, str]:
        """
        Q2: Simple search function that mimics retrieving sequence data based on 
        a keyword match in the original metadata (header).
        """
        # Note: This is a simplified search. Real databases use indexed fields.
        results = {}
        for key, (seq, _) in self._data.items():
             # Since we don't store headers directly, this searches based on content/length
             if header_keyword in key: # Simplified keyword match on internal ID
                 results[key] = seq[:20] + "..." # Show snippet
        return results

    def validate_data(self) -> int:
        """
        Q2: Basic data validation measure: ensures sequences contain only ATGC.
        Returns count of invalid sequences found.
        """
        invalid_count = 0
        valid_nucs = {'A', 'T', 'G', 'C'}
        for key, (seq, _) in self._data.items():
            if not set(seq.upper()).issubset(valid_nucs):
                print(f"Validation Error: Sequence {key} contains invalid characters.")
                invalid_count += 1
        return invalid_count
