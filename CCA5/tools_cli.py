import argparse
import sys
import os
import time
import concurrent.futures
from typing import List, Dict, Callable, Any

# Assuming these files are in the same directory for relative import
from cca5_file_handlers import parse_fasta_generator, write_fasta, SequenceDatabase
from cca5_advanced_algs import find_longest_common_subsequence, generate_consensus_sequence

# =================================================================
# PART C: Real-World Applications and Optimization (Q5, Q6)
# =================================================================

# --- Q5: Performance and Scalability Scaffolding ---

def process_chunk(sequences: List[Tuple[str, str]], analysis_func: Callable) -> List[Any]:
    """Helper function to process a chunk of sequences for parallel execution."""
    results = []
    for header, seq in sequences:
        try:
            # Example: Running LCS against a fixed reference sequence (e.g., first sequence in chunk)
            if analysis_func.__name__ == 'find_lcs':
                 results.append((header, analysis_func(seq, sequences[0][1])))
            else:
                 # Placeholder for generic analysis
                 results.append((header, analysis_func(seq)))
        except Exception as e:
            print(f"Error processing {header}: {e}", file=sys.stderr)
            results.append((header, "ERROR"))
    return results

def process_in_parallel(file_path: str, chunk_size: int, analysis_func: Callable) -> List[Any]:
    """
    Q5: Implements parallel processing for sequence analysis using ThreadPoolExecutor.
    Handles memory constraints by processing in chunks.

    :param file_path: Path to the FASTA file.
    :param chunk_size: Number of sequences to process per thread/core.
    :param analysis_func: The function to apply to each sequence/chunk.
    :return: List of all analysis results.
    """
    print(f"Starting parallel processing on {file_path} with chunk size {chunk_size}...")
    start_time = time.time()
    all_results = []
    current_chunk = []

    # Use a high number of workers (e.g., 4) or based on os.cpu_count()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        futures = []
        
        # Use generator to read file, ensuring low memory footprint (Q1, Q5)
        for header, seq in parse_fasta_generator(file_path):
            current_chunk.append((header, seq))
            
            if len(current_chunk) >= chunk_size:
                # Submit the chunk to a worker
                futures.append(executor.submit(process_chunk, current_chunk.copy(), analysis_func))
                current_chunk = []

        # Process the final, potentially smaller chunk
        if current_chunk:
            futures.append(executor.submit(process_chunk, current_chunk, analysis_func))

        # Collect results
        for future in concurrent.futures.as_completed(futures):
            all_results.extend(future.result())

    end_time = time.time()
    print(f"Parallel analysis complete. Time elapsed: {end_time - start_time:.2f} seconds.")
    return all_results

# --- Q6: Command-Line Interface and Integration ---

def main():
    """
    Q6: Creates a command-line interface (CLI) for your tools.
    Implements comprehensive error handling and logging.
    """
    parser = argparse.ArgumentParser(
        description="Bioinformatics Tool Suite (CCA5). Processes FASTA files for advanced analysis.",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('input_file', help="Path to the input FASTA file.")
    
    subparsers = parser.add_subparsers(dest='command', required=True, help='Available commands')

    # Subparser for Consensus Generation (Q3)
    parser_consensus = subparsers.add_parser('consensus', help='Generate consensus sequence from aligned sequences.')
    parser_consensus.add_argument('--output', required=True, help="Output file path for the consensus sequence.")

    # Subparser for Repeat Finding (Q4)
    parser_repeats = subparsers.add_parser('repeats', help='Find common repeats using the Suffix Array concept.')

    # Subparser for Database Management (Q2)
    parser_db = subparsers.add_parser('db_load', help='Load sequences into an in-memory database and perform basic validation.')

    args = parser.parse_args()
    
    # --- Error Handling & Execution ---
    if not os.path.exists(args.input_file):
        print(f"Error: Input file not found at '{args.input_file}'", file=sys.stderr)
        sys.exit(1)

    try:
        if args.command == 'consensus':
            # Load all sequences into memory (required for consensus)
            sequences = [seq for _, seq in parse_fasta_generator(args.input_file)]
            if not sequences: raise ValueError("Input file contained no sequences.")
            
            # --- Use Parallel Processing (Q5) for the calculation ---
            # NOTE: Consensus is usually a single calculation, but we use parallel 
            # to demonstrate the Q5 requirement if sequences were being compared pairwise.
            # Here, we will just use the direct function for simplicity.
            
            consensus, _ = generate_consensus_sequence(sequences)
            
            # Write only the consensus sequence to the output file
            write_fasta({"Consensus_Sequence": consensus}, args.output, line_width=80)
            
        elif args.command == 'repeats':
            # For repeat finding, we process sequence by sequence (only one sequence expected)
            header, seq = next(parse_fasta_generator(args.input_file))
            
            print(f"Analyzing repeats in: {header} (Length: {len(seq)})")
            
            # Run the conceptual suffix array analysis
            repeats = find_repeats_by_suffix_array_concept(seq)
            
            print("\nFound Repeats (Segment, Length):")
            for r, l in repeats:
                print(f"  {r} (Length {l})")

        elif args.command == 'db_load':
            db = SequenceDatabase()
            count = db.load_from_fasta(args.input_file)
            invalid_count = db.validate_data()
            
            print(f"\n--- Database Report ---")
            print(f"Loaded {count} sequences.")
            print(f"Found {invalid_count} sequences with invalid characters.")
            
            # Example search
            results = db.search_by_metadata('seq_1')
            print(f"Search result for 'seq_1': {len(results)} matches.")

    except ValueError as e:
        print(f"Processing Error (Validation): {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        # Q6: Comprehensive error handling
        print(f"An unexpected error occurred: {type(e).__name__}: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    # Usage Example (if run directly):
    # python cca5_tools_cli.py test.fasta consensus --output consensus.fasta
    # python cca5_tools_cli.py genome.fasta repeats
    main()
