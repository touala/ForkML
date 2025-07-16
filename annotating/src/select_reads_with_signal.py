#!/usr/bin/env python3.10

import sys
import os
import pysam
import tqdm
import numpy as np
import pandas as pd
from itertools import groupby
from concurrent.futures import ThreadPoolExecutor, as_completed
from concurrent.futures import ProcessPoolExecutor
import argparse
import queue

def smooth(y: np.ndarray, x: np.ndarray, sc: int) -> np.ndarray:
    """Smooth a given series."""
    x_round = np.around(x / sc, decimals=0) * sc
    return pd.Series(y).groupby(x_round).mean().to_numpy()

def find1(stro: str, ch: str) -> np.ndarray:
    """Find positions of a character in a string using a buffer."""
    npbuf = np.frombuffer(stro.encode('utf-8'), dtype=np.uint8)
    return np.where(npbuf == ord(ch))[0]

def convert_to_coordinate(seq: str, Ml: np.ndarray, Mm: np.ndarray, mod_type: str) -> np.ndarray: # which: str = "T"
    """Convert sequence, Ml and Mm arrays to coordinate system based on 'T'."""

    assert len(Ml) == len(Mm)
    Mm = Mm + 1
    cum = np.cumsum(Mm) - 1

    result = np.full(len(seq), np.nan)
    if mod_type[0] == "N":
        pos = np.array(range(1, len(seq)))
    else:
        pos = find1(seq, mod_type[0]) # which

    if len(cum) and cum[-1] <= len(pos) - 1:
        Ml = Ml[:len(cum)]
        result[pos[cum]] = np.array(Ml) / 255

    return result

def get_bam_tag(read, possible_tags):
    """Retrieve the first available tag from a list of possible tag names."""
    for tag in possible_tags:
        try:
            return read.get_tag(tag)
        except KeyError:
            continue  # Try the next tag if the current one is missing
    return None  # Return None if no valid tag is found

def parse_modifications(mm_tag: str, ml_values):
    """
    Extract modification positions and corresponding likelihoods from Mm/ML tags.
    
    - `mm_tag`: The Mm/MM tag from the BAM file (modification positions).
    - `ml_values`: The Ml/ML values (modification likelihoods) as a list or array.
    
    Returns:
        Dictionary {modification_type: (positions, likelihoods)}
    """
    modifications = {}
    ml_index = 0  # Tracks position in the likelihood array
    
    for mod_entry in mm_tag.split(";"):
        parts = mod_entry.split(",")
        if len(parts) > 1:
            mod_type = parts[0]  # Modification type (e.g., C, A, G, T)
            positions = np.array(parts[1:], dtype=int)

            # Extract corresponding likelihoods (same length as positions)
            likelihoods = np.array(ml_values[ml_index: ml_index + len(positions)])
            ml_index += len(positions)

            modifications[mod_type] = (positions, likelihoods)

    return modifications

def reverse_complement(seq: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    complement = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(complement)[::-1]

def extract_reference_sequence(reference, read):
    """Retrieve the reference sequence for a given read's alignment."""
    try:
        ref_seq = reference.fetch(read.reference_name, read.reference_start, read.reference_end).upper()
        # Reverse complement if read is on the reverse strand (flag 16)
        if read.flag == 16:
            ref_seq = reverse_complement(ref_seq)
        return ref_seq
    except Exception as e:
        print(f"Warning: Unable to fetch reference for {read.reference_name}:{read.reference_start}-{read.reference_end} ({e})")
        return None

def process_read_with_threshold(read_data):
    """Process a single read's data, apply threshold selection, and return the read_id if it passes the threshold."""
    read_name, seq, Ml_values, Mm_tag, filter_b, n_b, res, threshold, num_points = read_data

    # Extract modification positions and likelihoods
    modifications = parse_modifications(Mm_tag, Ml_values)

    # Process all modification types (C, A, G, T)
    for mod_type, (positions, likelihoods) in modifications.items():
        val = convert_to_coordinate(seq, likelihoods, positions, mod_type)

        if res != 1:
            val = smooth(val, np.arange(len(val)), res).astype(np.float16)

        # Apply threshold selection
        values_above_threshold = np.sum(val > threshold)
        consecutive_high = [sum(1 for _ in g) for k, g in groupby(val > threshold) if k]
        max_consecutive = max(consecutive_high) if consecutive_high else np.nan

        if values_above_threshold > 0 and max_consecutive > num_points:
            return read_name

    return None

def load_reads_bam(bam: str, ref_fasta: str, filter_b: float = 0.5, n_b: int = 0, res: int = 1, threshold: float = 0.1, num_points: int = 10, max_reads: int = None, chs: list = None, tqdm_enabled: bool = False, num_threads: int = 4) -> list:
    """Load reads from BAM file and process them in parallel with threshold selection using ProcessPoolExecutor."""
    valid_read_ids = []

    def generate_read_data():
        """Generator to yield relevant read data to worker processes."""
        with pysam.AlignmentFile(bam, "r") as samfile: # , pysam.FastaFile(ref_fasta) as reference
            for ir, read in tqdm.tqdm(enumerate(samfile), total=max_reads if max_reads else samfile.mapped, disable=not tqdm_enabled):
                if chs and read.reference_name not in chs:
                    continue

                Ml_tag = get_bam_tag(read, ["Ml", "ML"])
                Mm_tag = get_bam_tag(read, ["Mm", "MM"])

                if Ml_tag is None or Mm_tag is None:
                    continue  # Skip if required tags are missing

                if ref_fasta is not None:
                    seq = extract_reference_sequence(reference, read)
                else:
                    seq = read.get_forward_sequence()

                yield (
                    read.query_name,
                    seq,
                    np.array(Ml_tag),  # Convert likelihoods to NumPy array
                    Mm_tag,
                    filter_b, n_b, res, threshold, num_points
                )

                if max_reads and ir >= max_reads - 1:
                    break

    futures = []
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        # Use a queue to manage futures dynamically
        for read_data in generate_read_data():
            # Limit the number of outstanding futures to avoid memory issues
            if len(futures) >= num_threads * 2:
                # Process a batch of completed tasks
                for future in as_completed(futures[:num_threads]):  # Check only the first batch
                    result = future.result()
                    if result:
                        valid_read_ids.append(result)
                futures = futures[num_threads:]  # Keep only the pending tasks

            # Submit new task
            futures.append(executor.submit(process_read_with_threshold, read_data))

        # Process remaining tasks
        for future in as_completed(futures):
            result = future.result()
            if result:
                valid_read_ids.append(result)

    return valid_read_ids

def main():
    parser = argparse.ArgumentParser(description="Select reads with specific signal characteristics from BAM file.")
    parser.add_argument('--in_bam', type=str, required=True, help="Input BAM file path.")
    parser.add_argument('--ref_fasta', type=str, default=None, help="Reference genome FASTA file.")
    parser.add_argument('--res_step', type=float, required=True, help="Resolution step size.")
    parser.add_argument('--thre_level', type=float, required=True, help="Threshold level for signal selection.")
    parser.add_argument('--num_points', type=int, required=True, help="Minimum number of points above threshold.")
    parser.add_argument('--out_folder', type=str, required=True, help="Output folder for results.")
    parser.add_argument('--base_name', type=str, required=True, help="Output base name for results.")
    parser.add_argument('--size_thre', type=int, required=True, help="Minimum length of reads for selection.")
    parser.add_argument('--max_reads', type=int, default=None, help="Maximum number of reads to process.")
    parser.add_argument('--threads', type=int, default=4, help="Number of threads for parallel processing.")
    args = parser.parse_args()

    os.makedirs(args.out_folder, exist_ok=True)

    # Process the BAM file and return valid read IDs
    valid_read_ids = load_reads_bam(
        bam=args.in_bam, 
        ref_fasta=args.ref_fasta,  
        filter_b=0.5,
        n_b=0,
        res=args.res_step, 
        threshold=args.thre_level,
        num_points=args.num_points,
        max_reads=args.max_reads, 
        tqdm_enabled=True, 
        num_threads=args.threads
    )

    # Save results to file
    output_file = os.path.join(args.out_folder, f'listids_thre{args.thre_level}_points_{args.num_points}.{args.base_name}.txt')
    with open(output_file, 'w') as f:
        print(f'Saving {len(valid_read_ids)} read_ids.')
        for read_id in valid_read_ids:
            f.write(f"{read_id}\n")

if __name__ == "__main__":
    main()


