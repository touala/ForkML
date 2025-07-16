#!/usr/bin/env python3

import pysam
import numpy as np
from tqdm import tqdm
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed


def get_bam_tag(read, possible_tags):
    """Retrieve the first available tag from a list of possible tag names."""
    for tag in possible_tags:
        try:
            return read.get_tag(tag)
        except KeyError:
            continue  # Try the next tag if the current one is missing
    return None  # Return None if no valid tag is found

def parse_modifications(mm_tag: str, ml_values):
    """Extract modification positions and likelihoods from Mm/ML tags."""
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

def process_read(read_data):
    read_id, reference_name, reference_start, ml_tag, mm_tag = read_data

    if ml_tag is None or mm_tag is None:
        return None  # Skip if required tags are missing

    # Convert to numpy array
    ml_values = np.array(ml_tag)

    # Extract modification positions and likelihoods
    modifications = parse_modifications(mm_tag, ml_values)

    results = []
    for mod_type, (positions, likelihoods) in modifications.items():
        likelihoods = np.array(likelihoods) / 255.0  # Normalize likelihoods (0 to 1)
        mod_score = np.mean(np.abs(likelihoods - 0.5))  # Compute score
        ml_length = len(likelihoods)  # Get number of modification likelihoods

        results.append((read_id, mod_type, mod_score, reference_name, reference_start, ml_length))

    return results  # Return list of tuples

def extract_ml_tag_mean(bam, num_threads=4, max_reads=None, tqdm_enabled=True):
    results = []

    def generate_read_data():
        """Generator that extracts only needed attributes from BAM reads."""
        with pysam.AlignmentFile(bam, "r") as samfile:
            for ir, read in tqdm(enumerate(samfile), total=max_reads if max_reads else samfile.mapped, disable=not tqdm_enabled):
                ml_tag = get_bam_tag(read, ["Ml", "ML"])
                mm_tag = get_bam_tag(read, ["Mm", "MM"])

                if ml_tag is None or mm_tag is None:
                    continue  # Skip if missing tags

                yield (
                    read.query_name,
                    read.reference_name,
                    read.reference_start,
                    np.array(ml_tag),  # Convert likelihoods to NumPy array
                    mm_tag
                )

                if max_reads and ir >= max_reads - 1:
                    break

    futures = []
    with ProcessPoolExecutor(max_workers=num_threads) as executor:
        for read_data in generate_read_data():
            if len(futures) >= num_threads * 2:
                for future in as_completed(futures[:num_threads]):  # Process only first batch
                    result = future.result()
                    if result:
                        results.extend(result)  # Extend with list of tuples
                futures = futures[num_threads:]  # Keep pending tasks

            futures.append(executor.submit(process_read, read_data))

        # Process remaining tasks
        for future in as_completed(futures):
            result = future.result()
            if result:
                results.extend(result)

    return results

def write_results_to_file(tag_results, output_file):
    sorted_results = sorted(tag_results, key=lambda x: x[2])  # Sort by score

    with open(output_file, 'w') as f:
        f.write("read_id\tmod_type\tsignal_score\tchromosome\tmap_start\ttag_length\n")
        for read_id, mod_type, signal_score, chromosome, map_start, tag_length in sorted_results:
            f.write(f"{read_id}\t{mod_type}\t{signal_score:.6f}\t{chromosome}\t{map_start}\t{tag_length}\n")

    print(f"Results written to {output_file}")

def main():
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Extract and compute the mean of 'Ml' or 'ML' tags from a BAM file")
    parser.add_argument("bam_file", type=str, help="Path to the BAM file")
    parser.add_argument("output_file", type=str, help="Path to the output file where results will be saved")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing (default: 4)")
    parser.add_argument("--max_reads", type=int, default=None, help="Number of reads to process (default: All)")
    
    # Parse command line arguments
    args = parser.parse_args()
    bam_file_path = args.bam_file
    output_file = args.output_file
    num_threads = args.threads
    max_reads = args.max_reads

    # Process the BAM file and compute ML tag scores
    tag_scores = extract_ml_tag_mean(bam_file_path, num_threads=num_threads, max_reads=max_reads)

    # Write the results to the output file
    write_results_to_file(tag_scores, output_file)

if __name__ == "__main__":
    main()
