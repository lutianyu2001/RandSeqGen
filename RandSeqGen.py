#!/usr/bin/env python3
"""
Random DNA Sequence Generator

This program generates random DNA sequences or combines reference sequences with random bases
to create new sequences. It supports multiprocessing for efficient sequence generation.

Author: Tianyu Lu (tlu83@wisc.edu)
Date: 2024-11-27
"""

import os
import argparse
import random
from typing import List, Optional, Union, Callable
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

VERSION = "v0.3.0"
INFO = ("by Tianyu (Sky) Lu (tlu83@wisc.edu)\n"
        "released under LGPL-2.1 license")

# PRE-DEFINED PARAMETERS
NUCLEOTIDES = ["A", "T", "G", "C"]
PROCESS_CORE_RATIO = 1  # Number of processes = CPU cores * PROCESS_CORE_RATIO
DEFAULT_RANDOM_BASE_RATIO = 0.2
DEFAULT_ALLOCATED_CPU_CORES = os.cpu_count() - 2 if os.cpu_count() > 2 else 1
DEFAULT_OUTPUT_DIR_REL_PATH = "RandSeqGen-Result"

PROGRAM_ROOT_DIR_ABS_PATH = os.path.dirname(__file__)
DEFAULT_OUTPUT_DIR_ABS_PATH = os.path.join(os.getcwd(), DEFAULT_OUTPUT_DIR_REL_PATH)


# ======================================================================================================================


def process_length(length_str: str) -> int:
    """
    Process a length string to return the number of bases.

    This function converts various length formats (plain number, kb, mb) into
    the actual number of bases for sequence generation.

    Args:
        length_str (str): A string specifying sequence length. Can be:
            - A plain number (e.g., "100")
            - Kilobases with 'kb' suffix (e.g., "1kb")
            - Megabases with 'mb' suffix (e.g., "1mb")

    Returns:
        int: The number of bases represented by the input string.

    Raises:
        SystemExit: If the length string format is invalid.

    Examples:
        >>> process_length("100")
        100
        >>> process_length("1kb")
        1000
        >>> process_length("1mb")
        1000000
    """
    length_str = length_str.lower()

    # Check for kb/mb suffix
    if length_str.endswith("kb"):
        return int(float(length_str[:-2]) * 1000)
    elif length_str.endswith("mb"):
        return int(float(length_str[:-2]) * 1000000)

    # Check if pure number or b suffix
    if length_str.isdigit() or length_str.endswith("b"):
        return int(length_str.replace('b', ''))

    raise SystemExit("ERROR: Invalid length format. Use number, kb or mb (e.g. 100, 1kb, 1mb)")


def sort_multiple_lists(base: list, *lists: list, 
                        key: Optional[Callable] = None, reverse: bool = False) -> Union[list, tuple]:
    """
    Synchronously sorts multiple lists based on the sorting order of the base list `base`.

    Args:
        base (list): The base list used as the reference for sorting
        *lists (list): Other lists to be sorted synchronously
        key (callable, optional): A function to execute to decide the order. Default is None.
        reverse (bool, optional): Whether to sort in descending order. Default is False.
    
    Returns:
        tuple: A tuple containing all sorted lists. The first element is the sorted base list,
        the rest are the sorted lists in the same order as they were passed in.
    """
    if not base:
        raise ValueError("ERROR: Base list cannot be empty.")
    if not lists:
        return sorted(base, key=key, reverse=reverse)

    # Combine all lists and validate length consistency
    all_lists = [base] + list(lists)
    base_len = len(base)
    len_set = set(len(l) for l in all_lists)
    if len(len_set) != 1 and not (len(len_set) == 2 and 0 in len_set):
        raise ValueError("ERROR: All non-empty lists must have the same length as the base list.")

    # Generate sorting indices
    if key is None:
        sorted_idx = sorted(range(base_len), key=lambda i: base[i], reverse=reverse)
    else:
        sorted_idx = sorted(range(base_len), key=lambda i: key(base[i]), reverse=reverse)

    # Reorganize all lists
    sorted_lists = []
    for l in all_lists:
        if l:
            sorted_l = [l[i] for i in sorted_idx]
            sorted_lists.append(sorted_l)
        else:
            sorted_lists.append([])

    return tuple(sorted_lists)


def create_sequence_record(seq: str, id: str) -> SeqRecord:
    """
    Create a BioPython SeqRecord object for a sequence.

    Args:
        seq (str): The sequence string
        id (str): The id of the sequence

    Returns:
        SeqRecord: BioPython SeqRecord object for the sequence
    """
    return SeqRecord(Seq(seq), id=id, description="")


def generate_random_sequence(length: int) -> str:
    """
    Generate a random DNA sequence of specified length.

    Args:
        length (int): The desired length of the sequence

    Returns:
        str: A random DNA sequence composed of nucleotides from NUCLEOTIDES
    """
    return "".join(random.choices(NUCLEOTIDES, k=length))


def generate_random_segment_lengths(num_segments: int, total_length: int) -> List[int]:
    """
    Generate random segment lengths that sum to total_length.

    Args:
        num_segments (int): Number of segments to generate
        total_length (int): Total length to divide into segments

    Returns:
        List[int]: List of segment lengths that sum to total_length

    Note:
        This is used to determine the distribution of random nucleotides
        throughout the sequence in reference mode.
    """
    if num_segments <= 1:
        return [total_length]

    # Generate random points to split the sequence
    split_points: List[int] = sorted(random.sample(range(1, total_length), num_segments - 1))

    # Calculate segment lengths
    segments: List[int] = [split_points[0]]
    for i in range(1, len(split_points)):
        segments.append(split_points[i] - split_points[i - 1])
    segments.append(total_length - split_points[-1])

    return segments


def binary_search_suitable_refs_index(sorted_ref_len_list: List[int], len_limit: int) -> int:
    """
    Find the index of the first reference sequence that does not fit within length limit using binary search.

    Args:
        sorted_ref_sequences (List[str]): Pre-sorted list of reference sequences
        len_limit (int): Maximum length of reference sequences to find

    Returns:
        int: Index of the first reference sequence whose length > len_limit

    Note:
        Uses binary search because reference sequences are pre-sorted by length.
    """
    left = 0
    right = len(sorted_ref_len_list) - 1
    while left <= right:
        mid = (left + right) // 2
        if sorted_ref_len_list[mid] <= len_limit:
            left = mid + 1
        else:
            right = mid - 1

    return left


def save_single_batch(records: List[SeqRecord], batch: int, output_dir_path: str):
    """
    Save a single batch of sequences to a FASTA file.

    Args:
        records (List[SeqRecord]): List of BioPython SeqRecord objects to save
        batch (int): Current batch number
        output_dir_path (str): Directory path where the output files will be saved

    Note:
        Files are saved in FASTA format with names like 'sequences_batch_1.fa'
    """
    output_file = os.path.join(output_dir_path, f"sequences_batch_{batch}.fa")
    SeqIO.write(records, output_file, "fasta")
    print(f"Generated batch {batch}")


def load_sequences(path_list: Optional[List[str]], len_limit: Optional[int] = None) -> Optional[List[str]]:
    """
    Load reference sequences from reference library files.
    Only loads sequences that meet the length criteria into memory.
    Sequences are stored as strings for efficient processing.

    Args:
        path_list: List of paths to reference sequence files
        len_limit: Optional maximum length limit for sequences to load

    Returns:
        Optional[List[str]]: List of reference sequences as strings, or None if no reference paths provided
    """
    if not path_list:
        return None

    ref_sequences: List[str] = []
    for file_path in path_list:
        if len_limit is not None:
            ref_sequences.extend([str(record.seq.upper()) for record in SeqIO.parse(file_path, "fasta") if
                                  len(record) <= len_limit])
        else:
            ref_sequences.extend([str(record.seq.upper()) for record in SeqIO.parse(file_path, "fasta")])

    ref_sequences.sort(key=len)
    return ref_sequences


def _find_ref_lib_abs_path_list(path: str) -> Optional[List[str]]:
    """
    Locate the reference library files.

    This function attempts to find the given reference library directory or file in two locations:
    1. The directly specified path
    2. Under the default "lib" directory in the program root directory

    Args:
        path (str): The path to search for reference library. Can be:
                        - An absolute path to a directory or file
                        - A relative path to a directory or file

    Returns:
        Optional[List[str]]: A list containing the absolute paths to the found reference library files,
                             or None if path is None.

    Raises:
        SystemExit: If no valid reference library directory or file is found.

    Examples:
        >>> _find_ref_lib_abs_path_list("/abs/path/to/refs.fa")
        ["/abs/path/to/refs.fa"]
        >>> _find_ref_lib_abs_path_list("lib/rice")
        ["/path/to/program/lib/rice/rice_DTA.fa", "/path/to/program/lib/rice/rice_DTC.fa",
         "/path/to/program/lib/rice/rice_DTH.fa", "/path/to/program/lib/rice/rice_DTM.fa",
         "/path/to/program/lib/rice/rice_DTT.fa"]
    """
    if path is None:
        return None

    # Check if path exists and is a file
    if os.path.exists(path) and os.path.isfile(path):
        return [os.path.abspath(path)]

    # Check if path exists and is a directory
    if os.path.exists(path) and os.path.isdir(path):
        ref_lib_files: List[str] = []
        for file in os.listdir(path):
            if file.endswith((".fa", ".fasta")):
                ref_lib_files.append(os.path.join(path, file))
        if ref_lib_files:
            return list(map(os.path.abspath, ref_lib_files))

    # Check in built-in RefLib directory
    default_ref_lib_abs_path = os.path.join(PROGRAM_ROOT_DIR_ABS_PATH, path)
    if os.path.exists(default_ref_lib_abs_path):
        if os.path.isfile(default_ref_lib_abs_path):
            return [default_ref_lib_abs_path]
        elif os.path.isdir(default_ref_lib_abs_path):
            ref_lib_files: List[str] = []
            for file in os.listdir(default_ref_lib_abs_path):
                if file.endswith((".fa", ".fasta")):
                    ref_lib_files.append(os.path.join(default_ref_lib_abs_path, file))
            if ref_lib_files:
                return ref_lib_files

    raise SystemExit(f"ERROR: Specified reference library not found at {path} or built-in reference library!")


def _load_multiple_ref_libs(ref_lib_path_list: List[str], ref_lib_weight_list: Optional[List[float]] = None,
                            len_limit: Optional[int] = None) -> Optional[List[str]]:
    """
    Load multiple reference libraries by combining their sequences.

    Args:
        ref_lib_path_list (List[str]): List of paths to reference libraries
        ref_lib_weight_list (Optional[List[float]]): List of weights for each reference library
        len_limit (Optional[int]): Maximum length limit for sequences to load

    Returns:
        Optional[List[str]]: Combined list of reference sequences from all libraries,
                             or None if no paths are provided.
    """
    if not ref_lib_path_list:
        return None

    all_ref_sequences: List[str] = []
    all_ref_weights: List[float] = []
    if ref_lib_weight_list:
        for path_list, weight in zip(ref_lib_path_list, ref_lib_weight_list):
            single_ref_lib_abs_path_list = _find_ref_lib_abs_path_list(path_list)
            if single_ref_lib_abs_path_list:
                single_ref_lib_sequences = load_sequences(single_ref_lib_abs_path_list, len_limit)
                if single_ref_lib_sequences:
                    all_ref_sequences.extend(single_ref_lib_sequences)

                    all_ref_weights.extend([weight] * len(single_ref_lib_sequences))
    else:
        for path_list in ref_lib_path_list:
            single_ref_lib_abs_path_list = _find_ref_lib_abs_path_list(path_list)
            if single_ref_lib_abs_path_list:
                single_ref_lib_sequences = load_sequences(single_ref_lib_abs_path_list, len_limit)
                if single_ref_lib_sequences:
                    all_ref_sequences.extend(single_ref_lib_sequences)

    return sort_multiple_lists(all_ref_sequences, all_ref_weights, key=len)


class SeqGenerator:
    """
    Class for generating random DNA sequences with parallelization.

    This class handles both purely random sequence generation and reference-based
    sequence generation with random base insertion.
    """

    def __init__(self, length: str, number: int, batch: int, processors: int, output_dir: str,
                 ref_lib: Optional[List[str]] = None, ref_lib_weight: Optional[List[float]] = None,
                 random_ratio: float = DEFAULT_RANDOM_BASE_RATIO, flag_verbose: bool = False):
        self.seq_length: int = process_length(length)
        self.seq_number: int = number
        self.batch_number: int = batch
        self.processors: int = processors
        self.output_dir_abs_path: str = os.path.abspath(output_dir)
        self.ref_list: Optional[List[str]] = ref_lib
        self.ref_len_list: Optional[List[int]] = None
        self.ref_weight_list: Optional[List[float]] = ref_lib_weight
        self.random_ratio: float = random_ratio
        self.flag_verbose: bool = flag_verbose

        self._pre_check()

        # self.ref_sequences: List[str] = load_reference_sequences(self.ref_lib_abs_path,
        #                                                          self.seq_length // self.random_ratio)
        self.ref_list, self.ref_weight_list = _load_multiple_ref_libs(ref_lib, ref_lib_weight,
                                                                      self.seq_length // self.random_ratio)

        if ref_lib and not self.ref_list:
            print("WARN: No reference sequences found in the specified reference libraries! "
                  "Two possible reasons:\n"
                  "1. Length of sequence does not meet the requirement "
                  "(length of sequence must >= random_ratio * shortest reference sequence's length).\n"
                  "2. Specified reference library does not contain any reference sequences.\n"
                  "Referenced mode has been disabled and program will run in random mode.")
            self.ref_list = None
        
        if self.ref_list:
            self.ref_len_list = list(map(len, self.ref_list))

    def _pre_check(self):
        """
        Perform pre-execution checks.

        Raises:
            SystemExit: If any validation checks fail
        """
        if self.random_ratio < 0 or self.random_ratio > 1:
            raise ValueError("ERROR: Ratio of random bases must be between 0 and 1!")

        if self.processors < 1:
            raise ValueError("ERROR: Number of processors must be at least 1!")

        if not self.ref_list:
            print("INFO: Random mode execution.")
            if self.ref_weight_list:
                print("WARN: weights for each reference library is specified but will be ignored.")
        else:
            print("INFO: Referenced mode execution.")
            len_ref_lib = len(self.ref_list)
            print(f"INFO: {len_ref_lib} reference libraries are used.")
            if self.ref_weight_list:
                if len_ref_lib != len(self.ref_weight_list):
                    if len_ref_lib == 1:
                        print("WARN: Multiple weights are specified but only one reference library is specified. "
                              "Weights will be ignored.")
                        self.ref_weight_list = None
                    else:
                        raise ValueError("ERROR: Number of reference libraries and weights must be the same!")
                if sum(self.ref_weight_list) != 1.0:
                    print("WARN: The sum of weights for all reference libraries is not 1.0.")

    def _generate_single_referenced_sequence(self) -> str:
        """
        Generate a sequence combining reference segments with random nucleotides.

        This method:
        1. Calculates random segment distribution
        2. Interleaves random segments with reference sequences
        3. Maintains intact reference sequences
        4. Falls back to random nucleotides when no suitable reference is found

        Returns:
            str: Generated sequence combining reference and random segments
        """
        random_length = int(self.seq_length * self.random_ratio)
        num_random_segments = max(2, int(self.seq_length / 1000))

        random_segment_lengths = generate_random_segment_lengths(num_random_segments, random_length)
        # ref_length = self.seq_length - random_length

        sequence: List[str] = []
        current_length = 0
        random_idx: int = 0

        while current_length < self.seq_length:
            # Add random segment
            if random_idx < len(random_segment_lengths):
                sequence.append(generate_random_sequence(random_segment_lengths[random_idx]))
                current_length += random_segment_lengths[random_idx]
                random_idx += 1

            # Try to add reference segment
            if current_length < self.seq_length:
                remaining = self.seq_length - current_length
                suitable_refs_idx = binary_search_suitable_refs_index(self.ref_len_list, remaining)

                if suitable_refs_idx > 0:
                    ref_seq = random.choices(self.ref_list[:suitable_refs_idx], 
                                             weights=self.ref_weight_list[:suitable_refs_idx], k=1)[0]
                    sequence.append(ref_seq)
                    current_length += len(ref_seq)
                else:
                    sequence.append(generate_random_sequence(remaining))
                    current_length += remaining

        return "".join(sequence)

    def _generate_single_sequence(self, seq_id: str, verbose: bool = False) -> SeqRecord:
        """
        Generate a single sequence with appropriate mode.

        Args:
            seq_id (str): ID for the sequence to be generated

        Returns:
            SeqRecord: Generated sequence as a BioPython SeqRecord object
        """

        if self.ref_list:
            seq: str = self._generate_single_referenced_sequence()
        else:
            seq: str = generate_random_sequence(self.seq_length)

        if verbose:
            print(f"  Generated: {seq_id}")

        return create_sequence_record(seq, seq_id)

    def _generate_batch(self, batch: int, verbose: bool = False) -> List[SeqRecord]:
        """
        Generate a batch of sequences using multiprocessing.

        Args:
            batch (int): Batch number to generate

        Returns:
            List[SeqRecord]: List of generated sequences for this batch
        """
        print(f"\nGenerating batch {batch}...")

        mp_args_list = [(f"seq_{i}_batch_{batch}_len_{self.seq_length}", verbose) for i in range(self.seq_number)]

        with mp.Pool(self.processors) as pool:
            sequences = pool.starmap(self._generate_single_sequence, mp_args_list)

        return sequences

    def execute(self):
        """
        Execute the entire sequence generation process.

        This method:
        1. Sets up multiprocessing pool based on CPU cores
        2. Generates sequences in batches using parallel processing
        3. Saves each batch to a separate FASTA file
        4. Handles both random and reference-based sequence generation
        """
        batch_results: List[List[SeqRecord]] = []
        for batch in range(self.batch_number):
            batch_results.append(self._generate_batch(batch, self.flag_verbose))

        print()
        os.makedirs(self.output_dir_abs_path, exist_ok=True)

        with mp.Pool(self.processors) as pool:
            mp_args_list = [(records, batch, self.output_dir_abs_path) for batch, records in enumerate(batch_results)]
            pool.starmap(save_single_batch, mp_args_list)


# ======================================================================================================================


if __name__ == "__main__":
    """
    Main entry point for the sequence generator.
    """
    parser = argparse.ArgumentParser(prog="RandomSequenceGenerator")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}\n{INFO}")

    parser.add_argument("-l", "--length", help="Length of each sequence (e.g. 100, 1kb, 1mb)", required=True)
    parser.add_argument("-n", "--number", help="Number of sequences in each file (Optional, default: 1)",
                        type=int, default=1)
    parser.add_argument("-b", "--batch", help="Number of generated fasta files (Optional, default: 1)",
                        type=int, default=1)
    parser.add_argument("-p", "--processor", help="Number of processor cores allowed to be used by RandSeqGen "
                                                  f"(Optional, default: {DEFAULT_ALLOCATED_CPU_CORES}). "
                                                  "Please use the number of logical processor cores.",
                        type=int, default=DEFAULT_ALLOCATED_CPU_CORES)
    parser.add_argument("-o", "--output", help="Output directory path "
                                               f"(Optional, default: {DEFAULT_OUTPUT_DIR_ABS_PATH})",
                        default=DEFAULT_OUTPUT_DIR_ABS_PATH)

    # parser.add_argument("-r", "--reference", help="Path to reference library directory (Optional)."
    #                                               "If not specified, the program will generate complete "
    #                                               "random sequence.", default=None)
    
    parser.add_argument("-r", "--reference", help="Path to reference library directory / file (Optional)."
                                                  "If not specified, the program will generate complete "
                                                  "random sequence.", action="append", default=[])
    parser.add_argument("-w", "--weight", help="Weight for each reference library when multiple reference library "
                                                "are specified (Optional). The weights' sum is SUGGESTED to be 1.0 but "
                                                "not required. Will be ignored if only one reference library is used.",
                                                action="append", default=[])

    parser.add_argument("--ratio", help="Ratio of random bases in reference mode "
                                        f"(Optional, default: {DEFAULT_RANDOM_BASE_RATIO})",
                        type=float, default=DEFAULT_RANDOM_BASE_RATIO)

    parser.add_argument("--verbose", help="Verbose mode (Optional). Will show more execution details.",
                        action="store_true")

    parsed_args = parser.parse_args()

    length = parsed_args.length
    number = parsed_args.number
    batch = parsed_args.batch
    processors = parsed_args.processor
    output_dir_path = parsed_args.output
    ref_lib = parsed_args.reference
    ref_lib_weight = parsed_args.weight
    random_ratio = parsed_args.ratio
    flag_verbose = parsed_args.verbose

    SeqGenerator_instance = SeqGenerator(length=length, number=number, batch=batch, processors=processors,
                                         output_dir=output_dir_path, ref_lib=ref_lib, 
                                         ref_lib_weight=ref_lib_weight,
                                         random_ratio=random_ratio, flag_verbose=flag_verbose)

    SeqGenerator_instance.execute()
