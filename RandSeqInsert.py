#!/usr/bin/env python3
"""
Random DNA Sequence Generator with Insertion

This program takes an input DNA sequence and performs random insertions using sequences
from a reference library. It supports multiprocessing for efficient sequence generation.

Author: Tianyu Lu (tlu83@wisc.edu)
Date: 2024-11-27
"""

import os
import argparse
import itertools
import random
from typing import List, Optional, Union, Callable, Tuple, Dict
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from h5py import ref_dtype

VERSION = "v0.1.0"
INFO = ("by Tianyu (Sky) Lu (tlu83@wisc.edu)\n"
        "released under LGPL-2.1 license")

# PRE-DEFINED PARAMETERS
NUCLEOTIDES = ["A", "T", "G", "C"]
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


def generate_random_segment_lengths(num_segment: int, total_length: int) -> List[int]:
    """
    Generate random segment lengths that sum to total_length.

    Args:
        num_segment (int): Number of segments to generate
        total_length (int): Total length to divide into segments

    Returns:
        List[int]: List of segment lengths that sum to total_length

    Note:
        This is used to determine the distribution of random nucleotides
        throughout the sequence in reference mode.
    """
    if num_segment < 2:
        return [total_length]

    # Generate random points to split the sequence
    split_points: List[int] = sorted(random.sample(range(1, total_length), num_segment - 1))

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
        sorted_ref_len_list (List[int]): List of lengths of reference sequences
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


def save_multi_fasta_from_dict(records_dict: Dict[str, List[SeqRecord]], output_dir_path: str):
    """
    Save multiple FASTA files from a dictionary of sequence records.

    Args:
        records_dict (Dict[str, List[SeqRecord]]): Dictionary mapping file names to lists of SeqRecord objects    
        output_dir_path (str): Directory path where the output files will be saved
    """
    os.makedirs(output_dir_path, exist_ok=True)
    for file_name, records in records_dict.items():
        output_file_path: str = os.path.join(output_dir_path, file_name)
        SeqIO.write(records, output_file_path, "fasta")


def load_sequences(path_list: Optional[List[str]], len_limit: Optional[int] = None,
                   flag_filter_n: bool = False) -> Optional[List[str]]:
    """
    Load reference sequences from reference library files.
    Only loads sequences that meet the length criteria into memory.
    Sequences are stored as strings for efficient processing.

    Args:
        path_list: List of paths to reference sequence files
        len_limit: Optional maximum length limit for sequences to load
        flag_filter_n: Flag to filter out sequences containing N

    Returns:
        Optional[List[str]]: List of reference sequences as strings, or None if no reference paths provided
    """
    if not path_list:
        return None

    ref_sequences: List[str] = []
    for file_path in path_list:
            ref_sequences.extend([str(record.seq.upper()) for record in SeqIO.parse(file_path, "fasta") if
                                  (len(record) <= len_limit if len_limit else True) and
                                  ("N" not in record.seq.upper() if flag_filter_n else True)])

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


def _load_multiple_ref_libs(path_list: List[str], weight_list: Optional[List[float]] = None,
                            len_limit: Optional[int] = None,
                            flag_filter_n: bool = False) -> Optional[Tuple[List[str], List[int], List[float]]]:
    """
    Load multiple reference libraries by combining their sequences.

    Args:
        path_list (List[str]): List of paths to reference libraries
        weight_list (Optional[List[float]]): List of weights for each reference library
        len_limit (Optional[int]): Maximum length limit for sequences to load
        flag_filter_n (bool): Flag to filter out sequences containing N

    Returns:
        Tuple[List[str], List[int], List[float]] (Optional):
            - List[str]: Combined list of reference sequences from all libraries
            - List[int]: List of lengths of reference sequences
            - List[float]: List of weights for each reference sequence
    """
    if not path_list:
        return None

    all_ref_sequence_list: List[str] = []
    all_ref_len_list: List[int] = []
    all_ref_weight_list: List[float] = []

    weight_list: Union[List[float], List[None]] = weight_list or [None] * len(path_list)

    for path_list, weight in zip(path_list, weight_list):
        single_ref_lib_abs_path_list = _find_ref_lib_abs_path_list(path_list)
        if single_ref_lib_abs_path_list:
            single_ref_lib_sequences = load_sequences(single_ref_lib_abs_path_list, len_limit, flag_filter_n)
            if single_ref_lib_sequences:
                all_ref_sequence_list.extend(single_ref_lib_sequences)
                all_ref_len_list.extend([len(seq) for seq in single_ref_lib_sequences])
                if weight:
                    len_single_ref_lib_sequences = len(single_ref_lib_sequences)
                    all_ref_weight_list.extend([weight / len_single_ref_lib_sequences] * len_single_ref_lib_sequences)

    # return sort_multiple_lists(all_ref_sequence_list, all_ref_weight_list, key=len)
    all_ref_len_list, all_ref_sequence_list, all_ref_weight_list = sort_multiple_lists(
        all_ref_len_list, all_ref_sequence_list, all_ref_weight_list)
    return all_ref_sequence_list, all_ref_len_list, all_ref_weight_list


class SeqGenerator:
    """
    Class for generating sequences with random insertions from reference libraries.
    """

    def __init__(self, input_file: str, insertion: int, iteration: int, batch: int, processors: int, output_dir: str,
                 ref_lib: Optional[List[str]] = None, ref_lib_weight: Optional[List[float]] = None,
                 ref_len_limit: Optional[int] = None, flag_filter_n: bool = False, flag_track: bool = False):
        """
        Initialize the sequence generator.

        Args:
            input_file (str): Path to the input sequence file
            insertion (int): Number of insertions per sequence
            iteration (int): Number of times to iterate the insertion process
            batch (int): Number of sequences to generate in each batch
            processors (int): Number of processors to use
            output_dir (str): Directory to save output files
            ref_lib (List[str], optional): List of reference library file paths
            ref_lib_weight (List[float], optional): Weights for reference libraries
            ref_len_limit (int, optional): Maximum length limit for reference sequences to load
            flag_filter_n (bool): Whether to filter out sequences containing N
            flag_track (bool): Whether to track reference sequences used
        """
        self.input_file = input_file
        self.insertion = insertion
        self.iteration = iteration
        self.batch = batch
        self.processors = processors
        self.output_dir = output_dir
        # self.flag_verbose = flag_verbose
        self.ref_len_limit = ref_len_limit
        self.flag_filter_n = flag_filter_n
        self.flag_track = flag_track

        # Load input sequence
        try:
            self.input = list(SeqIO.parse(input_file, "fasta"))
            if not self.input:
                raise SystemExit(f"[ERROR] No sequences found in input file {input_file}")
        except Exception as e:
            raise SystemExit(f"[ERROR] Failed to load input file {input_file}: {str(e)}")

        # Load reference sequences
        if ref_lib:
            ref_lib_data = _load_multiple_ref_libs(ref_lib, ref_lib_weight, ref_len_limit, flag_filter_n)
            if ref_lib_data:
                self.ref_sequences, _, self.ref_weights = ref_lib_data
            else:
                raise SystemExit("ERROR: No valid reference sequences loaded")
        else:
            raise SystemExit("ERROR: Reference library is required for insertion")

    def _pre_check(self):
        """
        Perform pre-execution checks.
        """
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def _perform_single_insertion(self, sequence: str) -> Tuple[str, Optional[SeqRecord]]:
        """
        Perform a single insertion on the given sequence.

        Args:
            sequence (str): The sequence to insert into

        Returns:
            Tuple[str, Optional[SeqRecord]]: The modified sequence and the used reference if tracking
        """
        # Choose random insertion point
        insert_pos = random.randint(0, len(sequence))
        
        # Choose random reference sequence
        ref_seq = random.choices(self.ref_sequences, weights=self.ref_weights, k=1)[0]
        
        # Perform insertion
        new_sequence = sequence[:insert_pos] + ref_seq + sequence[insert_pos:]
        
        if self.flag_track:
            return new_sequence, create_sequence_record(ref_seq, f"ref_used_at_{insert_pos}")
        return new_sequence, None

    def _process_single_sequence(self, seq_record: SeqRecord, seq_num: int) -> Tuple[SeqRecord, Optional[List[SeqRecord]]]:
        """
        Process a single sequence through all iterations.

        Args:
            seq_record (SeqRecord): The input sequence record
            seq_num (int): The sequence number

        Returns:
            Tuple[SeqRecord, Optional[List[SeqRecord]]]: The processed sequence and used references if tracking
        """
        current_seq = str(seq_record.seq)
        used_refs = [] if self.flag_track else None
        
        for i in range(self.iteration):
            for j in range(self.insertion):
                current_seq, used_ref = self._perform_single_insertion(current_seq)
                if self.flag_track and used_ref:
                    used_refs.append(used_ref)
        
        new_id = f"{seq_record.id}_iter{self.iteration}_ins{self.insertion}"
        return create_sequence_record(current_seq, new_id), used_refs

    def _process_batch(self, batch_num: int) -> Tuple[List[SeqRecord], Optional[List[SeqRecord]]]:
        """
        Process a batch of sequences.

        Args:
            batch_num (int): The batch number

        Returns:
            Tuple[List[SeqRecord], Optional[List[SeqRecord]]]: The processed sequences and used references if tracking
        """
        start_idx = batch_num * self.batch
        end_idx = min(start_idx + self.batch, len(self.input))
        
        processed_seqs = []
        all_used_refs = [] if self.flag_track else None
        
        for i, seq_record in enumerate(self.input[start_idx:end_idx], start=start_idx):
            if self.flag_verbose:
                print(f"Processing sequence {i + 1}/{len(self.input)}")
            
            processed_seq, used_refs = self._process_single_sequence(seq_record, i + 1)
            processed_seqs.append(processed_seq)
            
            if self.flag_track and used_refs:
                all_used_refs.extend(used_refs)
        
        return processed_seqs, all_used_refs

    def execute(self):
        """
        Execute the sequence generation process.
        """
        self._pre_check()
        
        total_batches = (len(self.input) + self.batch - 1) // self.batch
        
        if self.processors > 1 and total_batches > 1:
            with mp.Pool(self.processors) as pool:
                results = pool.map(self._process_batch, range(total_batches))
        else:
            results = [self._process_batch(i) for i in range(total_batches)]
        
        # Combine results
        all_sequences = []
        all_refs = [] if self.flag_track else None
        
        for seqs, refs in results:
            all_sequences.extend(seqs)
            if self.flag_track and refs:
                all_refs.extend(refs)
        
        # Save results
        output_dict = {"sequences.fasta": all_sequences}
        if self.flag_track and all_refs:
            output_dict["used_refs.fasta"] = all_refs
        
        save_multi_fasta_from_dict(output_dict, self.output_dir)

def main():
    """Main function to handle command line arguments and execute the program."""
    parser = argparse.ArgumentParser(prog="RandomSequenceInsertion",
                                     description="RandSeqInsert is a high-performance Python tool for "
                                                 "randomly inserting genomic fragments into sequences")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}\n{INFO}")

    # TODO
    # 1. revise help information

    # Required arguments
    parser.add_argument("-i", "--input", help="Input sequence file in FASTA format",
                        type=str, required=True)
    parser.add_argument("-s", "--insertion", help="Number of insertions per sequence",
                        type=int, required=True)
    parser.add_argument("--iteration", type=int, required=True,
                       help="Number of times to iterate the insertion process")
    
    # Optional arguments
    parser.add_argument("-b", "--batch", type=int, default=1,
                       help="Number of sequences to process in each batch")
    parser.add_argument("-p", "--processors", type=int, default=DEFAULT_ALLOCATED_CPU_CORES,
                       help="Number of processors to use")
    parser.add_argument("-o", "--output", default=DEFAULT_OUTPUT_DIR_REL_PATH,
                       help="Output directory path")
    parser.add_argument("-r", "--reference", nargs="+",
                       help="Reference sequence library file paths")
    parser.add_argument("-w", "--weight", type=float, nargs="+",
                       help="Weights for reference libraries")
    parser.add_argument("-l", "--ref_len_limit", help="Weights for reference libraries", type=int,
                        default=None)

    # Optional arguments - flags
    # parser.add_argument("--verbose", action="store_true",
    #                    help="Print verbose output")
    parser.add_argument("--filter-n", action="store_true",
                       help="Filter out sequences containing N")
    parser.add_argument("--track", action="store_true", help="Track and save used reference sequences")

    parsed_args = parser.parse_args()

    input_file = parsed_args.input
    insertion = parsed_args.insertion
    iteration = parsed_args.iteration
    batch = parsed_args.batch
    processors = parsed_args.processors
    output_dir_path = parsed_args.output
    ref_lib = parsed_args.reference
    ref_lib_weight = parsed_args.weight
    ref_len_limit = parsed_args.len_limit
    flag_filter_n = parsed_args.filter_n
    flag_track = parsed_args.track

    generator = SeqGenerator(
        input_file=input_file,
        insertion=insertion,
        iteration=iteration,
        batch=batch,
        processors=processors,
        output_dir=output_dir_path,
        ref_lib=ref_lib,
        ref_lib_weight=ref_lib_weight,
        ref_len_limit=ref_len_limit,
        flag_filter_n=flag_filter_n,
        flag_track=flag_track
    )
    generator.execute()

if __name__ == "__main__":
    main()
