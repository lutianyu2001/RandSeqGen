#!/usr/bin/env python3
"""
Random DNA Sequence Generator

This program generates random DNA sequences or combines reference sequences with random bases
to create new sequences. It supports multiprocessing for efficient sequence generation.

Author: Tianyu Lu (tlu83@wisc.edu)
Date: 2025-03-19
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

VERSION = "v0.6.0"
INFO = ("by Tianyu (Sky) Lu (tlu83@wisc.edu)\n"
        "released under Apache-2.0 license")

# PRE-DEFINED PARAMETERS
NUCLEOTIDES = ["A", "T", "G", "C"]
DEFAULT_RANDOM_BASE_RATIO = 0.2
MINIMUM_RANDOM_SEGMENT_NUM = 4
DEFAULT_ALLOCATED_CPU_CORES = os.cpu_count() - 2 if os.cpu_count() > 2 else 1
DEFAULT_OUTPUT_DIR_REL_PATH = "RandSeqGen-Result"

PROGRAM_ROOT_DIR_ABS_PATH = os.path.dirname(__file__)
DEFAULT_OUTPUT_DIR_ABS_PATH = os.path.join(os.getcwd(), DEFAULT_OUTPUT_DIR_REL_PATH)


# ======================================================================================================================


def process_humanized_length(length_str: str) -> int:
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
        ValueError: If the length string format is invalid.

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

    raise ValueError("[ERROR] Invalid length format. Use number, kb or mb (e.g. 100, 1kb, 1mb)")


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
        raise ValueError("[ERROR] Base list cannot be empty.")
    if not lists:
        return sorted(base, key=key, reverse=reverse)

    # Combine all lists and validate length consistency
    all_lists = [base] + list(lists)
    base_len = len(base)
    len_set = set(len(l) for l in all_lists)
    if len(len_set) != 1 and not (len(len_set) == 2 and 0 in len_set):
        raise ValueError("[ERROR] All non-empty lists must have the same length as the base list.")

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


def binary_search_ascending_int_list_latest_not_greater_index(source: List[int], target: int) -> int:
    """
    Find the index of the last element in the source list that is not greater than the target using binary search.

    This function is particularly useful for finding the index of the first reference sequence that 
    does not fit within a length limit in our case.

    Args:
        source (List[int]): List of elements to search
        target (int): Target value to compare against

    Returns:
        int: Index of the last element <= target, or -1 if all elements > target or list is empty

    Precondition:
        The source list is pre-sorted in ascending order
    """
    if not source:
        return -1
    
    left = 0
    right = len(source) - 1
    
    while left <= right:
        mid = (left + right) // 2
        if source[mid] <= target:
            left = mid + 1
        else:
            right = mid - 1
    
    return left - 1


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


def _find_ref_lib_abs_path_list(path: Optional[str]) -> Optional[List[str]]:
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
    if not path:
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

    raise SystemExit(f"[ERROR] Specified reference library not found at {path} or built-in reference library!")


def _load_multiple_ref_libs(path_list: Optional[List[str]], weight_list: Optional[List[float]] = None,
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
    Class for generating random DNA sequences with parallelization.

    This class handles both purely random sequence generation and reference-based
    sequence generation with random base insertion.
    """

    def __init__(self, length: str, number: int, batch: int, processors: int, output_dir: str,
                 ref_lib: Optional[List[str]] = None, ref_lib_weight: Optional[List[float]] = None,
                 num_segment: Optional[int] = None, random_ratio: float = DEFAULT_RANDOM_BASE_RATIO,
                 flag_verbose: bool = False, flag_filter_n: bool = False, flag_track: bool = False):
        self.seq_length: int = process_humanized_length(length)
        self.seq_number: int = number
        self.batch_number: int = batch
        self.processors: int = processors
        self.output_dir_abs_path: str = os.path.abspath(output_dir)

        self.ref_list: Optional[List[str]] = ref_lib
        self.ref_len_list: Optional[List[int]] = None
        self.ref_weight_list: Optional[List[float]] = ref_lib_weight

        # Variables will be reused to store results after initially 
        # being used to pass parameters in _load_multiple_ref_libs()
        self.num_segment: Optional[int] = num_segment
        self.random_ratio: float = random_ratio
        self.random_length = int(self.seq_length * self.random_ratio)

        self.flag_verbose: bool = flag_verbose
        self.flag_filter_n: bool = flag_filter_n
        self.flag_track: bool = flag_track

        self._pre_check()

        # Load reference libraries
        referenced_length = self.seq_length - self.random_length
        ref_lib_data = _load_multiple_ref_libs(ref_lib, ref_lib_weight, referenced_length, flag_filter_n)
        if ref_lib_data:
            self.ref_list, self.ref_len_list, self.ref_weight_list = ref_lib_data
            if self.flag_track:
                print("[INFO] Used reference sequences will be tracked.")
            if not self.num_segment:
                self._choose_num_segment()
        else:
            self.ref_list = None
            self.ref_weight_list = None
            self.num_segment = None
            if ref_lib:
                print("[WARN] No reference sequences found in the specified reference libraries! "
                      "Two possible reasons:\n"
                      "1. Length of sequence does not meet the requirement "
                      "(length of sequence must >= random_ratio * shortest reference sequence's length).\n"
                      "2. Specified reference library does not contain any reference sequences.\n"
                      "Referenced mode has been disabled and program will run in random mode.")

    def _pre_check(self):
        """
        Perform pre-execution checks.

        Raises:
            ValueError: If any strong constraints are violated
        """
        if self.random_ratio < 0 or self.random_ratio > 1:
            raise ValueError("[ERROR] Ratio of random bases must be between 0 and 1!")

        if self.processors < 1:
            raise ValueError("[ERROR] Number of processors must be at least 1!")

        if not self.ref_list:
            print("[INFO] Will execute in random mode.")
            if self.ref_weight_list:
                print("[WARN] weights for each reference library is specified but will be ignored.")
        else:
            print("[INFO] Will execute in referenced mode.")
            len_ref_lib = len(self.ref_list)
            print(f"[INFO] {len_ref_lib} reference libraries are used.")
            if self.ref_weight_list:
                if len_ref_lib != len(self.ref_weight_list):
                    if len_ref_lib == 1:
                        print("[WARN] Multiple weights are specified but only one reference library is specified. "
                              "Weights will be ignored.")
                        self.ref_weight_list = None
                    else:
                        raise ValueError("[ERROR] Number of reference libraries and weights must be the same!")
                if sum(self.ref_weight_list) != 1.0:
                    print("[WARN] The sum of weights for all reference libraries is not 1.0.")

    def _choose_num_segment(self):
        """
        Choose the number of segments for the sequence.
        """
        self.num_segment = self.seq_length // self.ref_len_list[-1]
        len_ref_len_list = len(self.ref_len_list)
        if self.num_segment < MINIMUM_RANDOM_SEGMENT_NUM:
            self.num_segment = self.seq_length // self.ref_len_list[len_ref_len_list * 3 // 4]
        if self.num_segment < MINIMUM_RANDOM_SEGMENT_NUM:
            self.num_segment = self.seq_length // self.ref_len_list[len_ref_len_list // 2]
        if self.num_segment < MINIMUM_RANDOM_SEGMENT_NUM:
            self.num_segment = self.seq_length // self.ref_len_list[len_ref_len_list // 4]
        if self.num_segment < MINIMUM_RANDOM_SEGMENT_NUM:
            self.num_segment = MINIMUM_RANDOM_SEGMENT_NUM

    def _generate_single_referenced_sequence(self, seq_id: str) -> Tuple[SeqRecord, List[SeqRecord]]:
        """
        Generate a sequence combining reference segments with random nucleotides.

        This method:
        1. Calculates random segment distribution
        2. Interleaves random segments with reference sequences
        3. Maintains intact reference sequences
        4. Falls back to random nucleotides when no suitable reference is found

        Returns: 
            Tuple[SeqRecord, List[SeqRecord]]
                - SeqRecord: Generated sequence combining reference and random segments
                - List[SeqRecord]: List of reference sequences
        """
        random_segment_lengths: List[int] = generate_random_segment_lengths(self.num_segment, self.random_length)

        sequence: List[str] = []
        reference_sequences: List[SeqRecord] = []
        current_length: int = 0
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
                suitable_refs_idx = binary_search_ascending_int_list_latest_not_greater_index(self.ref_len_list, 
                                                                                              remaining)

                if suitable_refs_idx > -1:
                    if self.ref_weight_list:
                        ref_seq: str = random.choices(self.ref_list[:suitable_refs_idx + 1],
                                                      weights=self.ref_weight_list[:suitable_refs_idx + 1], k=1)[0]
                    else:
                        ref_seq: str = random.choice(self.ref_list[:suitable_refs_idx + 1])
                    sequence.append(ref_seq)
                    len_ref_seq = len(ref_seq)
                    new_length = current_length + len_ref_seq
                    if self.flag_track:
                        ref_seq_id = f"{seq_id}_{current_length}_{new_length}-+-{len_ref_seq}"
                        reference_sequences.append(create_sequence_record(ref_seq, ref_seq_id))
                    current_length = new_length
                else:
                    sequence.append(generate_random_sequence(remaining))
                    current_length += remaining

        return create_sequence_record("".join(sequence), seq_id), reference_sequences

    def _generate_single_sequence(self, seq_id: str, verbose: bool = False) -> Tuple[SeqRecord, 
                                                                                     Optional[List[SeqRecord]]]:
        """
        Generate a single sequence with appropriate mode.

        Args:
            seq_id (str): ID for the sequence to be generated

        Returns: 
            Tuple[SeqRecord, Optional[List[SeqRecord]]]
                - SeqRecord: Generated sequence as a BioPython SeqRecord object
                - List[SeqRecord] (Optional): List of reference sequences, empty if no reference sequences are used,
                                              or None if executed in random mode
        """

        if not self.ref_list:
            seq: SeqRecord = create_sequence_record(generate_random_sequence(self.seq_length), seq_id)
            return seq, None

        seq, ref_seqs = self._generate_single_referenced_sequence(seq_id)
        return seq, ref_seqs

    def _generate_single_batch(self, batch: int, verbose: bool = False) -> Tuple[List[SeqRecord], 
                                                                                 Optional[List[SeqRecord]]]:
        """
        Generate a batch of sequences using multiprocessing.

        Args:
            batch (int): Batch number to generate

        Returns:
            Tuple[List[SeqRecord], Optional[List[SeqRecord]]]
                - List[SeqRecord]: List of generated sequences for this batch
                - Optional[List[SeqRecord]]: List of reference sequences, empty if no reference sequences are used,
                                              or None if executed in random mode
        """
        print(f"Generating batch {batch}...")

        mp_args_list = [(f"seq_{i}_batch_{batch}_len_{self.seq_length}", verbose) for i in range(self.seq_number)]

        with mp.Pool(self.processors) as pool:
            sequences, ref_seqs_list = zip(*pool.starmap(self._generate_single_sequence, mp_args_list))
        
        sequences: List[SeqRecord] = list(sequences)

        if not self.ref_list or not self.flag_track:
            return sequences, None

        return sequences, list(itertools.chain.from_iterable(ref_seqs_list))

    def execute(self):
        """
        Execute the entire sequence generation process, 
        generates sequences in batches using parallel processing and saves them to separate fasta files
        """
        os.makedirs(self.output_dir_abs_path, exist_ok=True)

        for batch in range(1, self.batch_number + 1):
            sequences, ref_seqs = self._generate_single_batch(batch, self.flag_verbose)
            if not self.ref_list or not self.flag_track or not ref_seqs:
                output_file_path: str = os.path.join(self.output_dir_abs_path, f"sequences_batch_{batch}.fa")
                SeqIO.write(sequences, output_file_path, "fasta")
                print(f"Saved batch {batch}")
            else:
                batch_results_dict: Dict[str, List[SeqRecord]] = {
                    f"batch_{batch}.fa": sequences,
                    f"batch_{batch}_ref.fa": ref_seqs
                }
                save_multi_fasta_from_dict(batch_results_dict, self.output_dir_abs_path)
                print(f"Saved batch {batch} with tracked references")

# ======================================================================================================================

def main():
    """Main function to handle command line arguments and execute the program."""
    parser = argparse.ArgumentParser(prog="RandomSequenceGenerator", 
                                     description="RandSeqGen is a high-performance Python tool for "
                                                 "generating random DNA sequences")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}\n{INFO}")

    parser.add_argument("-l", "--length", help="Length of each sequence (e.g. 100, 1kb, 1mb)", required=True)
    parser.add_argument("-n", "--number", help="Number of sequences in each file (Optional, default: 1)",
                        type=int, default=1)
    parser.add_argument("-b", "--batch", help="Number of generated fasta files (Optional, default: 1)",
                        type=int, default=1)
    parser.add_argument("-p", "--processors", help="Number of processor cores allowed to be used by RandSeqGen "
                                                  f"(Optional, default: {DEFAULT_ALLOCATED_CPU_CORES}). "
                                                  "Please use the number of logical processor cores.",
                        type=int, default=DEFAULT_ALLOCATED_CPU_CORES)
    parser.add_argument("-o", "--output", help="Output directory path "
                                               f"(Optional, default: {DEFAULT_OUTPUT_DIR_ABS_PATH})",
                        type=str, default=DEFAULT_OUTPUT_DIR_ABS_PATH)
    parser.add_argument("-r", "--reference", help="Path to reference library directory or file (Optional)."
                                                  "If not specified, the program will generate complete "
                                                  "random sequence", type=str, nargs="+")
    parser.add_argument("-w", "--weight", help="Weight for each reference library when multiple reference library "
                                                "are specified (Optional). The weights' sum is SUGGESTED to be 1.0 but "
                                                "not required. Will be ignored if only one reference library is used",
                                                type=float, nargs="+")
    parser.add_argument("-s", "--segment", help="Number of randomly split segments for filling in reference sequences "
                                                "(Optional). Program will try to fill every segment with one suitable "
                                                "reference sequence but this is not guaranteed, i.e. the total number "
                                                "of reference sequences used <= number of randomly split segments",
                        type=int, default=None)
    parser.add_argument("--ratio", help="Ratio of random bases in reference mode "
                                        f"(Optional, default: {DEFAULT_RANDOM_BASE_RATIO})",
                        type=float, default=DEFAULT_RANDOM_BASE_RATIO)

    # parser.add_argument("--verbose", help="Verbose mode (Optional). Will show more execution details",
    #                     action="store_true")
    parser.add_argument("--filter_n", help="Avoid using reference sequence with N (Optional)",
                        action="store_true")
    parser.add_argument("--track", help="Track the used reference sequences and store them together with the "
                                        "generated fasta files (Optional). Will not store anything if no reference "
                                        "sequence is used eventually",
                        action="store_true")

    parsed_args = parser.parse_args()

    length = parsed_args.length
    number = parsed_args.number
    batch = parsed_args.batch
    processors = parsed_args.processors
    output_dir_path = parsed_args.output
    ref_lib = parsed_args.reference
    ref_lib_weight = parsed_args.weight
    num_segment = parsed_args.segment
    random_ratio = parsed_args.ratio
    # flag_verbose = parsed_args.verbose
    flag_filter_n = parsed_args.filter_n
    flag_track = parsed_args.track

    SeqGenerator_instance = SeqGenerator(
        length=length, 
        number=number, 
        batch=batch, 
        processors=processors,
        output_dir=output_dir_path, 
        ref_lib=ref_lib, 
        ref_lib_weight=ref_lib_weight,
        num_segment=num_segment, 
        random_ratio=random_ratio,
        flag_filter_n=flag_filter_n, 
        flag_track=flag_track
    )

    SeqGenerator_instance.execute()

if __name__ == "__main__":
    main()
