#!/usr/bin/env python3
"""
Random Sequence Insertion Generator

This program takes an input sequence and performs random insertions using sequences
from a reference library. It supports multiprocessing for efficient sequence generation.

Author: Tianyu Lu (tlu83@wisc.edu)
Date: 2024-11-27
"""

import os
import argparse
import random
from typing import List, Optional, Union, Callable, Tuple, Dict
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

VERSION = "v1.0.0"
INFO = ("by Tianyu (Sky) Lu (tianyu@lu.fm)\n"
        "For TEGE Research")

PROGRAM_ROOT_DIR_ABS_PATH = os.path.dirname(__file__)

# PRE-DEFINED PARAMETERS
DEFAULT_TSD_SNP_MUTATION_RATE = 0.05
DEFAULT_TSD_INDEL_MUTATION_RATE = 0.05

TREE_VISUAL_DIR_NAME = "tree_visual"
DEFAULT_ALLOCATED_CPU_CORES = os.cpu_count() - 2 if os.cpu_count() > 2 else 1

DEFAULT_OUTPUT_DIR_REL_PATH = "RandSeqInsert-Result"
DEFAULT_OUTPUT_DIR_ABS_PATH = os.path.join(os.getcwd(), DEFAULT_OUTPUT_DIR_REL_PATH)


# ======================================================================================================================
# TSD (Target Site Duplication) Functions


def add_snp_mutation(seq: str) -> str:
    """Add a single SNP mutation at a random position in the sequence.

    Args:
        seq (str): Input DNA sequence

    Returns:
        str: Mutated sequence (maintaining the same length)
    """
    if not seq:
        return seq

    # Randomly choose mutation position
    pos = random.randrange(len(seq))
    curr_base = seq[pos]
    # Choose a different base
    new_base = random.choice(list(set("ATGC") - {curr_base}))

    # Construct mutated sequence
    return seq[:pos] + new_base + seq[pos+1:]

def add_indel_mutation(seq: str) -> str:
    """Add a single indel mutation at a random position in the sequence.

    Args:
        seq (str): Input DNA sequence

    Returns:
        str: Mutated sequence (maintaining the same length)
    """
    if len(seq) < 2:
        return seq

    # 50% probability of insertion or deletion
    is_insertion = random.choice([True, False])

    # Choose mutation position (avoid the last position)
    pos = random.randrange(len(seq)-1)
    new_base = random.choice("ATGC")

    if is_insertion:
        # Insert new base and remove last base to maintain length
        return seq[:pos] + new_base + seq[pos:-1]
    else:
        # Delete base and append random base at the end
        return seq[:pos] + seq[pos+1:] + new_base


def generate_TSD(seq_slice: str,
                 length: int = float("inf"),
                 snp: float = DEFAULT_TSD_SNP_MUTATION_RATE,
                 indel: float = DEFAULT_TSD_INDEL_MUTATION_RATE) -> Tuple[str, str]:
    """Generate TSD sequences with potential mutations.

    Args:
        seq_slice (str): Target sequence fragment
        length (int): TSD length (defaults to seq_slice length)
        snp (float): SNP mutation probability
        indel (float): Indel mutation probability

    Returns:
        Tuple[str, str]: (5' TSD sequence, 3' TSD sequence)
    """
    # Initialize both ends of TSD
    tsd_length = min(len(seq_slice), length)
    tsd_5 = tsd_3 = seq_slice[:tsd_length]

    # Apply SNP mutation (randomly choose end)
    if random.random() < snp:
        if random.choice([True, False]):
            tsd_5 = add_snp_mutation(tsd_5)
        else:
            tsd_3 = add_snp_mutation(tsd_3)

    # Apply Indel mutation (randomly choose end)
    if random.random() < indel:
        if random.choice([True, False]):
            tsd_5 = add_indel_mutation(tsd_5)
        else:
            tsd_3 = add_indel_mutation(tsd_3)

    return tsd_5, tsd_3


# ======================================================================================================================


class SequenceNode:
    """
    Node in a sequence tree structure, used for efficient sequence insertion operations.
    Implemented as an AVL tree to maintain balance during insertions.
    """
    def __init__(self, content: str, is_reference: bool = False, metadata: dict = None):
        """
        Initialize a sequence node.

        Args:
            content (str): The sequence content
            is_reference (bool): Whether this node contains a reference sequence
            metadata: Additional information about the node (such as position)
        """
        self.content = content
        self.length = len(content)
        self.is_reference = is_reference
        self.metadata = metadata
        self.left = None
        self.right = None
        # Total length of the subtree for efficient traversal
        self.total_length = self.length
        # Height of the node for AVL balancing
        self.height = 1

    def __str__(self) -> str:
        """
        Convert the tree to a string by in-order traversal.

        Returns:
            str: The concatenated sequence
        """
        return "".join(self.collect_content_in_order_traversal())

    def to_graphviz(self, node_id_prefix: str = "node", abs_pos: int = 0) -> str:
        """
        Generate Graphviz DOT format string for visualizing the tree structure.
        
        Args:
            node_id_prefix (str): Prefix for node IDs in the graph
            abs_pos (int): Current absolute position in the concatenated sequence
            
        Returns:
            str: Graphviz DOT format string
        """
        # Initialize the DOT string with graph declaration
        dot_str = ['digraph SequenceTree {',
                   '  node [shape=box, style=filled];']

        # Generate nodes and edges through recursive traversal
        nodes, edges = self.__build_graphviz_nodes_edges(node_id_prefix, abs_pos)
        
        # Add all nodes and edges to the DOT string
        for node in nodes:
            dot_str.append(f'  {node}')
        for edge in edges:
            dot_str.append(f'  {edge}')
        
        dot_str.append('}')
        return '\n'.join(dot_str)
    
    def __build_graphviz_nodes_edges(self, node_id_prefix: str, abs_pos: int) -> tuple:
        """
        Recursively build nodes and edges for Graphviz visualization.
        
        Args:
            node_id_prefix (str): Prefix for node IDs
            abs_pos (int): Current absolute position in the concatenated sequence
            
        Returns:
            tuple: (nodes list, edges list)
        """
        nodes = []
        edges = []
        
        # Calculate positions
        left_length = self.left.total_length if self.left else 0
        
        # Calculate start and end positions
        start_pos = abs_pos + left_length
        end_pos = start_pos + self.length
        
        # Generate unique ID for this node
        node_id = f"{node_id_prefix}_{start_pos}_{end_pos}"
        
        # Determine node type and color
        node_type = "REF" if self.is_reference else "SRC"
        fill_color = "lightblue" if self.is_reference else "lightgreen"
        
        # Create node label with position information
        label = f"{node_type}\\nStart: {start_pos}\\nEnd: {end_pos}\\nLength: {self.length}"
        
        # Add the node to the nodes list
        nodes.append(f'{node_id} [label="{label}", fillcolor="{fill_color}"];')
        
        # Process left child if exists
        if self.left:
            left_nodes, left_edges = self.left.__build_graphviz_nodes_edges(f"{node_id_prefix}_L", abs_pos)
            nodes.extend(left_nodes)
            edges.extend(left_edges)
            
            # Find the ID of the left child's root node
            left_abs_pos = abs_pos
            left_start = left_abs_pos + (self.left.left.total_length if self.left.left else 0)
            left_end = left_start + self.left.length
            left_id = f"{node_id_prefix}_L_{left_start}_{left_end}"
            
            # Add edge from this node to left child
            edges.append(f'{node_id} -> {left_id} [label="L"];')
        
        # Process right child if exists
        if self.right:
            right_abs_pos = abs_pos + left_length + self.length
            right_nodes, right_edges = self.right.__build_graphviz_nodes_edges(f"{node_id_prefix}_R", right_abs_pos)
            nodes.extend(right_nodes)
            edges.extend(right_edges)
            
            # Find the ID of the right child's root node
            right_start = right_abs_pos + (self.right.left.total_length if self.right.left else 0)
            right_end = right_start + self.right.length
            right_id = f"{node_id_prefix}_R_{right_start}_{right_end}"
            
            # Add edge from this node to right child
            edges.append(f'{node_id} -> {right_id} [label="R"];')
        
        return nodes, edges

    def __update_total_length(self):
        """
        Update the total length of this subtree.

        This method recalculates the total length of the subtree rooted at this node
        by summing the lengths of the left subtree, the current node, and the right subtree.
        """
        left_length = self.left.total_length if self.left else 0
        right_length = self.right.total_length if self.right else 0
        self.total_length = left_length + self.length + right_length

    def __update_height(self):
        """
        Update the height of this node.
        """
        left_height = self.left.height if self.left else 0
        right_height = self.right.height if self.right else 0
        self.height = max(left_height, right_height) + 1

    def __get_balance(self):
        """
        Calculate the balance factor of this node.

        Returns:
            int: Balance factor (left height - right height)
        """
        left_height = self.left.height if self.left else 0
        right_height = self.right.height if self.right else 0
        return left_height - right_height

    def __rotate_right(self):
        """
        Perform a right rotation on this node.

        Returns:
            SequenceNode: The new root after rotation
        """
        # Store the left child as the new root
        new_root = self.left
        
        # The left child's right subtree becomes this node's left subtree
        self.left = new_root.right
        
        # This node becomes the new root's right child
        new_root.right = self
        
        # Update heights and total lengths
        self.__update_height()
        self.__update_total_length()
        new_root.__update_height()
        new_root.__update_total_length()
        
        return new_root

    def __rotate_left(self):
        """
        Perform a left rotation on this node.

        Returns:
            SequenceNode: The new root after rotation
        """
        # Store the right child as the new root
        new_root = self.right
        
        # The right child's left subtree becomes this node's right subtree
        self.right = new_root.left
        
        # This node becomes the new root's left child
        new_root.left = self
        
        # Update heights and total lengths
        self.__update_height()
        self.__update_total_length()
        new_root.__update_height()
        new_root.__update_total_length()
        
        return new_root

    def __balance(self):
        """
        Balance this node if needed.

        Returns:
            SequenceNode: The new root after balancing
        """
        # Update height
        self.__update_height()
        
        # Get balance factor
        balance = self.__get_balance()
        
        # Left-Left case
        if balance > 1 and (self.left and self.__get_left_balance() >= 0):
            return self.__rotate_right()
        
        # Left-Right case
        if balance > 1 and (self.left and self.__get_left_balance() < 0):
            self.left = self.left.__rotate_left()
            return self.__rotate_right()
        
        # Right-Right case
        if balance < -1 and (self.right and self.__get_right_balance() <= 0):
            return self.__rotate_left()
        
        # Right-Left case
        if balance < -1 and (self.right and self.__get_right_balance() > 0):
            self.right = self.right.__rotate_right()
            return self.__rotate_left()
        
        return self

    def __get_left_balance(self):
        """Get balance factor of left child"""
        return self.left.__get_balance() if self.left else 0
    
    def __get_right_balance(self):
        """Get balance factor of right child"""
        return self.right.__get_balance() if self.right else 0

    def insert(self, abs_position: int, ref_seq: str, ref_metadata) -> "SequenceNode":
        """
        Insert a reference sequence at the absolute position in the tree.

        Args:
            abs_position (int): Absolute position in the tree to insert at
            ref_seq (str): Reference sequence to insert
            ref_metadata: Metadata for the reference

        Returns:
            SequenceNode: Root node after insertion
        """
        # Calculate positions in tree
        left_length = self.left.total_length if self.left else 0

        # Fast path: Position is before this node
        if abs_position <= left_length:
            if self.left:
                self.left = self.left.insert(abs_position, ref_seq, ref_metadata)
            else:
                # Insert as left child
                self.left = SequenceNode(ref_seq, True, ref_metadata)
            self.__update_total_length()
            self.__update_height()
            return self.__balance()

        # Position is within this node
        node_start = left_length
        node_end = node_start + self.length

        # Fast path: Position is within this node
        if node_start < abs_position < node_end:
            # Split this node
            rel_pos = abs_position - node_start
            left_content = self.content[:rel_pos]
            right_content = self.content[rel_pos:]
            
            # Handle TSD generation if requested in metadata
            tsd_length = ref_metadata.get("tsd_length", 0) if ref_metadata else 0
            if tsd_length > 0:
                # Extract source TSD sequence from the original sequence
                # We extract from the right side of the split point (beginning of right_content)
                source_tsd_seq = right_content[:min(tsd_length, len(right_content))]
                
                # Generate TSD sequences (potentially with mutations)
                tsd_5, tsd_3 = generate_TSD(source_tsd_seq, tsd_length)
                
                # Remove source TSD from right_content as it will be duplicated
                if len(source_tsd_seq) > 0:
                    right_content = right_content[len(source_tsd_seq):]
                
                # Add TSD sequences to the left and right content
                left_content = left_content + tsd_5
                right_content = tsd_3 + right_content
            
            # Create new left child with left content
            new_left = SequenceNode(left_content, self.is_reference, self.metadata)
            if self.left:
                new_left.left = self.left
                new_left.__update_total_length()
                new_left.__update_height()

            # Create new right child with right content and original right child
            new_right = SequenceNode(right_content, self.is_reference, self.metadata)
            if self.right:
                new_right.right = self.right
                new_right.__update_total_length()
                new_right.__update_height()

            # Replace this node's content with the reference
            self.content = ref_seq
            self.length = len(ref_seq)
            self.is_reference = True
            self.metadata = ref_metadata

            # Set new children
            self.left = new_left
            self.right = new_right
            self.__update_total_length()
            self.__update_height()

            return self.__balance()

        # Fast path: Position is after this node
        if abs_position >= node_end:
            if self.right:
                self.right = self.right.insert(abs_position - node_end, ref_seq, ref_metadata)
            else:
                # Insert as right child
                self.right = SequenceNode(ref_seq, True, ref_metadata)
            self.__update_total_length()
            self.__update_height()
            return self.__balance()

        raise RuntimeError("[ERROR] Should not reach here")

    def collect_content_in_order_traversal(self):
        """
        Helper method to collect node content in in-order traversal.
        """
        result = []
        if self.left:
            result.extend(self.left.collect_content_in_order_traversal())

        result.append(self.content)

        if self.right:
            result.extend(self.right.collect_content_in_order_traversal())

        return result

    def collect_refs(self, seq_id: str, abs_position: int = 0):
        """
        Collect reference nodes and calculate their absolute positions.

        Args:
            seq_id (str): ID of the original sequence
            abs_position (int): Current absolute position in the concatenated sequence

        Returns:
            Tuple[List[SeqRecord], int]: A tuple containing:
                - List of reference records
                - The total length of this subtree
        """
        ref_records = []

        left_length = 0
        if self.left:
            left_refs, left_length = self.left.collect_refs(seq_id, abs_position)
            ref_records.extend(left_refs)

        current_position = abs_position + left_length

        if self.is_reference:
            # Create a reference record with absolute position
            start_index = current_position
            end_index = start_index + self.length
            ref_id = f"{seq_id}_{start_index}_{end_index}-+-{self.length}"
            ref_record = create_sequence_record(self.content, ref_id)
            ref_records.append(ref_record)

        right_length = 0
        if self.right:
            right_refs, right_length = self.right.collect_refs(seq_id, current_position + self.length)
            ref_records.extend(right_refs)

        return ref_records, left_length + self.length + right_length


# ======================================================================================================================


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

    for path_list, weight in zip(path_list, weight_list or [None] * len(path_list)):
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

    # Uniform weights by default: If no weight provided, set all weights to 1
    if not weight_list:
        all_ref_weight_list = [1] * len(all_ref_sequence_list)

    return all_ref_sequence_list, all_ref_len_list, all_ref_weight_list


class SeqGenerator:
    """
    Class for generating sequences with random insertions from reference libraries.
    """

    def __init__(self, input_file: str, insertion: int, batch: int, processors: int, output_dir: str,
                 ref_lib: Optional[List[str]] = None, ref_lib_weight: Optional[List[float]] = None,
                 ref_len_limit: Optional[int] = None, flag_filter_n: bool = False, flag_track: bool = False,
                 tsd_length: Optional[int] = None, flag_visual: bool = False):
        """
        Initialize the sequence generator.

        Args:
            input_file (str): Path to the input sequence file
            insertion (int): Number of insertions per sequence
            batch (int): Number of independent result files to generate
            processors (int): Number of processors to use
            output_dir (str): Directory to save output files
            ref_lib (List[str], optional): List of reference library file paths
            ref_lib_weight (List[float], optional): Weights for reference libraries
            ref_len_limit (int, optional): Maximum length limit for reference sequences to load
            flag_filter_n (bool): Whether to filter out sequences containing N
            flag_track (bool): Whether to track reference sequences used
            tsd_length (int, optional): Length of Target Site Duplication (TSD) to generate
            flag_visual (bool): Whether to generate graphviz visualization of the sequence tree
        """
        self.input_file = input_file
        self.insertion = insertion
        self.batch = batch
        self.processors = processors
        self.output_dir = output_dir
        self.ref_len_limit = ref_len_limit
        self.flag_filter_n = flag_filter_n
        self.flag_track = flag_track
        self.tsd_length = tsd_length
        self.flag_visual = flag_visual

        # Load input sequence
        self.input = list(SeqIO.parse(input_file, "fasta"))
        if not self.input:
            raise ValueError(f"[ERROR] No sequences found in input file {input_file}")

        # Load reference sequences
        if ref_lib:
            ref_lib_data = _load_multiple_ref_libs(ref_lib, ref_lib_weight, ref_len_limit, flag_filter_n)
            if ref_lib_data:
                self.ref_sequences, _, self.ref_weights = ref_lib_data
            else:
                raise SystemExit("[ERROR] No valid reference sequences loaded")
        else:
            raise SystemExit("[ERROR] Reference library is required for insertion")

    def __pre_check(self):
        """
        Perform pre-execution checks.
        """
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

    def __process_single_sequence(self, seq_record: SeqRecord) -> Tuple[SeqRecord, Optional[List[SeqRecord]]]:
        """
        Process a single sequence with random insertions.
        
        Args:
            seq_record (SeqRecord): The input sequence record
            
        Returns:
            Tuple[SeqRecord, Optional[List[SeqRecord]]]: 
                The processed sequence and used references if tracking
        """
        # Start with a single node containing the original sequence
        root = SequenceNode(str(seq_record.seq))
        
        # Get the current total sequence length
        total_length = root.total_length
        
        # Generate all insertion positions and reference sequences at once
        insert_positions = random.choices(range(total_length + 1), k=self.insertion)
        selected_refs = random.choices(self.ref_sequences, weights=self.ref_weights, k=self.insertion)
        
        # Sort positions in descending order to maintain position integrity during insertion
        # This ensures earlier insertions don't affect the positions of later insertions
        insert_positions, selected_refs = sort_multiple_lists(insert_positions, selected_refs, reverse=True)
        
        # Insert references directly into the tree in optimal order
        for pos, ref_seq in zip(insert_positions, selected_refs):
            metadata = {"tsd_length": self.tsd_length} if self.tsd_length else {}
            root = root.insert(pos, ref_seq, metadata)
        
        # Generate final sequence
        final_seq = str(root)
        new_id = f"{seq_record.id}_ins{self.insertion}"
        if self.tsd_length:
            new_id += f"_tsd{self.tsd_length}"
        new_seq_record = create_sequence_record(final_seq, new_id)
        
        # Generate graphviz visualization if requested
        if self.flag_visual:
            graphviz_str = root.to_graphviz(node_id_prefix=seq_record.id)
            visual_dir_path = os.path.join(self.output_dir, TREE_VISUAL_DIR_NAME)
            visual_file_path = os.path.join(visual_dir_path, f"{seq_record.id}_tree_visual.dot")
            with open(visual_file_path, "w") as f:
                f.write(graphviz_str)
            print(f"Generated tree visualization for {seq_record.id}")
        
        # Only collect reference sequences if tracking is enabled
        used_refs = None
        if self.flag_track:
            used_refs, _ = root.collect_refs(seq_record.id)
        
        return new_seq_record, used_refs

    def _process_chunk(self, chunk_idx: int, chunk_size: int, total_sequences: int) -> Tuple[List[SeqRecord], 
                                                                                          Optional[List[SeqRecord]]]:
        """
        Process a chunk of sequences.
        
        Args:
            chunk_idx (int): The chunk index
            chunk_size (int): Number of sequences in each chunk
            total_sequences (int): Total number of sequences to process
            
        Returns:
            Tuple[List[SeqRecord], Optional[List[SeqRecord]]]: 
                The processed sequences and used references if tracking
        """
        start_idx = chunk_idx * chunk_size
        end_idx = min(start_idx + chunk_size, total_sequences)
        
        processed_seqs = []
        used_refs = [] if self.flag_track else None
        
        for i in range(start_idx, end_idx):
            seq_record = self.input[i]
            processed_seq, refs = self.__process_single_sequence(seq_record)
            processed_seqs.append(processed_seq)
            
            if self.flag_track and refs:
                used_refs.extend(refs)
            
        return processed_seqs, used_refs

    def __process_batch_multiprocessing(self, total_sequences: int) -> Tuple[List[SeqRecord], 
                                                                              Optional[List[SeqRecord]]]:
        """
        Process a batch of sequences using multiprocessing.
        
        Args:
            total_sequences (int): Total number of sequences to process
            
        Returns:
            Tuple[List[SeqRecord], Optional[List[SeqRecord]]]: 
                All processed sequences and references
        """
        # Calculate proper chunk size for work division
        chunk_size = max(1, total_sequences // (self.processors * 4))
        
        # Calculate number of chunks
        num_chunks = (total_sequences + chunk_size - 1) // chunk_size
        print(f"Using multiprocessing with {self.processors} worker(s) across {num_chunks} chunk(s)")
        
        # Prepare arguments for each chunk
        chunk_args = [(i, chunk_size, total_sequences) for i in range(num_chunks)]
        
        # Process chunks in parallel
        all_sequences = []
        all_refs = [] if self.flag_track else None
        
        with mp.Pool(self.processors) as pool:
            results_iter = pool.starmap(self._process_chunk, chunk_args)
            
            for i, (seqs, refs) in enumerate(results_iter, 1):
                all_sequences.extend(seqs)
                if self.flag_track and refs:
                    all_refs.extend(refs)
                print(f"Completed chunk {i}/{num_chunks}")
        
        return all_sequences, all_refs

    def __process_batch_single_thread(self, total_sequences: int) -> Tuple[List[SeqRecord], 
                                                                           Optional[List[SeqRecord]]]:
        """
        Process a batch of sequences using a single thread.
        
        Args:
            total_sequences (int): Total number of sequences to process
            
        Returns:
            Tuple[List[SeqRecord], Optional[List[SeqRecord]]]: 
                All processed sequences and references
        """
        all_sequences = []
        all_refs = [] if self.flag_track else None
        
        for i, seq_record in enumerate(self.input):
            processed_seq, refs = self.__process_single_sequence(seq_record)
            all_sequences.append(processed_seq)
            
            if self.flag_track and refs:
                all_refs.extend(refs)
            
            if (i + 1) % 10 == 0 or i + 1 == total_sequences:
                print(f"Processed {i+1}/{total_sequences} sequences")
        
        return all_sequences, all_refs

    def __process_single_batch(self, batch_num: int) -> float:
        """
        Process a single batch of sequences.
        
        Args:
            batch_num (int): The batch number (1-based index)
            
        Returns:
            float: The time taken to process the batch in seconds
        """
        batch_start_time = time.time()
        
        # Set a different random seed for each batch for result diversity
        random.seed(int(time.time()) + batch_num)
        
        print(f"\nProcessing batch {batch_num}/{self.batch}")
        
        # Calculate total number of input sequences to process
        total_sequences = len(self.input)
        
        # Divide work based on number of processors
        if self.processors > 1 and total_sequences > 1:
            all_sequences, all_refs = self.__process_batch_multiprocessing(total_sequences)
        else:
            all_sequences, all_refs = self.__process_batch_single_thread(total_sequences)
        
        os.makedirs(self.output_dir, exist_ok=True)
        print(f"Output directory: {self.output_dir}")

        # Save results for this batch
        self.__save_batch_results(self.output_dir, all_sequences, all_refs, batch_num)
        
        batch_elapsed_time = time.time() - batch_start_time
        print(f"Batch {batch_num} completed in {batch_elapsed_time:.2g} seconds")
        
        return batch_elapsed_time

    def __save_batch_results(self, output_dir: str, sequences: List[SeqRecord], 
                              references: Optional[List[SeqRecord]], 
                              batch_num: int = 1):
        """
        Save the batch results to output files.
        
        Args:
            output_dir (str): Directory to save the results
            sequences (List[SeqRecord]): List of sequence records to save
            references (Optional[List[SeqRecord]]): List of reference records to save if tracking is enabled
            batch_num (int): The batch number (1-based index)
        """
        # Add batch suffix to filenames if multiple batches
        suffix = f"_batch{batch_num}" if self.batch > 1 else ""
        print(f"Saving {len(sequences)} processed sequences")
        output_dict = {f"sequences{suffix}.fasta": sequences}
        
        if self.flag_track and references:
            print(f"Saving {len(references)} reference sequence records")
            output_dict[f"used_refs{suffix}.fasta"] = references
        
        save_multi_fasta_from_dict(output_dict, output_dir)

    def __print_header(self):
        """
        Print the program header with basic information.
        """
        print(f"=== RandSeqInsert {VERSION} ===")
        print(f"Processing input file: {self.input_file}")
        print(f"Number of input sequences: {len(self.input)}")
        print(f"Insertion settings: {self.insertion} insertions per sequence")
        if self.tsd_length:
            print(f"TSD settings: Generating TSD of length {self.tsd_length} at insertion sites")
        print(f"Reference library: {len(self.ref_sequences)} sequences loaded")
        print(f"Generating {self.batch} independent result file(s)")
        if self.flag_visual:
            print(f"Tree visualization enabled: DOT files will be generated for each sequence")

    def __print_summary(self, total_elapsed_time: float):
        """
        Print a summary of the processing results.
        
        Args:
            total_elapsed_time (float): Total time taken for all batches
        """
        print(f"\nAll batches completed in {total_elapsed_time:.2g} seconds")
        
        if self.batch > 1:
            print(f"Results saved to {os.path.abspath(self.output_dir)} with filenames sequences_batch[1-{self.batch}].fasta")
        else:
            print(f"Results saved to {os.path.abspath(self.output_dir)}")

    def execute(self):
        """
        Execute the sequence generation process.
        """
        self.__print_header()
        
        start_time = time.time()
        
        # Ensure output directory exists for the base path
        self.__pre_check()
        
        # Process each batch (separate run)
        for batch_num in range(1, self.batch + 1):
            self.__process_single_batch(batch_num)
        
        total_elapsed_time = time.time() - start_time
        self.__print_summary(total_elapsed_time)


def main():
    """Main function to handle command line arguments and execute the program."""
    parser = argparse.ArgumentParser(prog="RandSeqInsert",
                                     description="RandomSequenceInsertion: A high-performance Python tool for randomly inserting "
                                               "genomic fragments from reference libraries into target sequences. "
                                               "It supports multiprocessing for efficient processing of large datasets.",
                                     epilog="")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}\n{INFO}")

    # Core Arguments
    core_group = parser.add_argument_group("Core Arguments")
    core_group.add_argument("-i", "--input",
                       help="Input sequence file in FASTA format. Contains the target sequences to be inserted into.",
                       type=str, required=True, metavar="FILE")
    core_group.add_argument("-o", "--output", default=DEFAULT_OUTPUT_DIR_ABS_PATH, metavar="DIR",
                       help=f"Output directory path. Generated sequences and related files will be saved here. Default: '{DEFAULT_OUTPUT_DIR_ABS_PATH}'")

    core_group.add_argument("-is", "--insertion", metavar="INT",
                       help="Number of insertions per sequence. Specifies how many reference sequence fragments to insert into each input sequence.",
                       type=int, required=True)

    # Reference Library Arguments
    ref_group = parser.add_argument_group("Reference Library Arguments")
    ref_group.add_argument("-r", "--reference", nargs="+", metavar="FILE/DIR", required=True,
                            help="Reference sequence library file or directory paths. Multiple FASTA format reference files can be specified. Sequences from these files will be used as insertion fragments.")
    ref_group.add_argument("-w", "--weight", type=float, nargs="+", metavar="FLOAT",
                       help="Weights for reference libraries. Controls the probability of selecting sequences from different reference libraries. The number of weights should match the number of reference libraries.")
    ref_group.add_argument("-l", "--limit", type=int, default=None, metavar="INT",
                       help="Reference sequence length limit. Only loads reference sequences with length less than or equal to this value. Default: no limit.")

    # Control Arguments
    ctrl_group = parser.add_argument_group("Control Arguments")
    ctrl_group.add_argument("-b", "--batch", type=int, default=1, metavar="INT",
                            help="Number of independent result files to generate. Runs the entire process multiple times with different random seeds to generate multiple output sets. Default: 1")
    ctrl_group.add_argument("-p", "--processors", type=int, default=DEFAULT_ALLOCATED_CPU_CORES, metavar="INT",
                       help=f"Number of processors to use for parallel processing. Default: {DEFAULT_ALLOCATED_CPU_CORES}")

    # Flags
    flag_group = parser.add_argument_group("Flags")
    flag_group.add_argument("--filter_n", action="store_true",
                                help="Filter out sequences containing N. Enable this option to exclude reference sequences containing N bases.")
    flag_group.add_argument("--track", action="store_true",
                       help="Track and save used reference sequences. Enable this option to generate an additional FASTA file in the output directory recording all used reference sequences and their insertion positions.")
    flag_group.add_argument("--tsd", type=int, metavar="LENGTH",
                       help="Enable Target Site Duplication (TSD) with specified length. When a reference sequence is inserted, TSD of this length will be generated at the insertion site.")
    flag_group.add_argument("--visual", action="store_true",
                       help="Generate Graphviz DOT files visualizing the tree structure of each sequence. Files will be named {seqid}_tree_visual.dot and saved in the output directory.")

    parsed_args = parser.parse_args()

    input_file = parsed_args.input
    insertion = parsed_args.insertion
    batch = parsed_args.batch
    processors = parsed_args.processors
    output_dir_path = parsed_args.output
    ref_lib = parsed_args.reference
    ref_lib_weight = parsed_args.weight
    ref_len_limit = parsed_args.limit
    flag_filter_n = parsed_args.filter_n
    flag_track = parsed_args.track
    tsd_length = parsed_args.tsd
    flag_visual = parsed_args.visual

    generator = SeqGenerator(
        input_file=input_file,
        insertion=insertion,
        batch=batch,
        processors=processors,
        output_dir=output_dir_path,
        ref_lib=ref_lib,
        ref_lib_weight=ref_lib_weight,
        ref_len_limit=ref_len_limit,
        flag_filter_n=flag_filter_n,
        flag_track=flag_track,
        tsd_length=tsd_length,
        flag_visual=flag_visual
    )
    generator.execute()

if __name__ == "__main__":
    main()
