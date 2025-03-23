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
import random
from typing import List, Optional, Union, Callable, Tuple, Dict
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

VERSION = "v0.1.0"
INFO = ("by Tianyu (Sky) Lu (tlu83@wisc.edu)\n"
        "")

# PRE-DEFINED PARAMETERS
NUCLEOTIDES = ["A", "T", "G", "C"]
DEFAULT_ALLOCATED_CPU_CORES = os.cpu_count() - 2 if os.cpu_count() > 2 else 1
DEFAULT_OUTPUT_DIR_REL_PATH = "RandSeqInsert-Result"

# 可视化参数
VIZ_MAX_WIDTH = 80  # 可视化输出的最大宽度（字符数）
VIZ_SEQ_PREVIEW_LENGTH = 4  # 在可视化中显示序列内容的前后碱基数量
VIZ_PAGE_SIZE_FACTOR = 25  # 每页显示的序列长度 = VIZ_MAX_WIDTH * VIZ_PAGE_SIZE_FACTOR
VIZ_MAX_SYMBOLS = 40  # 序列可视化符号的最大数量

PROGRAM_ROOT_DIR_ABS_PATH = os.path.dirname(__file__)
DEFAULT_OUTPUT_DIR_ABS_PATH = os.path.join(os.getcwd(), DEFAULT_OUTPUT_DIR_REL_PATH)


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


class SequenceNode:
    """
    Node in a sequence tree structure, used for efficient sequence insertion operations.
    """
    def __init__(self, content: str, is_reference: bool = False, metadata = None):
        """
        Initialize a sequence node.
        
        Args:
            content (str): The sequence content
            is_reference (bool): Whether this node contains a reference sequence
            metadata: Additional information about the node (such as iteration and position)
        """
        self.content = content
        self.length = len(content)
        self.is_reference = is_reference
        self.metadata = metadata
        self.left = None
        self.right = None
        # Total length of the subtree for efficient traversal
        self.total_length = self.length
    
    def __update_total_length(self):
        """
        Update the total length of this subtree.
        
        This method recalculates the total length of the subtree rooted at this node
        by summing the lengths of the left subtree, the current node, and the right subtree.
        """
        left_length = self.left.total_length if self.left else 0
        right_length = self.right.total_length if self.right else 0
        self.total_length = left_length + self.length + right_length
    
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
            return self
        
        # Position is within this node
        node_start = left_length
        node_end = node_start + self.length
        
        # Fast path: Position is within this node
        if node_start < abs_position < node_end:
            # Split this node
            rel_pos = abs_position - node_start
            left_content = self.content[:rel_pos]
            right_content = self.content[rel_pos:]
            
            # Create new left child with left content
            new_left = SequenceNode(left_content, self.is_reference, self.metadata)
            if self.left:
                new_left.left = self.left
                new_left.__update_total_length()
            
            # Create new right child with right content and original right child
            new_right = SequenceNode(right_content, self.is_reference, self.metadata)
            if self.right:
                new_right.right = self.right
                new_right.__update_total_length()
            
            # Replace this node's content with the reference
            self.content = ref_seq
            self.length = len(ref_seq)
            self.is_reference = True
            self.metadata = ref_metadata
            
            # Set new children
            self.left = new_left
            self.right = new_right
            self.__update_total_length()
            
            return self
        
        # Fast path: Position is after this node
        if abs_position >= node_end:
            if self.right:
                self.right = self.right.insert(abs_position - node_end, ref_seq, ref_metadata)
            else:
                # Insert as right child
                self.right = SequenceNode(ref_seq, True, ref_metadata)
            self.__update_total_length()
            return self
        
        return self  # Should not reach here
    
    def __str__(self) -> str:
        """
        Convert the tree to a string by in-order traversal.
        
        Returns:
            str: The concatenated sequence
        """
        result = []
        self.__collect_content(result)
        return "".join(result)
    
    def __collect_content(self, result: list):
        """
        Helper method to collect node content in in-order traversal.
        
        Args:
            result (list): List to store content strings
        """
        if self.left:
            self.left.__collect_content(result)
        
        result.append(self.content)
        
        if self.right:
            self.right.__collect_content(result)
    
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
    
    def _create_node_info(self, current_position: int) -> tuple:
        """
        Create node information.
        
        Args:
            current_position (int): The absolute position of the current node in the sequence
            
        Returns:
            tuple: Contains node label, node width, and centered node line
        """
        # Node type and position information
        node_type = "REF" if self.is_reference else "SEQ"
        end_position = current_position + self.length - 1  # -1 for 0-based index
        
        # Metadata information
        metadata_str = ""
        if self.metadata and self.is_reference:
            metadata_str = f" [iter:{self.metadata.get('iteration', '?')},pos:{self.metadata.get('rel_pos', '?')}]"
        
        # Format node label
        node_label = f"{node_type}|{current_position}-{end_position}|{self.length}{metadata_str}"
        
        # Calculate node width and ensure frame can accommodate text
        node_width = len(node_label) + 4  # Add frame and some padding
        
        # Center node label
        node_line = node_label.center(node_width - 2, ' ')  # -2 for frame characters
        
        return node_label, node_width, node_line
    
    def _create_node_box(self, node_line: str, node_width: int) -> list:
        """
        Create a node box.
        
        Args:
            node_line (str): Node content line
            node_width (int): Node width
            
        Returns:
            list: List of strings containing the node box
        """
        return [
            "┌" + "─" * (node_width - 2) + "┐",
            "│" + node_line + "│",
            "└" + "─" * (node_width - 2) + "┘"
        ]
    
    def _create_leaf_visual(self, node_width: int, node_line: str) -> tuple:
        """
        Create visual representation for a leaf node.
        
        Args:
            node_width (int): Node width
            node_line (str): Node content line
            
        Returns:
            tuple: Contains list of strings representing the visual, node length, and height
        """
        result = self._create_node_box(node_line, node_width)
        return result, self.length, 1
    
    def _create_fork_line(self, parent_center: int, left_conn_pos: int, right_conn_pos: int, total_width: int) -> str:
        """
        Create a fork line connecting two child nodes.
        
        Args:
            parent_center (int): Center position of the parent node
            left_conn_pos (int): Connection position for the left child
            right_conn_pos (int): Connection position for the right child
            total_width (int): Total width
            
        Returns:
            str: Fork line string
        """
        # Determine fork center position
        if left_conn_pos <= parent_center <= right_conn_pos:
            # Parent between children - ideal case
            fork_center = parent_center
        elif parent_center < left_conn_pos:
            # Parent to the left of both children
            fork_center = left_conn_pos
        else:
            # Parent to the right of both children
            fork_center = right_conn_pos
        
        # Build fork line
        parts = []
        
        # Add left padding
        if left_conn_pos > 0:
            parts.append(" " * (left_conn_pos - 1))
            
        # Add left fork
        parts.append("┌")
        
        # Add horizontal line to fork center
        if fork_center > left_conn_pos + 1:
            parts.append("─" * (fork_center - left_conn_pos - 1))
            
        # Add fork connection point
        parts.append("┴")
        
        # Add horizontal line to right connection
        if right_conn_pos > fork_center + 1:
            parts.append("─" * (right_conn_pos - fork_center - 1))
            
        # Add right fork
        parts.append("┐")
        
        fork_line = "".join(parts)
        
        # Add padding to match total width
        return fork_line.ljust(total_width)
    
    def _align_child_trees(self, left_lines: list, right_lines: list, left_box_width: int, right_box_width: int, 
                          right_box_pos: int, total_width: int) -> list:
        """
        Align and merge left and right subtree visualizations.
        
        Args:
            left_lines (list): List of lines for left subtree
            right_lines (list): List of lines for right subtree
            left_box_width (int): Width of left subtree box
            right_box_width (int): Width of right subtree box
            right_box_pos (int): Position of right subtree box
            total_width (int): Total width
            
        Returns:
            list: Merged subtree lines
        """
        result = []
        max_lines = max(len(left_lines), len(right_lines))
        
        for i in range(max_lines):
            # Add left subtree line
            if i < len(left_lines):
                line = left_lines[i]
            else:
                line = " " * left_box_width
                
            # Add spacing
            line = line.ljust(right_box_pos)
            
            # Add right subtree line
            if i < len(right_lines):
                line += right_lines[i]
            else:
                line += " " * right_box_width
                
            # Ensure line matches total width
            result.append(line.ljust(total_width))
            
        return result
    
    def _create_single_child_connector(self, parent_center: int, child_center: int, is_left_child: bool, 
                                      total_width: int) -> str:
        """
        Create a connector line to a single child node.
        
        Args:
            parent_center (int): Center position of the parent node
            child_center (int): Center position of the child node
            is_left_child (bool): Whether it's a left child
            total_width (int): Total width
            
        Returns:
            str: Connector line string
        """
        if parent_center == child_center:
            # Centers aligned, use simple connector
            return (" " * parent_center + "┴").ljust(total_width)
        
        # Need horizontal connector
        if is_left_child:
            # Left child connection
            if parent_center > child_center:
                # Parent to the right of child
                connector = " " * (child_center - 1) + "┌" + "─" * (parent_center - child_center - 1) + "┴"
            else:
                # Parent to the left of child
                connector = " " * parent_center + "┴" + "─" * (child_center - parent_center - 1) + "┐"
        else:
            # Right child connection
            if parent_center < child_center:
                # Parent to the left of child
                connector = " " * (parent_center - 1) + "┌" + "─" * (child_center - parent_center - 1) + "┐"
            else:
                # Parent to the right of child
                connector = " " * (child_center - 1) + "┌" + "─" * (parent_center - child_center - 1) + "┴"
                
        return connector.ljust(total_width)
    
    def generate_tree_visual(self, abs_position: int = 0):
        """
        Generate a visual representation of the tree structure including node information.
        
        Args:
            abs_position (int): Current absolute position in the concatenated sequence
            
        Returns:
            Tuple[List[str], int, int]: Contains:
                - List of strings representing the tree structure
                - Total length of this subtree
                - Height of this subtree
        """
        # First recursively process left and right subtrees
        left_lines, left_length, left_height = [], 0, 0
        if self.left:
            left_lines, left_length, left_height = self.left.generate_tree_visual(abs_position)
        
        current_position = abs_position + left_length
        
        right_lines, right_length, right_height = [], 0, 0
        if self.right:
            right_lines, right_length, right_height = self.right.generate_tree_visual(current_position + self.length)
        
        # Create node information
        node_label, node_width, node_line = self._create_node_info(current_position)
        
        # Calculate this subtree height
        height = max(left_height, right_height) + 1
        
        # Handle leaf node case
        if not self.left and not self.right:
            return self._create_leaf_visual(node_width, node_line)
        
        # Handle nodes with children
        # Create node box
        result = self._create_node_box(node_line, node_width)
        
        # Calculate subtree width
        left_box_width = len(left_lines[0]) if left_lines else 0
        right_box_width = len(right_lines[0]) if right_lines else 0
        parent_center = node_width // 2
        
        if self.left and self.right:
            # Both children exist
            # Calculate child node center positions
            left_child_center = left_box_width // 2
            right_child_center = right_box_width // 2
            
            # Calculate layout
            spacing = 8  # Left and right subtree spacing
            total_width = left_box_width + spacing + right_box_width
            right_box_pos = left_box_width + spacing
            
            # Calculate connection line positions
            left_conn_pos = left_child_center
            right_conn_pos = right_box_pos + right_child_center
            
            # Create and add fork line
            fork_line = self._create_fork_line(parent_center, left_conn_pos, right_conn_pos, total_width)
            result.append(fork_line)
            
            # Add and align child trees
            child_lines = self._align_child_trees(
                left_lines, right_lines, left_box_width, right_box_width, right_box_pos, total_width)
            result.extend(child_lines)
            
        elif self.left:
            # Only left child
            left_child_center = left_box_width // 2
            total_width = max(node_width, left_box_width)
            
            # Create connection line
            connector = self._create_single_child_connector(parent_center, left_child_center, True, total_width)
            result.append(connector)
            
            # Add left subtree, ensure width consistency
            result.extend([line.ljust(total_width) for line in left_lines])
            
        elif self.right:
            # Only right child
            right_child_center = right_box_width // 2
            right_offset = max(0, parent_center - right_child_center)
            total_width = max(node_width, right_offset + right_box_width)
            
            # Calculate right subtree connection center
            right_tree_center = right_offset + right_child_center
            
            # Create connection line
            connector = self._create_single_child_connector(parent_center, right_tree_center, False, total_width)
            result.append(connector)
            
            # Add right subtree, ensure correct alignment and width
            padded_right_lines = [" " * right_offset + line for line in right_lines]
            result.extend([line.ljust(total_width) for line in padded_right_lines])
        
        return result, left_length + self.length + right_length, height

    def collect_sequence_structure(self, abs_position=0):
        """
        收集序列结构信息，用于可视化
        
        Args:
            abs_position (int): 节点在完整序列中的绝对起始位置
            
        Returns:
            list: 包含(节点，起始位置，结束位置，深度)元组的列表
        """
        result = []
        
        # 处理左子树
        if self.left:
            result.extend(self.left.collect_sequence_structure(abs_position))
        
        # 计算当前节点位置
        left_length = self.left.total_length if self.left else 0
        start_pos = abs_position + left_length
        end_pos = start_pos + self.length - 1
        
        # 添加当前节点
        result.append((self, start_pos, end_pos, 0))
        
        # 处理右子树
        if self.right:
            right_pos = start_pos + self.length
            result.extend(self.right.collect_sequence_structure(right_pos))
        
        return result
    
    def create_sequence_visualization(self, output_file, max_width=80, seq_id=None):
        """
        创建整个序列的可视化展示
        
        Args:
            output_file: 输出文件对象
            max_width (int): 输出的最大宽度
            seq_id (str): 序列ID
        """
        # 收集序列结构信息
        structure = self.collect_sequence_structure()
        
        # 按照位置排序
        structure.sort(key=lambda x: (x[1], -x[2]))
        
        # 计算嵌套深度
        node_depths = {}
        for i, (node, start, end, _) in enumerate(structure):
            # 查找该节点的所有父节点
            depth = 0
            for parent_node, parent_start, parent_end, _ in structure:
                if parent_node != node and parent_start <= start and parent_end >= end:
                    depth += 1
            # 更新节点深度
            structure[i] = (node, start, end, depth)
        
        # 如果提供了序列ID，显示标题
        if seq_id:
            print(f"输入序列 {seq_id} 结构可视化:", file=output_file)
            print("=" * max_width, file=output_file)
            print("", file=output_file)
        
        # 输出可视化
        total_length = structure[-1][2] + 1  # 最后一个节点的结束位置 + 1
        
        # 确定每页的长度范围 - 使用全局常量来动态计算
        page_size = max_width * VIZ_PAGE_SIZE_FACTOR  # 动态计算每页显示的序列长度
        page_count = (total_length + page_size - 1) // page_size
        
        for page in range(page_count):
            page_start = page * page_size
            page_end = min((page + 1) * page_size - 1, total_length - 1)
            
            # 计算当前页的长度
            current_page_length = page_end - page_start + 1
            
            # 打印页标题
            if page_count > 1:
                print(f"序列部分 {page+1}/{page_count}: 位置 {page_start}-{page_end}", file=output_file)
                print("-" * max_width, file=output_file)
            
            # 筛选当前页范围内的节点
            page_nodes = [(node, start, end, depth) for node, start, end, depth in structure 
                         if not (end < page_start or start > page_end)]
            
            # 对于每个节点，生成其可视化行
            for node, start, end, depth in page_nodes:
                # 计算在页面中的起始和结束位置
                vis_start = max(start, page_start)
                vis_end = min(end, page_end)
                
                # 计算节点在当前页中的可见长度
                visible_length = vis_end - vis_start + 1
                
                # 计算节点可见长度占当前页长度的百分比
                length_percentage = visible_length / current_page_length
                
                # 缩进显示嵌套关系
                indent = "  " * depth
                
                # 确定节点类型
                node_type = "REF" if node.is_reference else "SEQ"
                
                # 获取序列预览 - 修复短序列问题
                if len(node.content) <= VIZ_SEQ_PREVIEW_LENGTH * 2:
                    # 如果序列很短，仅显示一半内容在开始和结束
                    preview_length = len(node.content) // 2
                    if preview_length == 0:  # 处理极短序列
                        preview_length = 1 if len(node.content) > 0 else 0
                    start_preview = node.content[:preview_length]
                    end_preview = node.content[-preview_length:] if preview_length > 0 else ""
                else:
                    # 正常情况
                    start_preview = node.content[:VIZ_SEQ_PREVIEW_LENGTH]
                    end_preview = node.content[-VIZ_SEQ_PREVIEW_LENGTH:]
                
                # 计算可用宽度，确保不超出最大宽度
                available_width = max_width - len(indent) - len(str(start)) - len(str(end)) - \
                                 len(start_preview) - len(end_preview) - len(node_type) - len(str(node.length)) - 20
                
                # 确保available_width不为负数
                available_width = max(1, available_width)
                
                # 基于百分比和可用宽度计算符号数量，使用全局常量限制
                symbols_count = max(1, min(available_width, int(length_percentage * VIZ_MAX_SYMBOLS)))
                
                # 使用不同的符号表示不同类型的节点
                symbols = "+" * symbols_count if node_type == "SEQ" else "-" * symbols_count
                
                # 生成节点可视化行，确保其不超过最大宽度
                line = f"{indent}| {start} {start_preview} {symbols} {node_type} ({node.length}) {symbols} {end_preview} {end} |"
                
                # 如果行太长，截断符号
                if len(line) > max_width:
                    excess = len(line) - max_width + 3  # +3 为了添加"..."
                    symbols_new_length = max(1, symbols_count - excess // 2)
                    symbols = "=" * symbols_new_length if node_type == "SEQ" else "-" * symbols_new_length
                    line = f"{indent}| {start} {start_preview} {symbols} {node_type} ({node.length}) {symbols} {end_preview} {end} |"
                    if len(line) > max_width:
                        line = line[:max_width-3] + "..."
                
                print(line, file=output_file)
            
            print("", file=output_file)
            print("-" * max_width, file=output_file)
            print("", file=output_file)


class SeqGenerator:
    """
    Class for generating sequences with random insertions from reference libraries.
    """

    def __init__(self, input_file: str, insertion: int, iteration: int, batch: int, processors: int, output_dir: str,
                 ref_lib: Optional[List[str]] = None, ref_lib_weight: Optional[List[float]] = None,
                 ref_len_limit: Optional[int] = None, flag_filter_n: bool = False, flag_track: bool = False,
                 flag_visual: bool = False):
        """
        Initialize the sequence generator.

        Args:
            input_file (str): Path to the input sequence file
            insertion (int): Number of insertions per sequence
            iteration (int): Number of times to iterate the insertion process
            batch (int): Number of independent result files to generate
            processors (int): Number of processors to use
            output_dir (str): Directory to save output files
            ref_lib (List[str], optional): List of reference library file paths
            ref_lib_weight (List[float], optional): Weights for reference libraries
            ref_len_limit (int, optional): Maximum length limit for reference sequences to load
            flag_filter_n (bool): Whether to filter out sequences containing N
            flag_track (bool): Whether to track reference sequences used
            flag_visual (bool): Whether to generate tree visualization files
        """
        self.input_file = input_file
        self.insertion = insertion
        self.iteration = iteration
        self.batch = batch
        self.processors = processors
        self.output_dir = output_dir
        self.ref_len_limit = ref_len_limit
        self.flag_filter_n = flag_filter_n
        self.flag_track = flag_track
        self.flag_visualize = flag_visual

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

    def __process_single_sequence(self, seq_record: SeqRecord, seq_num: int) -> Tuple[SeqRecord, 
                                                                                      Optional[List[SeqRecord]],
                                                                                      Optional[SequenceNode]]:
        """
        Process a single sequence through all iterations using a tree-based algorithm.
        Directly navigates the tree for insertions without rebuilding it.
        
        This optimized implementation maintains a balanced tree structure for efficient
        insertion and traversal operations, especially for large sequences with many 
        insertion points.
        
        Args:
            seq_record (SeqRecord): The input sequence record
            seq_num (int): The sequence number
            
        Returns:
            Tuple[SeqRecord, Optional[List[SeqRecord]], Optional[SequenceNode]]: 
                The processed sequence, used references if tracking, and the root node if tracking
        """
        # Start with a single node containing the original sequence
        root = SequenceNode(str(seq_record.seq))
        
        for i in range(self.iteration):
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
                metadata = {'iteration': i + 1, 'rel_pos': pos}
                root = root.insert(pos, ref_seq, metadata)
        
        # 最终生成序列
        final_seq = str(root)
        new_id = f"{seq_record.id}_iter{self.iteration}_ins{self.insertion}"
        new_seq_record = create_sequence_record(final_seq, new_id)
        
        # 生成可视化（如果启用）
        if self.flag_visualize:
            self.visualize_sequence(seq_record.id, root, self.output_dir)
        
        # Only collect reference sequences if tracking is enabled
        used_refs = None
        if self.flag_track:
            used_refs, _ = root.collect_refs(seq_record.id)
            return new_seq_record, used_refs, root
        
        return new_seq_record, used_refs, None

    def _process_chunk(self, chunk_idx: int, chunk_size: int, total_sequences: int) -> Tuple[List[SeqRecord], 
                                                                                             Optional[List[SeqRecord]],
                                                                                             Optional[Dict[str, List[str]]]]:
        """
        Process a chunk of sequences.
        
        Args:
            chunk_idx (int): The chunk index
            chunk_size (int): Number of sequences in each chunk
            total_sequences (int): Total number of sequences to process
            
        Returns:
            Tuple[List[SeqRecord], Optional[List[SeqRecord]], Optional[Dict[str, List[str]]]]: 
                The processed sequences, used references if tracking, and tree visualizations if tracking
        """
        start_idx = chunk_idx * chunk_size
        end_idx = min(start_idx + chunk_size, total_sequences)
        
        processed_seqs = []
        used_refs = [] if self.flag_track else None
        tree_visuals = {} if self.flag_track else None
        
        for i in range(start_idx, end_idx):
            seq_record = self.input[i]
            processed_seq, refs, root_node = self.__process_single_sequence(seq_record, i + 1)
            processed_seqs.append(processed_seq)
            
            if self.flag_track:
                if refs:
                    used_refs.extend(refs)
                if root_node:
                    tree_lines, _, _ = root_node.generate_tree_visual()
                    tree_visuals[seq_record.id] = tree_lines
        
        return processed_seqs, used_refs, tree_visuals

    def __process_batch_multiprocessing(self, total_sequences: int) -> Tuple[List[SeqRecord], 
                                                                             Optional[List[SeqRecord]],
                                                                             Optional[Dict[str, List[str]]]]:
        """
        Process a batch of sequences using multiprocessing.
        
        Args:
            total_sequences (int): Total number of sequences to process
            
        Returns:
            Tuple[List[SeqRecord], Optional[List[SeqRecord]], Optional[Dict[str, List[str]]]]: 
                All processed sequences, references, and tree visualizations
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
        all_tree_visuals = {} if self.flag_track else None
        
        with mp.Pool(self.processors) as pool:
            results_iter = pool.starmap(self._process_chunk, chunk_args)
            
            for i, (seqs, refs, tree_visuals) in enumerate(results_iter, 1):
                all_sequences.extend(seqs)
                if self.flag_track:
                    if refs:
                        all_refs.extend(refs)
                    if tree_visuals:
                        all_tree_visuals.update(tree_visuals)
                print(f"Completed chunk {i}/{num_chunks}")
        
        return all_sequences, all_refs, all_tree_visuals

    def __process_batch_single_thread(self, total_sequences: int) -> Tuple[List[SeqRecord], 
                                                                          Optional[List[SeqRecord]],
                                                                          Optional[Dict[str, List[str]]]]:
        """
        Process a batch of sequences using a single thread.
        
        Args:
            total_sequences (int): Total number of sequences to process
            
        Returns:
            Tuple[List[SeqRecord], Optional[List[SeqRecord]], Optional[Dict[str, List[str]]]]: 
                All processed sequences, references, and tree visualizations
        """
        all_sequences = []
        all_refs = [] if self.flag_track else None
        all_tree_visuals = {} if self.flag_track else None
        
        for i, seq_record in enumerate(self.input):
            processed_seq, refs, root_node = self.__process_single_sequence(seq_record, i + 1)
            all_sequences.append(processed_seq)
            
            if self.flag_track:
                if refs:
                    all_refs.extend(refs)
                if root_node:
                    tree_lines, _, _ = root_node.generate_tree_visual()
                    all_tree_visuals[seq_record.id] = tree_lines
            
            if (i + 1) % 10 == 0 or i + 1 == total_sequences:
                print(f"Processed {i+1}/{total_sequences} sequences")
        
        return all_sequences, all_refs, all_tree_visuals

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
            all_sequences, all_refs, all_tree_visuals = self.__process_batch_multiprocessing(total_sequences)
        else:
            all_sequences, all_refs, all_tree_visuals = self.__process_batch_single_thread(total_sequences)
        
        os.makedirs(self.output_dir, exist_ok=True)
        print(f"Output directory: {self.output_dir}")

        # Save results for this batch
        self.__save_batch_results(self.output_dir, all_sequences, all_refs, all_tree_visuals, batch_num)
        
        batch_elapsed_time = time.time() - batch_start_time
        print(f"Batch {batch_num} completed in {batch_elapsed_time:.2g} seconds")
        
        return batch_elapsed_time

    def __save_batch_results(self, output_dir: str, sequences: List[SeqRecord], 
                             references: Optional[List[SeqRecord]], 
                             tree_visuals: Optional[Dict[str, List[str]]], 
                             batch_num: int = 1):
        """
        Save the batch results to output files.
        
        Args:
            output_dir (str): Directory to save the results
            sequences (List[SeqRecord]): List of sequence records to save
            references (Optional[List[SeqRecord]]): List of reference records to save if tracking is enabled
            tree_visuals (Optional[Dict[str, List[str]]]): Dictionary mapping sequence IDs to tree visualization lines
            batch_num (int): The batch number (1-based index)
        """
        # Add batch suffix to filenames if multiple batches
        suffix = f"_batch{batch_num}" if self.batch > 1 else ""
        print(f"Saving {len(sequences)} processed sequences")
        output_dict = {f"sequences{suffix}.fasta": sequences}
        
        if self.flag_track:
            if references:
                print(f"Saving {len(references)} reference sequence records")
                output_dict[f"used_refs{suffix}.fasta"] = references
            
            if tree_visuals:
                print(f"Saving tree structure visualizations")
                tree_file_path = os.path.join(output_dir, f"ins_tree_visual{suffix}.txt")
                with open(tree_file_path, 'w', encoding='utf-8') as tree_file:
                    for seq_id, lines in tree_visuals.items():
                        tree_file.write(f"=== Tree Structure Visualization: {seq_id} ===\n")
                        tree_file.write("Node Format: Type|Start-End|Length[Metadata]\n")
                        tree_file.write("SEQ: Original sequence node, REF: Inserted reference sequence node\n\n")
                        tree_file.write("\n".join(lines))
                        tree_file.write("\n\n" + "="*80 + "\n\n")
        
        save_multi_fasta_from_dict(output_dict, output_dir)

    def __print_header(self):
        """
        Print the program header with basic information.
        """
        print(f"=== RandSeqInsert {VERSION} ===")
        print(f"Processing input file: {self.input_file}")
        print(f"Number of input sequences: {len(self.input)}")
        print(f"Insertion settings: {self.insertion} insertions per sequence × {self.iteration} iterations")
        print(f"Reference library: {len(self.ref_sequences)} sequences loaded")
        print(f"Generating {self.batch} independent result file(s)")
        if self.flag_visualize:
            print(f"Visualization: Enabled (output to {os.path.join(self.output_dir, 'seq_visualization')})")

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
            
        if self.flag_visualize:
            print(f"Sequence visualizations saved to {os.path.join(os.path.abspath(self.output_dir), 'seq_visualization')}")

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

    def visualize_sequence(self, seq_id, root_node, output_dir):
        """
        为指定序列创建可视化文件
        
        Args:
            seq_id (str): 序列ID
            root_node (SequenceNode): 序列的根节点
            output_dir (str): 输出目录
        """
        # 创建可视化文件目录
        viz_dir = os.path.join(output_dir, "seq_visualization")
        os.makedirs(viz_dir, exist_ok=True)
        
        # 创建可视化文件
        viz_filename = f"{seq_id}_visualization.txt"
        viz_path = os.path.join(viz_dir, viz_filename)
        
        with open(viz_path, 'w') as f:
            # 添加标题
            f.write(f"序列可视化: {seq_id}\n")
            f.write("="*VIZ_MAX_WIDTH + "\n\n")
            f.write("说明:\n")
            f.write("- SEQ: 原始序列节点\n")
            f.write("- REF: 插入的参考序列节点\n")
            f.write("- 数字表示序列在完整序列中的位置\n")
            f.write(f"- 左右两侧显示序列的实际内容（前后各{VIZ_SEQ_PREVIEW_LENGTH}个碱基）\n")
            f.write("- 缩进表示嵌套层级\n\n")
            f.write("="*VIZ_MAX_WIDTH + "\n\n")
            
            # 使用新的可视化方法
            root_node.create_sequence_visualization(f, max_width=VIZ_MAX_WIDTH, seq_id=seq_id)


def main():
    """Main function to handle command line arguments and execute the program."""
    parser = argparse.ArgumentParser(prog="RandSeqInsert",
                                     description="RandomSequenceInsertion: A high-performance Python tool for randomly inserting "
                                               "genomic fragments from reference libraries into target sequences. "
                                               "It supports multiprocessing for efficient processing of large datasets.",
                                     epilog="")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}\n{INFO}")

    # TODO Revise the help information

    # Core Arguments
    core_group = parser.add_argument_group("Core Arguments")
    core_group.add_argument("-i", "--input",
                       help="Input sequence file in FASTA format. Contains the target sequences to be inserted into.",
                       type=str, required=True, metavar="FILE")
    core_group.add_argument("-o", "--output", default=DEFAULT_OUTPUT_DIR_REL_PATH, metavar="DIR",
                       help=f"Output directory path. Generated sequences and related files will be saved here. Default: '{DEFAULT_OUTPUT_DIR_REL_PATH}'")

    core_group.add_argument("-is", "--insertion", metavar="INT",
                       help="Number of insertions per sequence. Specifies how many reference sequence fragments to insert into each input sequence per iteration.",
                       type=int, required=True)
    core_group.add_argument("-it", "--iteration", metavar="INT", type=int, default=1,
                       help="Number of times to iterate the insertion process. Each iteration builds upon the results of the previous one.")

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
    flag_group.add_argument("--visual", action="store_true",
                       help="生成序列插入可视化文件，以文本形式展示原始序列和参考序列的嵌套关系")

    parsed_args = parser.parse_args()

    input_file = parsed_args.input
    insertion = parsed_args.insertion
    iteration = parsed_args.iteration
    batch = parsed_args.batch
    processors = parsed_args.processors
    output_dir_path = parsed_args.output
    ref_lib = parsed_args.reference
    ref_lib_weight = parsed_args.weight
    ref_len_limit = parsed_args.limit
    flag_filter_n = parsed_args.filter_n
    flag_track = parsed_args.track
    flag_visual = parsed_args.visual

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
        flag_track=flag_track,
        flag_visual=flag_visual
    )
    generator.execute()

if __name__ == "__main__":
    main()
