#!/usr/bin/env python3
"""
Random Sequence Insertion Generator

This program takes acceptor sequences and performs random insertions using sequences
from a donor library. It supports multiprocessing for efficient sequence generation.

Author: Tianyu Lu (tianyu@lu.fm)
Date: 2024-11-27
"""

import os
import argparse
import random
import re
from typing import Tuple, List, Dict, Iterable, Sequence, Union, Optional, Callable
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

VERSION = "v1.0.0"
INFO = ("by Tianyu (Sky) Lu (tianyu@lu.fm) "
        "released under GPLv3")

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
    new_base = random.choice("ATGC")

    if is_insertion:
        # Choose insertion position (n + 1 positions)
        pos = random.randrange(len(seq) + 1)
        return seq[:pos] + new_base + seq[pos:]
    else:
        # Choose deletion position (n positions)
        pos = random.randrange(len(seq))
        return seq[:pos] + seq[pos+1:]

def generate_TSD(seq_slice: str,
                 length: int = float("inf"),
                 snp: float = DEFAULT_TSD_SNP_MUTATION_RATE,
                 indel: float = DEFAULT_TSD_INDEL_MUTATION_RATE) -> Tuple[str, str]:
    """Generate TSD sequences with potential mutations.

    Args:
        seq_slice (str): Acceptor sequence fragment
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
    def __init__(self, data: str, is_donor: bool = False, attrs: dict = None, uid: int = None):
        """
        Initialize a sequence node.

        Args:
            data (str): The sequence string
            is_donor (bool): Whether this node contains a donor sequence
            attrs: Additional information about the node (such as position)
            uid (int): Unique identifier for this node
        """
        self.data = data
        self.length = len(data)
        self.is_donor = is_donor
        self.attrs = {} if attrs is None else attrs.copy()
        self.left = None
        self.right = None
        self.uid = uid

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
        return "".join(self.collect_data_in_order_traversal())

    def update_total_length(self):
        """
        Update the total length of this subtree.

        This method recalculates the total length of the subtree rooted at this node
        by summing the lengths of the left subtree, the current node, and the right subtree.
        """
        left_length = self.left.total_length if self.left else 0
        right_length = self.right.total_length if self.right else 0
        self.total_length = left_length + self.length + right_length

    def update_height(self):
        """
        Update the height of this node.
        """
        left_height = self.left.height if self.left else 0
        right_height = self.right.height if self.right else 0
        self.height = max(left_height, right_height) + 1

    def get_balance(self):
        """
        Calculate the balance factor of this node.

        Returns:
            int: Balance factor (left height - right height)
        """
        left_height = self.left.height if self.left else 0
        right_height = self.right.height if self.right else 0
        return left_height - right_height

    def rotate_right(self):
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
        self.update_height()
        self.update_total_length()
        new_root.update_height()
        new_root.update_total_length()

        return new_root

    def rotate_left(self):
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
        self.update_height()
        self.update_total_length()
        new_root.update_height()
        new_root.update_total_length()

        return new_root

    def balance(self):
        """
        Balance this node if needed.

        Returns:
            SequenceNode: The new root after balancing
        """
        # Update height
        self.update_height()

        # Get balance factor
        balance = self.get_balance()

        # Left-Left case
        if balance > 1 and (self.left and self.__get_left_balance() >= 0):
            return self.rotate_right()

        # Left-Right case
        if balance > 1 and (self.left and self.__get_left_balance() < 0):
            self.left = self.left.rotate_left()
            return self.rotate_right()

        # Right-Right case
        if balance < -1 and (self.right and self.__get_right_balance() <= 0):
            return self.rotate_left()

        # Right-Left case
        if balance < -1 and (self.right and self.__get_right_balance() > 0):
            self.right = self.right.rotate_right()
            return self.rotate_left()

        return self

    def __get_left_balance(self):
        """Get balance factor of left child"""
        return self.left.get_balance() if self.left else 0

    def __get_right_balance(self):
        """Get balance factor of right child"""
        return self.right.get_balance() if self.right else 0

    def collect_data_in_order_traversal(self):
        """
        Helper method to collect node data in in-order traversal.
        """
        result = []
        if self.left:
            result.extend(self.left.collect_data_in_order_traversal())

        result.append(self.data)

        if self.right:
            result.extend(self.right.collect_data_in_order_traversal())

        return result


class SequenceTree:
    """
    Manages a tree of SequenceNode objects, with its own UID management system.
    Provides high-level operations for sequence insertion and traversal.
    """
    def __init__(self, initial_seq: str, base_uid: int = 0):
        """
        Initialize a new sequence tree with a root node containing the initial sequence.

        Args:
            initial_seq (str): The initial sequence to store in the root node
            base_uid (int): Base UID for this tree's UID management system
        """
        # Initialize UID management system
        self.next_uid = base_uid
        self.available_uids = []
        self.node_dict = {}

        # Create the root node
        self.root = self._create_node(initial_seq, False)

    def __str__(self) -> str:
        """
        Convert the tree to a string.

        Returns:
            str: The concatenated sequence
        """
        return str(self.root) if self.root else ""

    def _get_next_uid(self) -> int:
        """Get the next available UID from the UID management system"""
        if self.available_uids:
            return self.available_uids.pop(0)

        uid = self.next_uid
        self.next_uid += 1
        return uid

    def _release_uid(self, uid: int):
        """Release an UID back to the UID management system"""
        if uid in self.node_dict:
            del self.node_dict[uid]
        self.available_uids.append(uid)

    def _create_node(self, data: str, is_donor: bool = False, attrs: dict = None) -> SequenceNode:
        """
        Create a new SequenceNode with a unique UID.

        Args:
            data (str): Sequence data
            is_donor (bool): Whether this node contains a donor sequence
            attrs (dict): Additional attributes

        Returns:
            SequenceNode: The newly created node
        """
        uid = self._get_next_uid()
        node = SequenceNode(data, is_donor, attrs, uid)
        self.node_dict[uid] = node
        return node

    def insert(self, abs_position: int, donor_seq: str, donor_attrs: dict) -> None:
        """
        Insert a donor sequence at the specified position.

        Args:
            abs_position (int): Absolute position for insertion
            donor_seq (str): Donor sequence to insert
            donor_attrs (dict): Attributes for the donor
        """
        self.root = self._insert_iterative(self.root, abs_position, donor_seq, donor_attrs)

    def _insert_iterative(self, node: SequenceNode, abs_position: int, donor_seq: str, donor_attrs: dict) -> SequenceNode:
        """
        Iteratively insert a donor sequence at the absolute position in the tree.

        Args:
            node (SequenceNode): Current root node
            abs_position (int): Absolute position in the tree to insert at
            donor_seq (str): Donor sequence to insert
            donor_attrs (dict): Attributes for the donor

        Returns:
            SequenceNode: New root node after insertion
        """
        current = node
        parent_stack = []
        path_directions = []  # Record the path direction from root to current node ('left' or 'right')

        # Iteratively find insertion position
        while True:
            left_length = current.left.total_length if current.left else 0
            node_start = left_length
            node_end = node_start + current.length

            # Case 1: Position is in the current node's left subtree
            if abs_position <= left_length:
                if current.left:
                    # Record parent node and direction for backtracking
                    parent_stack.append(current)
                    path_directions.append('left')
                    current = current.left
                else:
                    # Create new left child node
                    current.left = self._create_node(donor_seq, True, donor_attrs)
                    # Use _update_node helper method instead of directly calling private methods
                    self._update_node(current)
                    break

            # Case 2: Position is inside the current node
            elif node_start < abs_position < node_end:
                # Calculate relative position
                rel_pos = abs_position - node_start
                left_data = current.data[:rel_pos]
                right_data = current.data[rel_pos:]

                # Handle TSD generation
                tsd_length = donor_attrs.get("tsd_length", 0)
                if tsd_length > 0:
                    # Extract source TSD sequence from the original sequence
                    source_tsd_seq = right_data[:min(tsd_length, len(right_data))]

                    # Generate TSD sequences (potentially with mutations)
                    tsd_5, tsd_3 = generate_TSD(source_tsd_seq, tsd_length)

                    # Remove source TSD from right_data as it will be duplicated
                    if len(source_tsd_seq) > 0:
                        right_data = right_data[len(source_tsd_seq):]

                    # Add TSD sequences to the left and right data
                    left_data = left_data + tsd_5
                    right_data = tsd_3 + right_data

                # Get current node attributes
                current_attrs = current.attrs.copy()
                current_is_donor = current.is_donor
                current_uid = current.uid

                # Create new left child with left data
                new_left = self._create_node(left_data, current_is_donor, current_attrs.copy())
                if current.left:
                    new_left.left = current.left
                    # Use _update_node helper method
                    self._update_node(new_left)

                # Create new right child with right data and original right child
                new_right = self._create_node(right_data, current_is_donor, current_attrs.copy())
                if current.right:
                    new_right.right = current.right
                    # Use _update_node helper method
                    self._update_node(new_right)

                # Record nested insertion relationships
                if current_is_donor:
                    # Record that the donor was cut
                    # Define the half attribute (left or right)
                    new_left.attrs["half"] = ["L"] if "half" not in new_left.attrs else new_left.attrs["half"] + ["L"]
                    new_right.attrs["half"] = ["R"] if "half" not in new_right.attrs else new_right.attrs["half"] + ["R"]

                    # Prepare donor attrs to record that it's nested inside this donor
                    donor_attrs_copy = donor_attrs.copy()
                    if "nested_in" not in donor_attrs_copy:
                        donor_attrs_copy["nested_in"] = []
                    donor_attrs_copy["nested_in"].append(current_uid)

                    # Record which donor cut this node
                    if "cut_by" not in new_left.attrs:
                        new_left.attrs["cut_by"] = []
                    if "cut_by" not in new_right.attrs:
                        new_right.attrs["cut_by"] = []

                    # Create a new donor node and get its UID
                    new_donor_uid = self._get_next_uid()

                    new_left.attrs["cut_by"].append(new_donor_uid)
                    new_right.attrs["cut_by"].append(new_donor_uid)

                    # Release current node's UID and set new values
                    self._release_uid(current.uid)
                    current.data = donor_seq
                    current.length = len(donor_seq)
                    current.is_donor = True
                    current.attrs = donor_attrs_copy
                    current.uid = new_donor_uid
                    self.node_dict[new_donor_uid] = current
                else:
                    # Replace current node's data with the donor sequence
                    current.data = donor_seq
                    current.length = len(donor_seq)
                    current.is_donor = True
                    current.attrs = donor_attrs.copy()

                # Set new children
                current.left = new_left
                current.right = new_right
                # Use _update_node helper method
                self._update_node(current)
                break

            # Case 3: Position is in the current node's right subtree
            elif abs_position >= node_end:
                if current.right:
                    # Record parent node and direction for backtracking
                    parent_stack.append(current)
                    path_directions.append('right')
                    # Adjust absolute position to fit the right subtree's relative position
                    abs_position -= node_end
                    current = current.right
                else:
                    # Create new right child node
                    current.right = self._create_node(donor_seq, True, donor_attrs)
                    # Use _update_node helper method
                    self._update_node(current)
                    break
            else:
                raise RuntimeError("[ERROR] Should not reach here")

        # Backtrack and update node heights and total lengths while executing balance
        while parent_stack:
            parent = parent_stack.pop()
            direction = path_directions.pop()

            # Update parent node's child reference
            if direction == 'left':
                # Balance current node using helper method
                current = current.balance()
                parent.left = current
            else:  # 'right'
                # Balance current node using helper method
                current = current.balance()
                parent.right = current

            # Update parent node
            self._update_node(parent)

            # Balance parent node
            current = parent.balance()

        return current if not parent_stack else node.balance()

    def insert_recursive(self, abs_position: int, donor_seq: str, donor_attrs: dict) -> None:
        """
        Recursively insert a donor sequence at the specified position.

        Args:
            abs_position (int): Absolute position for insertion
            donor_seq (str): Donor sequence to insert
            donor_attrs (dict): Attributes for the donor
        """
        self.root = self._insert_recursive(self.root, abs_position, donor_seq, donor_attrs)

    def _insert_recursive(self, node: SequenceNode, abs_position: int, donor_seq: str, donor_attrs: dict) -> SequenceNode:
        """
        Recursively insert a donor sequence at the absolute position in the tree.

        Args:
            node (SequenceNode): Current node
            abs_position (int): Absolute position in the tree to insert at
            donor_seq (str): Donor sequence to insert
            donor_attrs (dict): Attributes for the donor

        Returns:
            SequenceNode: New node after insertion
        """
        # Calculate positions in tree
        left_length = node.left.total_length if node.left else 0

        # Case 1: Position is in the current node's left subtree
        if abs_position <= left_length:
            if node.left:
                node.left = self._insert_recursive(node.left, abs_position, donor_seq, donor_attrs)
            else:
                # Insert as left child
                node.left = self._create_node(donor_seq, True, donor_attrs)
            # Use _update_node helper method
            self._update_node(node)
            return node.balance()

        # Case 2: Position is inside the current node
        node_start = left_length
        node_end = node_start + node.length
        if node_start < abs_position < node_end:
            # Split this node
            rel_pos = abs_position - node_start
            left_data = node.data[:rel_pos]
            right_data = node.data[rel_pos:]

            # Handle TSD generation if requested in attributes
            tsd_length = donor_attrs.get("tsd_length", 0)
            if tsd_length > 0:
                # Extract source TSD sequence from the original sequence
                # We extract from the right side of the split point (beginning of right_data)
                source_tsd_seq = right_data[:min(tsd_length, len(right_data))]

                # Generate TSD sequences (potentially with mutations)
                tsd_5, tsd_3 = generate_TSD(source_tsd_seq, tsd_length)

                # Remove source TSD from right_data as it will be duplicated
                if len(source_tsd_seq) > 0:
                    right_data = right_data[len(source_tsd_seq):]

                # Add TSD sequences to the left and right data
                left_data = left_data + tsd_5
                right_data = tsd_3 + right_data

            # Get current node attributes
            current_attrs = node.attrs.copy()
            current_is_donor = node.is_donor
            current_uid = node.uid

            # Create new left child with left data
            new_left = self._create_node(left_data, current_is_donor, current_attrs.copy())
            if node.left:
                new_left.left = node.left
                # Use _update_node helper method
                self._update_node(new_left)
                # Balance left subtree
                new_left = new_left.balance()

            # Create new right child with right data and original right child
            new_right = self._create_node(right_data, current_is_donor, current_attrs.copy())
            if node.right:
                new_right.right = node.right
                # Use _update_node helper method
                self._update_node(new_right)
                # Balance right subtree
                new_right = new_right.balance()

            # Record nested insertion relationships
            if current_is_donor:
                # Record that the donor was cut
                # Define the half attribute (left or right)
                new_left.attrs["half"] = ["L"] if "half" not in new_left.attrs else new_left.attrs["half"] + ["L"]
                new_right.attrs["half"] = ["R"] if "half" not in new_right.attrs else new_right.attrs["half"] + ["R"]

                # Prepare donor attrs to record that it's nested inside this donor
                donor_attrs_copy = donor_attrs.copy()
                if "nested_in" not in donor_attrs_copy:
                    donor_attrs_copy["nested_in"] = []
                donor_attrs_copy["nested_in"].append(current_uid)

                # Record which donor cut this node
                if "cut_by" not in new_left.attrs:
                    new_left.attrs["cut_by"] = []
                if "cut_by" not in new_right.attrs:
                    new_right.attrs["cut_by"] = []

                # Create a new donor node and get its UID
                new_donor_uid = self._get_next_uid()

                new_left.attrs["cut_by"].append(new_donor_uid)
                new_right.attrs["cut_by"].append(new_donor_uid)

                # Release current node's UID and set new values
                self._release_uid(node.uid)
                node.data = donor_seq
                node.length = len(donor_seq)
                node.is_donor = True
                node.attrs = donor_attrs_copy
                node.uid = new_donor_uid
                self.node_dict[new_donor_uid] = node
            else:
                # Replace this node's data with the donor sequence
                node.data = donor_seq
                node.length = len(donor_seq)
                node.is_donor = True
                node.attrs = donor_attrs.copy()

            # Set new children
            node.left = new_left
            node.right = new_right
            # Use _update_node helper method
            self._update_node(node)

            return node.balance()

        # Case 3: Position is in the current node's right subtree
        if abs_position >= node_end:
            if node.right:
                node.right = self._insert_recursive(node.right, abs_position - node_end, donor_seq, donor_attrs)
            else:
                # Insert as right child
                node.right = self._create_node(donor_seq, True, donor_attrs)
            # Use _update_node helper method
            self._update_node(node)
            return node.balance()

        raise RuntimeError("[ERROR] Should not reach here")

    def collect_donors(self, seq_id: str) -> Tuple[List[SeqRecord], List[SeqRecord]]:
        """
        Collect all donor nodes and reconstruct nested donors.

        Args:
            seq_id (str): ID of the original sequence

        Returns:
            Tuple[List[SeqRecord], List[SeqRecord]]: 
                - Regular donor records
                - Reconstructed donor records
        """
        donor_records, _ = self._collect_donors(self.root, seq_id)
        reconstructed_donors = self._reconstruct_donors(donor_records, seq_id)

        # Ensure we return two lists, even if reconstructed_donors is None
        if reconstructed_donors is None:
            reconstructed_donors = []

        return donor_records, reconstructed_donors

    def _collect_donors(self, node: SequenceNode, seq_id: str, abs_position: int = 0) -> Tuple[List[SeqRecord], int]:
        """
        Recursively collect donor nodes from the tree.

        Args:
            node (SequenceNode): Current node
            seq_id (str): ID of the original sequence
            abs_position (int): Current absolute position

        Returns:
            Tuple[List[SeqRecord], int]: Donor records and total length
        """
        if not node:
            return [], 0

        donor_records = []

        # Process left subtree
        left_donors, left_length = self._collect_donors(node.left, seq_id, abs_position)
        donor_records.extend(left_donors)

        # Process current node
        current_position = abs_position + left_length
        if node.is_donor:
            # Create a donor record with absolute position
            start_index = current_position
            end_index = start_index + node.length

            # Add information about nesting in the ID
            nested_info = ""
            if "nested_in" in node.attrs and node.attrs["nested_in"]:
                nested_info = f"-nested_in_{'_'.join(map(str, node.attrs['nested_in']))}"

            donor_id = f"{seq_id}_{start_index}_{end_index}-+-{node.length}{nested_info}"
            donor_record = create_sequence_record(node.data, donor_id)

            # Add the node's UID and other attributes as annotations
            donor_record.annotations["uid"] = node.uid

            if "cut_by" in node.attrs and node.attrs["cut_by"]:
                donor_record.annotations["cut_by"] = node.attrs["cut_by"]

            if "half" in node.attrs and node.attrs["half"]:
                donor_record.annotations["half"] = node.attrs["half"]

            donor_records.append(donor_record)

        # Process right subtree
        right_donors, right_length = self._collect_donors(node.right, seq_id, current_position + node.length)
        donor_records.extend(right_donors)

        # Return all donors and total length
        return donor_records, left_length + node.length + right_length

    @staticmethod
    def _reconstruct_donors(donor_records: List[SeqRecord], seq_id: str) -> List[SeqRecord]:
        """
        Reconstruct nested donors from collected donor records.

        Args:
            donor_records (List[SeqRecord]): Collected donor records
            seq_id (str): ID of the original sequence

        Returns:
            List[SeqRecord]: Reconstructed donor records
        """
        # Create a dictionary to organize donors by their UID
        donor_dict = {}
        for record in donor_records:
            if "uid" in record.annotations:
                donor_dict[record.annotations["uid"]] = record

        # Create an index to efficiently find matching halves
        half_index = {}
        for record in donor_records:
            if ("uid" not in record.annotations or 
                "cut_by" not in record.annotations or 
                "half" not in record.annotations):
                continue

            uid = record.annotations["uid"]
            cut_by_list = record.annotations["cut_by"]
            half_list = record.annotations["half"]

            # Ensure cut_by_list and half_list have matching lengths
            if len(cut_by_list) != len(half_list):
                print(f"Warning: Record {record.id} has mismatched cut_by and half lists, skipping")
                continue

            # Build index for each cutting relationship
            for cut_by_uid, half in zip(cut_by_list, half_list):
                key = (cut_by_uid, half)
                if key not in half_index:
                    half_index[key] = []
                half_index[key].append(uid)

        # Reconstruct cut donors
        reconstructed_donors = []
        processed_uids = set()

        for record in donor_records:
            if ("uid" not in record.annotations or 
                "cut_by" not in record.annotations or 
                "half" not in record.annotations or 
                record.annotations["uid"] in processed_uids):
                continue

            uid = record.annotations["uid"]
            cut_by_list = record.annotations["cut_by"]
            half_list = record.annotations["half"]

            # Ensure cut_by_list and half_list have matching lengths
            if len(cut_by_list) != len(half_list):
                continue  # Already checked above

            # Process each cutting relationship
            for cut_by_uid, half in zip(cut_by_list, half_list):
                # Use the index to quickly find matching halves
                opposite_half = "R" if half == "L" else "L"
                opposite_key = (cut_by_uid, opposite_half)

                if opposite_key in half_index and half_index[opposite_key]:
                    # Found matching halves
                    for matching_uid in half_index[opposite_key]:
                        if matching_uid in processed_uids or matching_uid == uid:
                            continue

                        matching_record = donor_dict.get(matching_uid)
                        if not matching_record:
                            continue

                        # Get cutting donor
                        cutting_donor = donor_dict.get(cut_by_uid)
                        if cutting_donor:
                            cutting_seq = str(cutting_donor.seq)
                        else:
                            # If cutting donor is missing, mark it
                            cutting_seq = "[MISSING_CUTTING_DONOR]"
                            print(f"Warning: Could not find cutting donor with UID {cut_by_uid}")

                        # Determine which half is left and which is right
                        left_record = record if half == "L" else matching_record
                        right_record = matching_record if half == "L" else record

                        # Create reconstructed donor sequence (with cutting donor)
                        reconstructed_seq = ""
                        # First add the left fragment
                        reconstructed_seq += str(left_record.seq)
                        # Add the cutting donor
                        reconstructed_seq += cutting_seq
                        # Add the right fragment
                        reconstructed_seq += str(right_record.seq)

                        # Create record for the reconstructed donor
                        rec_id = f"{seq_id}_reconstructed_cut_donor_{uid}_{cut_by_uid}"
                        rec_record = create_sequence_record(reconstructed_seq, rec_id)
                        rec_record.annotations["original_fragments"] = [uid, matching_uid]
                        rec_record.annotations["cutting_donor"] = cut_by_uid

                        reconstructed_donors.append(rec_record)

                        # Also create a clean reconstructed donor (without cutting donor)
                        clean_reconstructed_seq = str(left_record.seq) + str(right_record.seq)
                        clean_rec_id = f"{seq_id}_clean_reconstructed_donor_{uid}_{cut_by_uid}"
                        clean_rec_record = create_sequence_record(clean_reconstructed_seq, clean_rec_id)
                        clean_rec_record.annotations["original_fragments"] = [uid, matching_uid]

                        reconstructed_donors.append(clean_rec_record)

                        # Mark both fragments as processed
                        processed_uids.add(uid)
                        processed_uids.add(matching_uid)

                        # Successfully found a match, proceed to next record
                        break

        return reconstructed_donors

    def to_graphviz_dot(self, node_id_prefix: str = "node") -> str:
        """
        Generate a Graphviz DOT representation of the tree structure for visualization.

        Args:
            node_id_prefix (str): Prefix for node IDs in the graph

        Returns:
            str: Graphviz DOT format string
        """
        if not self.root:
            return "digraph SequenceTree { }"

        # Initialize the DOT string with graph declaration
        dot_str = ["digraph SequenceTree {",
                   "  bgcolor=\"#FFF\"",
                   "  node [fontcolor=\"#000\", shape=box, style=filled];",
                   "  edge [fontcolor=\"#000\"];"]

        # Generate nodes and edges through recursive traversal
        nodes, edges = self._build_graphviz_dot_nodes_edges(self.root, node_id_prefix)

        # Add all nodes and edges to the DOT string
        for node in nodes:
            dot_str.append(f"  {node}")
        for edge in edges:
            dot_str.append(f"  {edge}")

        dot_str.append('}')
        return '\n'.join(dot_str)

    def _build_graphviz_dot_nodes_edges(self, node: SequenceNode, node_id_prefix: str, abs_pos: int = 0) -> tuple:
        """
        Recursively build nodes and edges for Graphviz visualization.

        Args:
            node (SequenceNode): Current node
            node_id_prefix (str): Prefix for node IDs
            abs_pos (int): Current absolute position

        Returns:
            tuple: (nodes list, edges list)
        """
        if not node:
            return [], []

        nodes = []
        edges = []

        # Calculate positions
        left_length = node.left.total_length if node.left else 0

        # Calculate start and end positions
        start_pos = abs_pos + left_length
        end_pos = start_pos + node.length

        # Generate unique ID for this node
        node_id = f"{node_id_prefix}_{start_pos}_{end_pos}"

        # Determine node type and color
        node_type = "Donor" if node.is_donor else "Acceptor"
        fill_color = "lightblue" if node.is_donor else "lightgreen"
        nest_cut = node.attrs.get("nested_in", []) or node.attrs.get("half", []) or []
        additional_info = ','.join(nest_cut) if nest_cut else ""
        if nest_cut:
            if {"L", "R"} & set(nest_cut):  # node is cut by cutting donor
                fill_color = "lightpink"
            else:
                fill_color = "yellow"

        # Create node label with position information
        label = "".join([node_type, " | ", str(node.uid), "\\n",
                         str(start_pos), "\\l",
                         str(end_pos), "\\l",
                         "Length: ", str(node.length), "\\n",
                         additional_info])

        # Add the node to the nodes list
        nodes.append(f'{node_id} [label="{label}", fillcolor="{fill_color}"];')

        # Process left child if exists
        if node.left:
            # Left child should start at the same absolute position as its parent
            left_abs_pos = abs_pos
            left_nodes, left_edges = self._build_graphviz_dot_nodes_edges(
                node.left, f"{node_id_prefix}_L", left_abs_pos
            )
            nodes.extend(left_nodes)
            edges.extend(left_edges)

            # Get the first node ID from left nodes (should be the root of left subtree)
            # If it's empty (unlikely), create a placeholder ID
            left_id = None
            for left_node in left_nodes:
                # Extract node ID from the DOT node definition
                match = re.match(r'([^\s\[]+)\s*\[', left_node)
                if match:
                    left_id = match.group(1)
                    break

            # If we couldn't find a proper ID, create one based on our knowledge of the left child
            if not left_id:
                left_start = left_abs_pos + (node.left.left.total_length if node.left.left else 0)
                left_end = left_start + node.left.length
                left_id = f"{node_id_prefix}_L_{left_start}_{left_end}"

            # Add edge from this node to left child
            edges.append(f'{node_id} -> {left_id} [label="L"];')

        # Process right child if exists
        if node.right:
            # Right child starts at the end position of the current node
            right_abs_pos = abs_pos + left_length + node.length
            right_nodes, right_edges = self._build_graphviz_dot_nodes_edges(
                node.right, f"{node_id_prefix}_R", right_abs_pos
            )
            nodes.extend(right_nodes)
            edges.extend(right_edges)

            # Get the first node ID from right nodes (should be the root of right subtree)
            # If it's empty (unlikely), create a placeholder ID
            right_id = None
            for right_node in right_nodes:
                # Extract node ID from the DOT node definition
                match = re.match(r'([^\s\[]+)\s*\[', right_node)
                if match:
                    right_id = match.group(1)
                    break

            # If we couldn't find a proper ID, create one based on our knowledge of the right child
            if not right_id:
                right_start = right_abs_pos + (node.right.left.total_length if node.right.left else 0)
                right_end = right_start + node.right.length
                right_id = f"{node_id_prefix}_R_{right_start}_{right_end}"

            # Add edge from this node to right child
            edges.append(f'{node_id} -> {right_id} [label="R"];')

        return nodes, edges

    @staticmethod
    def _update_node(node: SequenceNode) -> None:
        """
        Helper method to update a node's height and total length.
        Correctly calls the node's internal methods.

        Args:
            node (SequenceNode): The node to update
        """
        node.update_height()
        node.update_total_length()

# ======================================================================================================================

def convert_humanized_number(text: str, base: Union[float, int], units: Sequence[Union[str, Iterable[str]]]) -> float:
    """
    Convert human-readable numbers with unit suffixes to numerical values.

    Parses strings containing numbers with optional unit suffixes and converts them
    to raw numerical values based on the provided base and unit definitions.

    :param text: Input string containing the number and optional unit suffix.
                Supports scientific notation.
    :type text: str
    :param base: Numerical base for unit conversion (e.g., 1000 or 1024).
    :type base: Union[float, int]
    :param units: Sequence of unit suffixes/aliases ordered by increasing scale.
                  Can contain strings or iterables of aliases for each unit.
                  Must start from the base unit. Use empty string for unitless case.
    :type units: Sequence[Union[str, Iterable[str]]]
    :return: The converted numerical value.
    :rtype: float
    :raises ValueError: If the input format is invalid or contains unrecognized units.

    Examples:
        >>> convert_humanized_number("1.14 kbp", 1000, ["bp", ("kbp", "kb"), ("Mbp", "Mb")])
        1140.0
        >>> convert_humanized_number("5.14", 1024, ('B', "KiB", "MiB", "GiB"))
        5.14
        >>> convert_humanized_number("1.919e-3M", 1000, ("", 'K', 'M', 'B', 'T'))
        1919.0
    """
    match = re.fullmatch(
        r"\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*([a-zA-Z]*)",
        text.strip()
    )
    if not match:
        raise ValueError(f"Invalid format: \"{text}\"")

    num_str, raw_unit = match.groups()

    try:
        number = float(num_str)
    except ValueError:
        raise ValueError(f"Invalid number: \"{num_str}\"")

    if not raw_unit:
        return number

    target_unit = raw_unit.lower()
    try:
        unit_index = next(index for index, entry in enumerate(units)
                          if any(unit.lower() == target_unit for unit in
                                 (entry if isinstance(entry, Iterable) and
                                           not isinstance(entry, str) else (entry,))))
    except StopIteration:
        raise ValueError(f"Unknown unit \"{raw_unit}\"")

    return number * (base ** unit_index)

def convert_humanized_int(number: Union[str, int]) -> int:
    if isinstance(number, int):
        return number
    return int(convert_humanized_number(number, 1000, ("", 'K', 'M', 'B', 'T')))


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
    Create a BioPython SeqRecord object from a sequence string and ID.

    This function takes a biological sequence string and an identifier,
    and returns a SeqRecord object from the BioPython library.

    :param seq: The biological sequence as a string
    :type seq: str
    :param id: The unique identifier for the sequence
    :type id: str
    :return: A SeqRecord object containing the sequence and its metadata
    :rtype: SeqRecord

    Example::

        >>> record = create_sequence_record("ATGC", "seq1")
        >>> print(record.id)
        seq1
    """
    return SeqRecord(Seq(seq), id=id, description="")

def save_multi_fasta_from_dict(records_dict: Dict[str, List[SeqRecord]], output_dir_path: str):
    """
    Save multiple FASTA files from a dictionary of sequence records.

    Args:
        records_dict (Dict[str, List[SeqRecord]]): Dictionary mapping file names to lists of SeqRecord objects
        output_dir_path (str): Directory path where the output files will be saved
    """
    if not records_dict:
        return

    os.makedirs(output_dir_path, exist_ok=True)
    for file_name, records in records_dict.items():
        if records:  # Skip empty record lists
            output_file_path: str = os.path.join(output_dir_path, file_name)
            SeqIO.write(records, output_file_path, "fasta")

def load_sequences(path_list: Optional[List[str]], len_limit: Optional[int] = None,
                   flag_filter_n: bool = False) -> Optional[List[str]]:
    """
    Load donor sequences from donor library files.
    Only loads sequences that meet the length criteria into memory.
    Sequences are stored as strings for efficient processing.

    Args:
        path_list: List of paths to donor sequence files
        len_limit: Optional maximum length limit for sequences to load
        flag_filter_n: Flag to filter out sequences containing N

    Returns:
        Optional[List[str]]: List of donor sequences as strings, or None if no donor paths provided
    """
    if not path_list:
        return None

    donor_sequences: List[str] = []
    for file_path in path_list:
        donor_sequences.extend([str(record.seq.upper()) for record in SeqIO.parse(file_path, "fasta") if
                              (len(record) <= len_limit if len_limit else True) and
                              ("N" not in record.seq.upper() if flag_filter_n else True)])

    donor_sequences.sort(key=len)
    return donor_sequences

def _find_donor_lib_abs_path_list(path: Optional[str]) -> Optional[List[str]]:
    """
    Locate the donor library files.

    This function attempts to find the given donor library directory or file in two locations:
    1. The directly specified path
    2. Under the default "lib" directory in the program root directory

    Args:
        path (str): The path to search for donor library. Can be:
                        - An absolute path to a directory or file
                        - A relative path to a directory or file

    Returns:
        Optional[List[str]]: A list containing the absolute paths to the found donor library files,
                             or None if path is None.

    Raises:
        SystemExit: If no valid donor library directory or file is found.

    Examples:
        >>> _find_donor_lib_abs_path_list("/abs/path/to/donors.fa")
        ["/abs/path/to/donors.fa"]
        >>> _find_donor_lib_abs_path_list("lib/rice")
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
        donor_lib_files: List[str] = []
        for file in os.listdir(path):
            if file.endswith((".fa", ".fa")):
                donor_lib_files.append(os.path.join(path, file))
        if donor_lib_files:
            return list(map(os.path.abspath, donor_lib_files))

    # Check in built-in DonorLib directory
    default_donor_lib_abs_path = os.path.join(PROGRAM_ROOT_DIR_ABS_PATH, path)
    if os.path.exists(default_donor_lib_abs_path):
        if os.path.isfile(default_donor_lib_abs_path):
            return [default_donor_lib_abs_path]
        elif os.path.isdir(default_donor_lib_abs_path):
            donor_lib_files: List[str] = []
            for file in os.listdir(default_donor_lib_abs_path):
                if file.endswith(".fa"):
                    donor_lib_files.append(os.path.join(default_donor_lib_abs_path, file))
            if donor_lib_files:
                return donor_lib_files

    raise SystemExit(f"[ERROR] Specified donor library not found at {path} or built-in donor library!")

def _load_multiple_donor_libs(path_list: List[str], weight_list: Optional[List[float]] = None,
                            len_limit: Optional[int] = None,
                            flag_filter_n: bool = False) -> Optional[Tuple[List[str], List[int], List[float]]]:
    """
    Load multiple donor libraries by combining their sequences.

    Args:
        path_list (List[str]): List of paths to donor libraries
        weight_list (Optional[List[float]]): List of weights for each donor library
        len_limit (Optional[int]): Maximum length limit for sequences to load
        flag_filter_n (bool): Flag to filter out sequences containing N

    Returns:
        Tuple[List[str], List[int], List[float]] (Optional):
            - List[str]: Combined list of donor sequences from all libraries
            - List[int]: List of lengths of donor sequences
            - List[float]: List of weights for each donor sequence
    """
    if not path_list:
        return None

    all_donor_sequence_list: List[str] = []
    all_donor_len_list: List[int] = []
    all_donor_weight_list: List[float] = []

    for path_list, weight in zip(path_list, weight_list or [None] * len(path_list)):
        single_donor_lib_abs_path_list = _find_donor_lib_abs_path_list(path_list)
        if single_donor_lib_abs_path_list:
            single_donor_lib_sequences = load_sequences(single_donor_lib_abs_path_list, len_limit, flag_filter_n)
            if single_donor_lib_sequences:
                all_donor_sequence_list.extend(single_donor_lib_sequences)
                all_donor_len_list.extend([len(seq) for seq in single_donor_lib_sequences])
                if weight:
                    len_single_donor_lib_sequences = len(single_donor_lib_sequences)
                    all_donor_weight_list.extend([weight / len_single_donor_lib_sequences] * len_single_donor_lib_sequences)

    # return sort_multiple_lists(all_donor_sequence_list, all_donor_weight_list, key=len)
    all_donor_len_list, all_donor_sequence_list, all_donor_weight_list = sort_multiple_lists(
        all_donor_len_list, all_donor_sequence_list, all_donor_weight_list)

    # Uniform weights by default: If no weight provided, set all weights to 1
    if not weight_list:
        all_donor_weight_list = [1] * len(all_donor_sequence_list)

    return all_donor_sequence_list, all_donor_len_list, all_donor_weight_list

class SeqGenerator:
    def __init__(self, input_file: str, insertion: Union[str, int], batch: int, processors: int, output_dir_path: str,
                 donor_lib: Optional[List[str]] = None, donor_lib_weight: Optional[List[float]] = None,
                 donor_len_limit: Optional[int] = None, flag_filter_n: bool = False, flag_track: bool = False,
                 tsd_length: Optional[int] = None, flag_visual: bool = False, flag_recursive: bool = False):
        """
        Initialize the sequence generator.

        Args:
            input_file (str): Path to the input sequence file
            insertion (Union[str, int]): Number of insertions per sequence
            batch (int): Number of independent result files to generate
            processors (int): Number of processors to use
            output_dir_path (str): Directory to save output files
            donor_lib (List[str], optional): List of donor library file paths
            donor_lib_weight (List[float], optional): Weights for donor libraries
            donor_len_limit (int, optional): Maximum length limit for donor sequences to load
            flag_filter_n (bool): Whether to filter out sequences containing N
            flag_track (bool): Whether to track donor sequences used
            tsd_length (int, optional): Length of Target Site Duplication (TSD) to generate
            flag_visual (bool): Whether to generate graphviz visualization of the sequence tree
            flag_recursive (bool): Whether to use recursive insertion method
        """
        self.input_file: str = input_file
        self.insertion: int = convert_humanized_int(insertion)
        self.batch: int = batch
        self.processors: int = processors
        self.output_dir_path: str = output_dir_path
        self.donor_len_limit: Optional[int] = donor_len_limit
        self.flag_filter_n: bool = flag_filter_n
        self.flag_track: bool = flag_track
        self.tsd_length: Optional[int] = tsd_length
        self.flag_visual: bool = flag_visual
        self.flag_recursive: bool = flag_recursive

        # Load input sequence
        try:
            self.input: List[SeqRecord] = list(SeqIO.parse(input_file, "fasta"))
            if not self.input:
                raise ValueError(f"No sequences found in input file {input_file}")
        except (FileNotFoundError, IOError) as e:
            raise ValueError(f"Error loading input file {input_file}: {str(e)}")

        # Load donor sequences
        if not donor_lib:
            raise ValueError("Donor library is required for insertion. Please provide at least one donor library path using the --donor parameter.")

        try:
            donor_lib_data = _load_multiple_donor_libs(donor_lib, donor_lib_weight, donor_len_limit, flag_filter_n)
            if not donor_lib_data or not donor_lib_data[0]:
                raise ValueError("No valid donor sequences loaded. Check if donor files exist and contain valid sequences.")
            # Assign donor sequences and weights, ignore lengths as they're not needed
            self.donor_sequences, _, self.donor_weights = donor_lib_data
            # Normalize weights if they were provided
            if donor_lib_weight and sum(self.donor_weights) > 0:
                weight_sum = sum(self.donor_weights)
                self.donor_weights = [w / weight_sum for w in self.donor_weights]
        except Exception as e:
            # Convert any exception during donor loading to a more informative error
            raise ValueError(f"Error loading donor libraries: {str(e)}")

    def __pre_check(self):
        """
        Perform pre-execution checks.
        """
        if not os.path.exists(self.output_dir_path):
            os.makedirs(self.output_dir_path)

    def __process_batch_multiprocessing(self, total_sequences: int) -> Tuple[List[SeqRecord], 
                                                                         Optional[List[SeqRecord]],
                                                                         Optional[List[SeqRecord]]]:
        """
        Process a batch of sequences using multiprocessing.

        Args:
            total_sequences (int): Total number of sequences to process

        Returns:
            Tuple[List[SeqRecord], Optional[List[SeqRecord]], Optional[List[SeqRecord]]]: 
                All processed sequences, donor records, and reconstructed donor records
        """
        if total_sequences == 0:
            return [], None, None

        item_indices = list(range(total_sequences))

        print(f"Using multiprocessing with {self.processors} worker processes to process {total_sequences} sequences")

        all_sequences = []
        all_donors = [] if self.flag_track else None
        all_reconstructed_donors = [] if self.flag_track else None

        with mp.Pool(self.processors) as pool:
            results = pool.map(self._process_single_sequence, item_indices)

            for i, (processed_seq, donors, reconstructed) in enumerate(results, 1):
                all_sequences.append(processed_seq)
                if self.flag_track:
                    if donors:
                        all_donors.extend(donors)
                    if reconstructed:
                        all_reconstructed_donors.extend(reconstructed)

                if i % 10 == 0 or i == total_sequences:
                    print(f"Completed {i}/{total_sequences} sequences")

        return all_sequences, all_donors, all_reconstructed_donors

    def _process_single_sequence(self, idx_or_record) -> Tuple[SeqRecord, Optional[List[SeqRecord]], Optional[List[SeqRecord]]]:
        """
        Process a single sequence with random insertions.

        Args:
            idx_or_record: Either an index (int) to look up in self.input, or a SeqRecord to process directly

        Returns:
            Tuple[SeqRecord, Optional[List[SeqRecord]], Optional[List[SeqRecord]]]: 
                The processed sequence, used donors, and reconstructed donors if tracking
        """
        # If idx_or_record is an integer, retrieve the sequence record
        # If it's already a SeqRecord, use it directly
        if isinstance(idx_or_record, int):
            # Set random seed for multiprocessing
            random.seed(os.getpid() + int(time.time()))
            seq_record = self.input[idx_or_record]
        else:
            seq_record = idx_or_record

        # Check if there are donor sequences to insert
        if not self.insertion or not self.donor_sequences:
            # If no insertions requested or no donor sequences available,
            # return the original sequence with modified ID
            # May be useful for future random sequence generation function
            new_id = f"{seq_record.id}_ins0"
            return create_sequence_record(str(seq_record.seq), new_id), None, None

        # Create a new sequence tree for this sequence
        seq_tree = SequenceTree(str(seq_record.seq), 0)

        # Get the current total sequence length
        total_length = len(str(seq_record.seq))

        # Generate insertion positions and donor sequences
        insert_positions = random.choices(range(total_length + 1), k=self.insertion)
        selected_donors = random.choices(self.donor_sequences, weights=self.donor_weights, k=self.insertion)

        # Sort positions in descending order to maintain position integrity during insertion
        insert_positions, selected_donors = sort_multiple_lists(insert_positions, selected_donors, reverse=True)

        # Insert donors into the tree
        for pos, donor_seq in zip(insert_positions, selected_donors):
            attrs = {"tsd_length": self.tsd_length} if self.tsd_length else {}
            if self.flag_recursive:
                seq_tree.insert_recursive(pos, donor_seq, attrs)
            else:
                seq_tree.insert(pos, donor_seq, attrs)

        # Generate final sequence
        new_id = f"{seq_record.id}_ins{self.insertion}"
        if self.tsd_length:
            new_id += f"_tsd{self.tsd_length}"
        new_seq_record = create_sequence_record(str(seq_tree), new_id)

        # Generate visualization if requested
        if self.flag_visual:
            graphviz_str = seq_tree.to_graphviz_dot(node_id_prefix=seq_record.id)
            visual_dir_path = os.path.join(self.output_dir_path, TREE_VISUAL_DIR_NAME)
            os.makedirs(visual_dir_path, exist_ok=True)
            with open(os.path.join(visual_dir_path, f"{seq_record.id}_tree_visual.dot"), "w") as f:
                f.write(graphviz_str)
            print(f"Generated tree visualization for {seq_record.id}")

        # Collect donor sequences if tracking is enabled
        used_donors = None
        reconstructed_donors = None
        if self.flag_track:
            used_donors, reconstructed_donors = seq_tree.collect_donors(seq_record.id)

        return new_seq_record, used_donors, reconstructed_donors

    def execute(self):
        """
        Execute the sequence generation process.
        """
        self.__print_header()

        start_time = time.time()
        self.__pre_check()

        # Process each batch
        for batch_num in range(1, self.batch + 1):
            self.__process_single_batch(batch_num)

        # Print summary
        total_elapsed_time = time.time() - start_time
        self.__print_summary(total_elapsed_time)

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
        print(f"Donor library: {len(self.donor_sequences)} sequences loaded")
        print(f"Generating {self.batch} independent result file(s)")
        if self.flag_recursive:
            print(f"Using recursive insertion method")
        if self.flag_visual:
            print(f"Tree visualization enabled: DOT files will be generated for each sequence")

    def __print_summary(self, total_elapsed_time: float):
        """
        Print a summary of the processing results.

        Args:
            total_elapsed_time (float): Total time taken for all batches
        """
        print(f"\nAll batches completed in {total_elapsed_time:.2g} seconds")
        print(f"Results saved to \"{os.path.abspath(self.output_dir_path)}\"")

    def __process_single_batch(self, batch_num: int) -> float:
        """
        Process a single batch of sequences.

        Args:
            batch_num (int): Batch number (starting from 1)

        Returns:
            float: Time taken to process this batch (seconds)
        """
        batch_start_time = time.time()

        # Set random seed for this batch
        random.seed(int(time.time()) + batch_num)

        print(f"\nProcessing batch {batch_num}/{self.batch}")

        total_sequences = len(self.input)

        # Choose processing method
        all_sequences, all_donors, all_reconstructed_donors = self.__process_batch_multiprocessing(total_sequences)

        os.makedirs(self.output_dir_path, exist_ok=True)
        print(f"Output directory: {self.output_dir_path}")

        # Save results
        self.__save_batch_results(self.output_dir_path, all_sequences, all_donors, all_reconstructed_donors, batch_num)

        batch_elapsed_time = time.time() - batch_start_time
        print(f"Batch {batch_num} completed in {batch_elapsed_time:.2g} seconds")

        return batch_elapsed_time

    def __save_batch_results(self, output_dir: str, sequences: List[SeqRecord], 
                             donors: Optional[List[SeqRecord]], 
                             reconstructed_donors: Optional[List[SeqRecord]],
                             batch_num: int = 1):
        """
        Save the batch results to output files.

        Args:
            output_dir (str): Directory to save the results
            sequences (List[SeqRecord]): List of sequence records to save
            donors (Optional[List[SeqRecord]]): List of donor records to save if tracking is enabled
            reconstructed_donors (Optional[List[SeqRecord]]): List of reconstructed donor records
            batch_num (int): The batch number (1-based index)
        """
        if not sequences:
            print("No sequences to save")
            return

        # Add batch suffix to filenames if multiple batches
        suffix = f"_batch{batch_num}" if self.batch > 1 else ""

        print(f"Saving {len(sequences)} processed sequences")
        output_dict = {f"sequences{suffix}.fa": sequences}

        if self.flag_track and donors:
            print(f"Saving {len(donors)} donor sequence records")
            output_dict[f"used_donors{suffix}.fa"] = donors

            if reconstructed_donors:
                # Group by reconstruction type
                full_recs = [d for d in reconstructed_donors if "reconstructed_cut_donor" in d.id]
                clean_recs = [d for d in reconstructed_donors if "clean_reconstructed_donor" in d.id]

                if full_recs:
                    print(f"Saving {len(full_recs)} full reconstructed donor records")
                    output_dict[f"reconstructed_donors{suffix}.fa"] = full_recs

                if clean_recs:
                    print(f"Saving {len(clean_recs)} clean reconstructed donor records")
                    output_dict[f"cut_donors{suffix}.fa"] = clean_recs

        save_multi_fasta_from_dict(output_dict, output_dir)

def main():
    """Main function to handle command line arguments and execute the program."""
    parser = argparse.ArgumentParser(prog="RandSeqInsert",
                                     description="RandomSequenceInsertion: A high-performance Python tool for randomly inserting "
                                               "genomic fragments from donor libraries into acceptor sequences. "
                                               "It supports multiprocessing for efficient processing of large datasets.",
                                     epilog="")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {VERSION}\n{INFO}")

    # Core Arguments
    core_group = parser.add_argument_group("Core Arguments")
    core_group.add_argument("-i", "--input",
                       help="Input sequence file in FASTA format. Contains the acceptor sequences to be inserted into.",
                       type=str, required=True, metavar="FILE")
    core_group.add_argument("-o", "--output", default=DEFAULT_OUTPUT_DIR_ABS_PATH, metavar="DIR",
                       help=f"Output directory path. Generated sequences and related files will be saved here. Default: \"{DEFAULT_OUTPUT_DIR_ABS_PATH}\"")

    core_group.add_argument("-is", "--insert", metavar="INT/STR",
                       help="Number of insertions per sequence. Accepts either a plain number (e.g., 100) or a string with k/m suffix (e.g., 1k, 1m).",
                       type=str, required=True)

    # Donor Library Arguments
    donor_group = parser.add_argument_group("Donor Library Arguments")
    donor_group.add_argument("-d", "--donor", nargs="+", metavar="FILE/DIR", required=True,
                            help="Donor sequence library file or directory paths. Multiple FASTA format donor files can be specified. Sequences from these files will be selected and inserted into the acceptor sequences.")
    donor_group.add_argument("-w", "--weight", type=float, nargs="+", metavar="FLOAT",
                       help="Weights for donor libraries. Controls the probability of selecting sequences from different donor libraries. The number of weights should match the number of donor libraries.")
    donor_group.add_argument("-l", "--limit", type=int, default=None, metavar="INT",
                       help="Donor sequence length limit. Only loads donor sequences with length less than or equal to this value. Default: no limit.")

    # Control Arguments
    ctrl_group = parser.add_argument_group("Control Arguments")
    ctrl_group.add_argument("-b", "--batch", type=int, default=1, metavar="INT",
                            help="Number of independent result files to generate. Runs the entire process multiple times with different random seeds to generate multiple output sets. Default: 1")
    ctrl_group.add_argument("-p", "--processors", type=int, default=DEFAULT_ALLOCATED_CPU_CORES, metavar="INT",
                       help=f"Number of processors to use for parallel processing. Default: {DEFAULT_ALLOCATED_CPU_CORES}")

    # Flags
    flag_group = parser.add_argument_group("Flags")
    flag_group.add_argument("--filter_n", action="store_true",
                                help="Filter out donor sequences containing N.")
    flag_group.add_argument("--track", action="store_true",
                       help="Track and save used donor sequences. Enable this option to generate an additional FASTA file in the output directory recording all used donor sequences and their insertion positions.")
    flag_group.add_argument("--tsd", type=int, metavar="LENGTH",
                       help="Enable Target Site Duplication (TSD) with specified length. When a donor sequence is inserted, TSD of this length will be generated at the insertion site.")
    flag_group.add_argument("--visual", action="store_true",
                       help="Generate Graphviz DOT files visualizing the tree structure of each sequence. Files will be named {seqid}_tree_visual.dot and saved in the output directory.")
    flag_group.add_argument("--recursive", action="store_true",
                       help="Use recursive insertion method instead of iterative one.")

    parsed_args = parser.parse_args()

    generator = SeqGenerator(
        input_file=parsed_args.input,
        insertion=parsed_args.insert,
        batch=parsed_args.batch,
        processors=parsed_args.processors,
        output_dir_path=parsed_args.output,
        donor_lib=parsed_args.donor,
        donor_lib_weight=parsed_args.weight,
        donor_len_limit=parsed_args.limit,
        flag_filter_n=parsed_args.filter_n,
        flag_track=parsed_args.track,
        tsd_length=parsed_args.tsd,
        flag_visual=parsed_args.visual,
        flag_recursive=parsed_args.recursive
    )
    generator.execute()

if __name__ == "__main__":
    main()
