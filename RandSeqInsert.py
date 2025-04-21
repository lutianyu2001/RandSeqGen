#!/usr/bin/env python3
"""
Random Sequence Insertion Generator

This program takes acceptor sequences and performs random insertions using sequences
from a donor library. It supports multiprocessing for efficient sequence generation.

Note: This program uses 1-based indexing for all sequence positions (both input and output).

Author: Tianyu Lu (tianyu@lu.fm)
Date: 2024-11-27
"""

import os
import argparse
import random
import re
from typing import Tuple, List, Dict, Iterable, Sequence, Union, Optional, Callable, Set
import multiprocessing as mp
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

VERSION = "v1.1.0"
INFO = ("by Tianyu (Sky) Lu (tianyu@lu.fm) "
        "released under GPLv3")

PROGRAM_ROOT_DIR_ABS_PATH = os.path.dirname(__file__)

# PRE-DEFINED PARAMETERS
DEFAULT_TSD_SNP_MUTATION_RATE = 0.05
DEFAULT_TSD_INDEL_MUTATION_RATE = 0.05

VISUAL_DIR_NAME = "visualization"
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

def generate_TSD(seq_slice_right: str,
                 length: Optional[int] = None,
                 snp_mutation_rate: float = DEFAULT_TSD_SNP_MUTATION_RATE,
                 indel_mutation_rate: float = DEFAULT_TSD_INDEL_MUTATION_RATE) -> Tuple[str, str]:
    """Generate TSD sequences with potential mutations.

    Args:
        seq_slice_right (str): sequence fragment on right side (3' direction)
        length (int): TSD length (defaults to seq_slice_right length)
        snp_mutation_rate (float): SNP mutation probability
        indel_mutation_rate (float): InDel mutation probability

    Returns:
        Tuple[str, str]: (5' TSD sequence, 3' TSD sequence)
    """
    # Initialize both ends of TSD
    tsd_length = min(len(seq_slice_right), length or len(seq_slice_right))
    tsd_5 = tsd_3 = seq_slice_right[:tsd_length]

    # Apply SNP mutation
    if random.random() < snp_mutation_rate:
        if random.random() < 0.5:
            tsd_5 = add_snp_mutation(tsd_5)
        else:
            tsd_3 = add_snp_mutation(tsd_3)

    # Apply InDel mutation
    if random.random() < indel_mutation_rate:
        if random.random() < 0.5:
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
    def __init__(self, data: str, is_donor: bool = False, donor_id: str = None, uid: int = None):
        """
        Initialize a sequence node.

        Args:
            data (str): The sequence string
            is_donor (bool): Whether this node contains a donor sequence
            donor_id (str): Identifier for the donor sequence (if is_donor is True)
            uid (int): Unique identifier for this node
        """
        self.data = data
        self.length = len(data)
        self.is_donor = is_donor
        self.donor_id = donor_id
        self.left = None
        self.right = None
        self.uid = uid

        # Total length of the subtree for efficient traversal
        self.total_length = self.length
        # Height of the node for AVL balancing
        self.height = 1

    def __iter__(self):
        """
        Implement in-order traversal of the tree, iterating all nodes in left-root-right order
        Yields:
            SequenceNode: Each node in the tree
        """
        yield from SequenceNode._inorder_traversal(self)

    @staticmethod
    def _inorder_traversal(node):
        if not node:
            return
        yield from SequenceNode._inorder_traversal(node.left)
        yield node
        yield from SequenceNode._inorder_traversal(node.right)

    def __str__(self) -> str:
        """
        Convert the tree to a string by in-order traversal.

        Returns:
            str: The concatenated sequence
        """
        return "".join([node.data for node in self])

    def update_height(self):
        """
        Update the height of this node.
        """
        left_height = self.left.height if self.left else 0
        right_height = self.right.height if self.right else 0
        self.height = max(left_height, right_height) + 1

    def update_total_length(self):
        """
        Update the total length of this subtree.

        This method recalculates the total length of the subtree rooted at this node
        by summing the lengths of the left subtree, the current node, and the right subtree.
        """
        left_length = self.left.total_length if self.left else 0
        right_length = self.right.total_length if self.right else 0
        self.total_length = left_length + self.length + right_length

    def update(self):
        self.update_height()
        self.update_total_length()

    def _get_balance_factor(self):
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

        self.update()
        new_root.update()
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

        self.update()
        new_root.update()
        return new_root

    def balance(self):
        """
        Balance this node if needed.

        Returns:
            SequenceNode: The new root after balancing
        """
        self.update_height()
        balance = self._get_balance_factor()

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
        return self.left._get_balance_factor() if self.left else 0

    def __get_right_balance(self):
        """Get balance factor of right child"""
        return self.right._get_balance_factor() if self.right else 0


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
        
        # 创建嵌套关系图
        self.nesting_graph = DonorNestingGraph()
        
        # 将根节点添加到嵌套关系图
        self.nesting_graph.add_node(self.root.uid, initial_seq)

    def __str__(self) -> str:
        """
        Convert the tree to a string.

        Returns:
            str: The concatenated sequence
        """
        return str(self.root) if self.root else ""

    def __iter__(self):
        yield from self.root

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

    def _create_node(self, data: str, is_donor: bool = False, donor_id: str = None) -> SequenceNode:
        """
        Create a new SequenceNode with a unique UID.

        Args:
            data (str): Sequence data
            is_donor (bool): Whether this node contains a donor sequence
            donor_id (str): Donor ID for tracking and visualization

        Returns:
            SequenceNode: The newly created node
        """
        uid = self._get_next_uid()
        node = SequenceNode(data, is_donor, donor_id, uid)
        self.node_dict[uid] = node
        return node

    def insert(self, abs_position: int, donor_seq: str, donor_id: str = None, tsd_length: int = 0, recursive: bool = False) -> None:
        """
        Insert a donor sequence at the specified position.

        Args:
            abs_position (int): Absolute position for insertion (1-based)
            donor_seq (str): Donor sequence to insert
            donor_id (str): Identifier for the donor sequence
            tsd_length (int): Length of Target Site Duplication (TSD) to generate
            recursive (bool): Whether to use recursive insertion method
        """
        # Skip insertion if donor sequence is empty
        if not donor_seq:
            return

        # Convert 1-based position to 0-based for internal processing
        zero_based_position = abs_position - 1

        if not recursive:
            self.root = self._insert_iterative(self.root, zero_based_position, donor_seq, donor_id, tsd_length)
        else:
            self.root = self._insert_recursive(self.root, zero_based_position, donor_seq, donor_id, tsd_length)

    def _insert_recursive(self, node: SequenceNode, abs_position: int, donor_seq: str, 
                      donor_id: str = None, tsd_length: int = 0) -> SequenceNode:
        """
        Recursively insert a donor sequence at the absolute position in the tree.

        Args:
            node (SequenceNode): Current node
            abs_position (int): Absolute position in the tree to insert at (0-based)
            donor_seq (str): Donor sequence to insert
            donor_id (str): Identifier for the donor sequence
            tsd_length (int): Length of Target Site Duplication (TSD) to generate

        Returns:
            SequenceNode: New node after insertion
        """
        # Skip insertion if donor sequence is empty
        if not donor_seq:
            return node

        # 创建donor节点并添加到嵌套关系图
        donor_node_uid = self._get_next_uid()  # 提前创建UID以在关系图中使用
        self.nesting_graph.add_node(donor_node_uid, donor_seq, donor_id, abs_position)

        # Calculate positions in tree
        node_start = node.left.total_length if node.left else 0
        node_end = node_start + node.length

        # Case 1: Position is in the current node's left subtree
        if abs_position <= node_start:
            if node.left:
                node.left = self._insert_recursive(node.left, abs_position, donor_seq, donor_id, tsd_length)
            else:
                # Insert as left child
                new_donor = self._create_node(donor_seq, True, donor_id)
                new_donor.uid = donor_node_uid  # 使用预先创建的UID
                self.node_dict[donor_node_uid] = new_donor
                node.left = new_donor
                
                # 如果当前节点是donor，记录嵌套关系
                if node.is_donor:
                    self.nesting_graph.add_nesting_relation(node.uid, donor_node_uid)
                    
            node.update()
            return node.balance()

        # Case 2: Position is inside the current node
        if node_start < abs_position < node_end:
            # Split this node
            rel_pos = abs_position - node_start
            left_data = node.data[:rel_pos]
            right_data = node.data[rel_pos:]

            # Handle TSD generation
            tsd_5 = tsd_3 = ""
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

            # 确定当前节点是否是donor
            is_current_donor = node.is_donor
            current_donor_id = node.donor_id

            # 创建左、右片段节点
            left_node_uid = self._get_next_uid()
            right_node_uid = self._get_next_uid()
            
            # 添加左右片段节点到关系图
            self.nesting_graph.add_node(left_node_uid, left_data, current_donor_id, node_start)
            self.nesting_graph.add_node(right_node_uid, right_data, current_donor_id, abs_position + len(donor_seq))

            # Create new left child with left data
            new_left = self._create_node(left_data, is_current_donor, current_donor_id)
            new_left.uid = left_node_uid
            self.node_dict[left_node_uid] = new_left
            
            if node.left:
                new_left.left = node.left
                new_left.update()
                # Balance left subtree
                new_left = new_left.balance()

            # Create new right child with right data and original right child
            new_right = self._create_node(right_data, is_current_donor, current_donor_id)
            new_right.uid = right_node_uid
            self.node_dict[right_node_uid] = new_right
            
            if node.right:
                new_right.right = node.right
                new_right.update()
                # Balance right subtree
                new_right = new_right.balance()

            # 如果当前节点是donor，记录切割关系
            if is_current_donor:
                # 添加切割关系到嵌套关系图
                self.nesting_graph.add_cut_relation(
                    donor_node_uid,  # 切割者donor的UID
                    node.uid,        # 被切割donor的UID
                    left_node_uid,   # 切割后左片段的UID
                    right_node_uid   # 切割后右片段的UID
                )
                
                # 将切割者donor添加为两个片段的容器（建立嵌套关系）
                self.nesting_graph.add_nesting_relation(donor_node_uid, left_node_uid)
                self.nesting_graph.add_nesting_relation(donor_node_uid, right_node_uid)
                
                # 如果当前节点已经在嵌套关系中，需要更新嵌套关系
                containers = self.nesting_graph.get_containers(node.uid)
                for container_uid in containers:
                    # 将嵌套关系从父节点转移到两个片段
                    self.nesting_graph.add_nesting_relation(container_uid, left_node_uid)
                    self.nesting_graph.add_nesting_relation(container_uid, right_node_uid)
            
            # Replace current node's data with the donor sequence
            node.data = donor_seq
            node.length = len(donor_seq)
            node.is_donor = True
            node.donor_id = donor_id
            node.uid = donor_node_uid
            self.node_dict[donor_node_uid] = node

            # Set new children
            node.left = new_left
            node.right = new_right
            node.update()

            return node.balance()

        # Case 3: Position is in the current node's right subtree
        if abs_position >= node_end:
            if node.right:
                node.right = self._insert_recursive(node.right, abs_position - node_end, donor_seq, donor_id, tsd_length)
            else:
                # Insert as right child
                new_donor = self._create_node(donor_seq, True, donor_id)
                new_donor.uid = donor_node_uid  # 使用预先创建的UID
                self.node_dict[donor_node_uid] = new_donor
                node.right = new_donor
                
                # 如果当前节点是donor，记录嵌套关系
                if node.is_donor:
                    self.nesting_graph.add_nesting_relation(node.uid, donor_node_uid)
                    
            node.update()
            return node.balance()

        raise RuntimeError("[ERROR] Should not reach here")

    def _insert_iterative(self, node: SequenceNode, abs_position: int, donor_seq: str, 
                          donor_id: str = None, tsd_length: int = 0) -> SequenceNode:
        """
        Iteratively insert a donor sequence at the absolute position in the tree.

        Args:
            node (SequenceNode): Current root node
            abs_position (int): Absolute position in the tree to insert at (0-based)
            donor_seq (str): Donor sequence to insert
            donor_id (str): Identifier for the donor sequence
            tsd_length (int): Length of Target Site Duplication (TSD) to generate

        Returns:
            SequenceNode: New root node after insertion
        """
        # Skip insertion if donor sequence is empty
        if not donor_seq:
            return node

        current = node
        parent_stack = []
        path_directions = []  # Record the path direction from root to current node ('left' or 'right')

        # 创建donor节点并添加到嵌套关系图
        donor_node_uid = self._get_next_uid()  # 提前创建UID以在关系图中使用
        self.nesting_graph.add_node(donor_node_uid, donor_seq, donor_id, abs_position)

        # Iteratively find insertion position
        while True:
            node_start = current.left.total_length if current.left else 0
            node_end = node_start + current.length

            # Case 1: Position is in the current node's left subtree
            if abs_position <= node_start:
                if current.left:
                    # Record parent node and direction for backtracking
                    parent_stack.append(current)
                    path_directions.append('left')
                    current = current.left
                else:
                    # Create new left child node
                    new_donor = self._create_node(donor_seq, True, donor_id)
                    new_donor.uid = donor_node_uid  # 使用预先创建的UID
                    self.node_dict[donor_node_uid] = new_donor
                    current.left = new_donor
                    
                    # 如果当前节点是donor，记录嵌套关系
                    if current.is_donor:
                        self.nesting_graph.add_nesting_relation(current.uid, donor_node_uid)
                        
                    current.update()
                    break

            # Case 2: Position is inside the current node
            elif node_start < abs_position < node_end:
                # Calculate relative position
                rel_pos = abs_position - node_start
                left_data = current.data[:rel_pos]
                right_data = current.data[rel_pos:]

                # Handle TSD generation
                tsd_5 = tsd_3 = ""
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

                # 确定当前节点是否是donor
                is_current_donor = current.is_donor
                current_donor_id = current.donor_id

                # 创建左、右片段节点
                left_node_uid = self._get_next_uid()
                right_node_uid = self._get_next_uid()
                
                # 添加左右片段节点到关系图
                self.nesting_graph.add_node(left_node_uid, left_data, current_donor_id, node_start)
                self.nesting_graph.add_node(right_node_uid, right_data, current_donor_id, abs_position + len(donor_seq))

                # Create new left child with left data
                new_left = self._create_node(left_data, is_current_donor, current_donor_id)
                new_left.uid = left_node_uid
                self.node_dict[left_node_uid] = new_left
                
                if current.left:
                    new_left.left = current.left
                    new_left.update()

                # Create new right child with right data and original right child
                new_right = self._create_node(right_data, is_current_donor, current_donor_id)
                new_right.uid = right_node_uid
                self.node_dict[right_node_uid] = new_right
                
                if current.right:
                    new_right.right = current.right
                    new_right.update()

                # 如果当前节点是donor，记录切割关系
                if is_current_donor:
                    # 添加切割关系到嵌套关系图
                    self.nesting_graph.add_cut_relation(
                        donor_node_uid,  # 切割者donor的UID
                        current.uid,     # 被切割donor的UID
                        left_node_uid,   # 切割后左片段的UID
                        right_node_uid   # 切割后右片段的UID
                    )
                    
                    # 将切割者donor添加为两个片段的容器（建立嵌套关系）
                    self.nesting_graph.add_nesting_relation(donor_node_uid, left_node_uid)
                    self.nesting_graph.add_nesting_relation(donor_node_uid, right_node_uid)
                    
                    # 如果当前节点已经在嵌套关系中，需要更新嵌套关系
                    containers = self.nesting_graph.get_containers(current.uid)
                    for container_uid in containers:
                        # 将嵌套关系从父节点转移到两个片段
                        self.nesting_graph.add_nesting_relation(container_uid, left_node_uid)
                        self.nesting_graph.add_nesting_relation(container_uid, right_node_uid)
                elif not is_current_donor:
                    # 如果不是donor，仍然要添加节点到图中，但不记录嵌套和切割关系
                    # 因为我们需要在图中保持完整的树结构，目前只有添加节点相关代码
                    pass

                # Replace current node's data with the donor sequence
                current.data = donor_seq
                current.length = len(donor_seq)
                current.is_donor = True
                current.donor_id = donor_id
                current.uid = donor_node_uid
                self.node_dict[donor_node_uid] = current

                # Set new children
                current.left = new_left
                current.right = new_right
                current.update()
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
                    new_donor = self._create_node(donor_seq, True, donor_id)
                    new_donor.uid = donor_node_uid  # 使用预先创建的UID
                    self.node_dict[donor_node_uid] = new_donor
                    current.right = new_donor
                    
                    # 如果当前节点是donor，记录嵌套关系
                    if current.is_donor:
                        self.nesting_graph.add_nesting_relation(current.uid, donor_node_uid)
                        
                    current.update()
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

            parent.update()

            # Balance parent node
            current = parent.balance()

        return current

    def donors(self, seq_id: str) -> Tuple[List[SeqRecord], List[SeqRecord]]:
        """
        收集所有donor节点并重建嵌套donor。
        利用嵌套关系图执行高效重建。

        Args:
            seq_id (str): 原始序列ID

        Returns:
            Tuple[List[SeqRecord], List[SeqRecord]]: 
                - 普通donor记录（不包括被重建donor覆盖的）
                - 重建的donor记录
        """
        # 收集所有donor节点
        donor_records = self._collect_donor_records(seq_id)
        
        # 使用嵌套关系图重建donor
        reconstructed_donors, excluded_uids = self.nesting_graph.reconstruct_donors(seq_id)
        
        # 过滤掉被重建donor覆盖的记录
        if excluded_uids:
            donor_records = [record for record in donor_records 
                             if record.annotations.get("uid") not in excluded_uids]

        return donor_records, reconstructed_donors

    def _collect_donor_records(self, seq_id: str) -> List[SeqRecord]:
        """
        从树中收集所有donor节点，生成SeqRecord记录
        
        Args:
            seq_id (str): 原始序列ID
            
        Returns:
            List[SeqRecord]: Donor记录列表
        """
        donor_records = []
        abs_position_map = self._calculate_absolute_positions()
        
        # 遍历所有节点
        for node in self:
            if node.is_donor:
                # 获取节点的绝对位置信息
                start_pos = abs_position_map.get(node.uid, 0)
                end_pos = start_pos + node.length
                
                # 转换为1-based索引用于输出
                start_pos_1based = start_pos + 1
                end_pos_1based = end_pos
                
                # 创建donor ID
                donor_id = f"{seq_id}_{start_pos_1based}_{end_pos_1based}-+-{node.length}"
                if node.donor_id:
                    donor_id += f"-{node.donor_id}"
                    
                # 创建record
                donor_record = create_sequence_record(node.data, donor_id)
                donor_record.annotations["uid"] = node.uid
                donor_record.annotations["position"] = start_pos_1based
                donor_record.annotations["length"] = node.length
                
                # 添加到结果列表
                donor_records.append(donor_record)
                
        return donor_records

    def _calculate_absolute_positions(self) -> Dict[int, int]:
        """
        计算树中每个节点的绝对位置
        
        Returns:
            Dict[int, int]: 节点UID到绝对位置的映射
        """
        positions = {}
        
        def _traverse(node, current_pos=0):
            if not node:
                return current_pos
                
            # 处理左子树
            left_end_pos = _traverse(node.left, current_pos)
            
            # 计算当前节点位置
            node_pos = left_end_pos
            positions[node.uid] = node_pos
            
            # 处理右子树
            right_end_pos = _traverse(node.right, node_pos + node.length)
            
            return right_end_pos
            
        _traverse(self.root, 0)
        return positions

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
                   "  bgcolor=\"#FFFFFF\"",
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
            node_id_prefix (str): Prefix for node IDs (not used when using uid as node ID)
            abs_pos (int): Current absolute position (0-based)

        Returns:
            tuple: (nodes list, edges list)
        """
        if not node:
            return [], []

        nodes = []
        edges = []

        # Calculate positions
        left_length = node.left.total_length if node.left else 0

        # Calculate start and end positions for display
        start_pos = abs_pos + left_length
        end_pos = start_pos + node.length

        # Convert to 1-based positions for display
        start_pos_1based = start_pos + 1
        end_pos_1based = end_pos

        # Use uid directly as node ID
        node_id = f"node_{node.uid}"

        # Determine node type and color
        node_type = "Donor" if node.is_donor else "Acceptor"
        fill_color = "lightblue" if node.is_donor else "lightgreen"
        
        # Process fragment and nesting information
        nested_in = ""
        cut_half = ""
        
        # Check if this is a fragment of a cut donor
        if self.nesting_graph.is_fragment(node.uid):
            # Get fragment info (original_uid, position, is_left)
            fragment_info = self.nesting_graph.fragments.get(node.uid)
            if fragment_info:
                orig_uid, _, is_left = fragment_info
                half_type = "L" if is_left else "R"
                cut_half = f"Cut: {half_type}\\n"
                fill_color = "lightpink"  # Cut fragments shown in pink
        
        # Check for nesting relationships
        containers = self.nesting_graph.get_containers(node.uid)
        if containers:
            nested_in = "Nest: " + ','.join(map(str, containers)) + "\\n"
            fill_color = "yellow"  # Nested nodes shown in yellow
        
        # If both nested and cut, use a distinctive color
        if nested_in and cut_half:
            fill_color = "plum"

        # Create node label with position information
        label = "".join([node_type, " | ", str(node.uid), "\\n",
                        str(start_pos_1based), "\\l",
                        str(end_pos_1based), "\\l",
                        "Length: ", str(node.length), "\\n",
                        nested_in,
                        cut_half])

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

            # Add edge from this node to left child using uid
            left_id = f"node_{node.left.uid}"
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

            # Add edge from this node to right child using uid
            right_id = f"node_{node.right.uid}"
            edges.append(f'{node_id} -> {right_id} [label="R"];')

        return nodes, edges

class DonorNestingGraph:
    """
    表示donor序列嵌套关系的图数据结构。
    用于高效追踪和重建嵌套的donor序列，替代之前在SequenceNode中的属性方式。
    """
    def __init__(self):
        # 节点信息：UID -> {sequence, length, donor_id, positions}
        self.nodes = {}
        
        # 嵌套关系：container_uid -> [nested_uid1, nested_uid2, ...]
        # 记录哪些donor嵌套在指定donor内
        self.nestings = {}
        
        # 反向嵌套关系：nested_uid -> [container_uid1, container_uid2, ...]
        # 记录指定donor嵌套在哪些donor内
        self.nested_in = {}
        
        # 切割关系：cutter_uid -> [(cut_uid, left_uid, right_uid), ...]
        # 记录哪个donor切割了哪个donor，产生了哪两个新片段
        self.cuts = {}
        
        # 反向切割关系：cut_uid -> [(cutter_uid, left_uid, right_uid), ...]
        # 记录哪个donor被哪些donor切割，产生了哪些片段对
        self.cut_by = {}
        
        # 片段关系：fragment_uid -> (original_uid, position, is_left)
        # 记录每个片段来自哪个原始donor，位于哪个位置，是左半部分还是右半部分
        self.fragments = {}
        
        # 返回父片段映射：left_uid/right_uid -> original_uid
        # 记录左/右片段对应的原始donor
        self.fragment_to_original = {}

    def add_node(self, uid: int, sequence: str, donor_id: str = None, position: int = None):
        """
        添加一个donor节点到图中
        
        Args:
            uid: 节点唯一标识符
            sequence: 节点包含的序列
            donor_id: donor序列标识符
            position: 插入位置
        """
        self.nodes[uid] = {
            'sequence': sequence,
            'length': len(sequence),
            'donor_id': donor_id,
            'position': position
        }
        
        # 初始化关系映射
        if uid not in self.nestings:
            self.nestings[uid] = []
        if uid not in self.nested_in:
            self.nested_in[uid] = []
            
        return self

    def add_nesting_relation(self, container_uid: int, nested_uid: int):
        """
        添加嵌套关系
        
        Args:
            container_uid: 容器donor的UID
            nested_uid: 嵌套在其中的donor的UID
        """
        # 确保节点存在
        if container_uid not in self.nodes or nested_uid not in self.nodes:
            return self
        
        # 添加正向和反向映射
        if container_uid not in self.nestings:
            self.nestings[container_uid] = []
        if nested_uid not in self.nestings[container_uid]:
            self.nestings[container_uid].append(nested_uid)
            
        if nested_uid not in self.nested_in:
            self.nested_in[nested_uid] = []
        if container_uid not in self.nested_in[nested_uid]:
            self.nested_in[nested_uid].append(container_uid)
            
        return self
        
    def add_cut_relation(self, cutter_uid: int, cut_uid: int, left_uid: int, right_uid: int):
        """
        添加切割关系
        
        Args:
            cutter_uid: 切割者donor的UID
            cut_uid: 被切割donor的UID
            left_uid: 切割后左片段的UID
            right_uid: 切割后右片段的UID
        """
        # 确保节点存在
        if not all(uid in self.nodes for uid in [cutter_uid, cut_uid, left_uid, right_uid]):
            return self
            
        # 添加切割关系
        if cutter_uid not in self.cuts:
            self.cuts[cutter_uid] = []
        self.cuts[cutter_uid].append((cut_uid, left_uid, right_uid))
        
        # 添加反向切割关系
        if cut_uid not in self.cut_by:
            self.cut_by[cut_uid] = []
        self.cut_by[cut_uid].append((cutter_uid, left_uid, right_uid))
        
        # 记录片段关系
        cut_position = self.nodes[cutter_uid]['position']
        self.fragments[left_uid] = (cut_uid, cut_position, True)  # True表示左片段
        self.fragments[right_uid] = (cut_uid, cut_position, False)  # False表示右片段
        
        # 添加返回父片段的映射
        self.fragment_to_original[left_uid] = cut_uid
        self.fragment_to_original[right_uid] = cut_uid
        
        return self
    
    def is_fragment(self, uid: int) -> bool:
        """判断指定UID的节点是否为片段"""
        return uid in self.fragments
    
    def get_original_donor(self, uid: int) -> int:
        """获取片段对应的原始donor UID"""
        return self.fragment_to_original.get(uid, uid)
    
    def get_all_fragments(self, original_uid: int) -> list:
        """获取原始donor的所有片段"""
        return [uid for uid, (orig, _, _) in self.fragments.items() if orig == original_uid]
    
    def get_nested_donors(self, container_uid: int) -> list:
        """获取嵌套在指定donor中的所有donor"""
        return self.nestings.get(container_uid, [])
    
    def get_containers(self, nested_uid: int) -> list:
        """获取包含指定donor的所有donor"""
        return self.nested_in.get(nested_uid, [])
    
    def get_cut_by(self, cut_uid: int) -> list:
        """获取切割了指定donor的所有donor信息"""
        return self.cut_by.get(cut_uid, [])
    
    def get_cuts(self, cutter_uid: int) -> list:
        """获取被指定donor切割的所有donor信息"""
        return self.cuts.get(cutter_uid, [])
    
    def reconstruct_donors(self, seq_id: str) -> tuple:
        """
        基于嵌套和切割关系重建donor序列
        
        Args:
            seq_id: 序列标识符，用于构建输出记录ID
            
        Returns:
            tuple: (重建的donor记录列表, 排除的UID集合)
        """
        MAX_RECURSION_DEPTH = 50
        reconstructed = []
        excluded = set()
        processed_fragments = set()
        
        # 首先查找所有需要递归重建的片段
        def process_fragment_tree(fragment_uid, depth=0):
            if fragment_uid in processed_fragments or depth > MAX_RECURSION_DEPTH:
                return
            processed_fragments.add(fragment_uid)
            
            # 检查此片段是否被进一步切割
            for cutter_uid, left_uid, right_uid in self.cut_by.get(fragment_uid, []):
                # 递归处理子片段
                process_fragment_tree(left_uid, depth+1)
                process_fragment_tree(right_uid, depth+1)
        
        # 处理所有被切割的donor，从底层往上重建
        for cut_uid in list(self.cut_by.keys()):
            process_fragment_tree(cut_uid)
        
        # 按照处理顺序排序，确保先处理深层片段
        fragments_to_process = sorted(list(processed_fragments), 
                                key=lambda uid: len(self.cut_by.get(uid, [])), reverse=True)
        
        # 处理所有片段
        for fragment_uid in fragments_to_process:
            if fragment_uid in excluded:
                continue
                
            for cutter_uid, left_uid, right_uid in self.cut_by.get(fragment_uid, []):
                # 跳过已处理的
                if any(uid in excluded for uid in [fragment_uid, cutter_uid, left_uid, right_uid]):
                    continue
                    
                # 检查我们是否有所有需要的节点数据
                if not all(uid in self.nodes for uid in [fragment_uid, cutter_uid, left_uid, right_uid]):
                    continue
                
                # 获取序列数据
                try:
                    left_seq = self.nodes[left_uid]['sequence']
                    right_seq = self.nodes[right_uid]['sequence']
                    cutter_seq = self.nodes[cutter_uid]['sequence']
                except KeyError:
                    continue
                
                # 重建两种形式的donor：带切割donor和不带切割donor
                
                # 1. 包含切割donor的完整重建
                full_seq = left_seq + cutter_seq + right_seq
                full_id = f"{seq_id}_reconstructed_{fragment_uid}_{cutter_uid}"
                full_rec = create_sequence_record(full_seq, full_id)
                full_rec.annotations["reconstruction_type"] = "full"
                full_rec.annotations["original_uid"] = fragment_uid
                full_rec.annotations["cutter_uid"] = cutter_uid
                full_rec.annotations["left_uid"] = left_uid
                full_rec.annotations["right_uid"] = right_uid
                full_rec.annotations["reconstruction_depth"] = fragments_to_process.index(fragment_uid)
                reconstructed.append(full_rec)
                
                # 2. 不包含切割donor的清洁重建
                clean_seq = left_seq + right_seq
                clean_id = f"{seq_id}_clean_reconstructed_{fragment_uid}"
                clean_rec = create_sequence_record(clean_seq, clean_id)
                clean_rec.annotations["reconstruction_type"] = "clean"
                clean_rec.annotations["original_uid"] = fragment_uid
                clean_rec.annotations["left_uid"] = left_uid
                clean_rec.annotations["right_uid"] = right_uid
                clean_rec.annotations["reconstruction_depth"] = fragments_to_process.index(fragment_uid)
                reconstructed.append(clean_rec)
                
                # 添加到排除集
                excluded.update([fragment_uid, cutter_uid, left_uid, right_uid])
        
        return reconstructed, excluded
    
    def to_graphviz_dot(self, prefix: str = "donor_graph") -> str:
        """
        生成Graphviz DOT格式可视化
        
        Args:
            prefix: 图形文件名前缀
            
        Returns:
            str: Graphviz DOT格式字符串
        """
        lines = ["digraph DonorNestingGraph {", 
                 "  bgcolor=\"#FFFFFF\";", 
                 "  node [shape=box, style=filled, fontsize=10];"]
                 
        # 获取有关系的节点UID列表
        related_uids = set()
        
        # 添加所有参与嵌套关系的节点
        for container_uid, nested_list in self.nestings.items():
            related_uids.add(container_uid)
            related_uids.update(nested_list)
            
        # 添加所有参与切割关系的节点
        for cutter_uid, cut_list in self.cuts.items():
            related_uids.add(cutter_uid)
            for cut_uid, left_uid, right_uid in cut_list:
                related_uids.update([cut_uid, left_uid, right_uid])
                
        # 添加所有片段节点
        related_uids.update(self.fragments.keys())
        
        # 添加所有在cut_by中的节点
        for cut_uid in self.cut_by:
            related_uids.add(cut_uid)
            for cutter_uid, left_uid, right_uid in self.cut_by[cut_uid]:
                related_uids.update([cutter_uid, left_uid, right_uid])
                
        # 添加所有在nested_in中的节点
        for nested_uid in self.nested_in:
            related_uids.add(nested_uid)
            related_uids.update(self.nested_in[nested_uid])
                 
        # 添加所有有关系的节点
        for uid in related_uids:
            if uid not in self.nodes:
                continue
                
            info = self.nodes[uid]
            node_type = "Donor"
            fragment_info = ""
            if self.is_fragment(uid):
                orig_uid, pos, is_left = self.fragments[uid]
                side = "Left" if is_left else "Right"
                fragment_info = f"\\nFragment of {orig_uid} ({side})"
                node_type = "Fragment"
                
            # 安全获取节点长度
            length = info.get('length', 0)
            donor_id = info.get('donor_id', '')
            
            # 创建标签，包含donor_id信息
            label = f"{node_type}\\nUID: {uid}\\nLen: {length}"
            if donor_id:
                label += f"\\nID: {donor_id}"
            label += fragment_info
            
            # 设置颜色
            color = "#AAFFAA"  # 默认绿色
            if fragment_info:
                color = "#FFAAAA"  # 片段为红色
            elif node_type == "Donor":
                color = "#AAAAFF"  # donor为蓝色
                
            lines.append(f'  node_{uid} [label="{label}", fillcolor="{color}"];')
            
        # 添加嵌套关系边
        for container_uid, nested_list in self.nestings.items():
            for nested_uid in nested_list:
                lines.append(f'  node_{container_uid} -> node_{nested_uid} [label="contains", color="green"];')
                
        # 添加切割关系边
        for cutter_uid, cut_list in self.cuts.items():
            for cut_uid, left_uid, right_uid in cut_list:
                lines.append(f'  node_{cutter_uid} -> node_{cut_uid} [label="cuts", color="red"];')
                lines.append(f'  node_{cut_uid} -> node_{left_uid} [label="left", style="dashed", color="blue"];')
                lines.append(f'  node_{cut_uid} -> node_{right_uid} [label="right", style="dashed", color="blue"];')
                
        lines.append("}")
        return '\n'.join(lines)

# ======================================================================================================================
# Utility Functions

def convert_humanized_number(text: str, base: Union[float, int], units: Sequence[Union[str, Iterable[str]]]) -> float:
    """
    Convert human-readable numbers with unit suffixes to numerical values.

    Parses strings containing numbers with optional unit suffixes and converts them
    to raw numerical values based on the provided base and unit definitions.

    Args:
        text: Input string containing the number and optional unit suffix
        base: Numerical base for unit conversion (e.g., 1000 or 1024)
        units: Sequence of unit suffixes ordered by increasing scale

    Returns:
        float: The converted numerical value
    
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
    """
    Convert a human-readable number string (e.g., '1k', '2M') to an integer.
    
    Args:
        number: Either an integer or a string with optional k/m/b/t suffix
        
    Returns:
        int: The converted integer value
    """
    if isinstance(number, int):
        return number
    return int(convert_humanized_number(number, 1000, ("", 'K', 'M', 'B', 'T')))

def sort_multiple_lists(base: list, *lists: list,
                        key: Optional[Callable] = None, reverse: bool = False) -> Union[list, tuple]:
    """
    Synchronously sorts multiple lists based on the sorting order of the base list.

    Args:
        base: The base list used as the reference for sorting
        *lists: Other lists to be sorted synchronously
        key: A function to execute to decide the order
        reverse: Whether to sort in descending order

    Returns:
        tuple: A tuple containing all sorted lists
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

    Args:
        seq: The biological sequence as a string
        id: The unique identifier for the sequence
        
    Returns:
        SeqRecord: A SeqRecord object containing the sequence and its metadata
    """
    return SeqRecord(Seq(seq), id=id, description="")

def save_multi_fasta_from_dict(records_dict: Dict[str, List[SeqRecord]], output_dir_path: str):
    """
    Save multiple FASTA files from a dictionary of sequence records.

    Args:
        records_dict: Dictionary mapping file names to lists of SeqRecord objects
        output_dir_path: Directory path where the output files will be saved
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

    Args:
        path: The path to search for donor library

    Returns:
        Optional[List[str]]: A list of absolute paths to the found donor library files,
                           or None if path is None
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
            if file.endswith((".fa", ".fasta")):
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
                if file.endswith((".fa", ".fasta")):
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
        path_list: List of paths to donor libraries
        weight_list: List of weights for each donor library
        len_limit: Maximum length limit for sequences to load
        flag_filter_n: Flag to filter out sequences containing N

    Returns:
        Tuple containing:
            - Combined list of donor sequences from all libraries
            - List of lengths of donor sequences
            - List of weights for each donor sequence
    """
    if not path_list:
        return None

    all_donor_sequence_list: List[str] = []
    all_donor_len_list: List[int] = []
    all_donor_weight_list: List[float] = []

    for lib_path, weight in zip(path_list, weight_list or [None] * len(path_list)):
        single_donor_lib_abs_path_list = _find_donor_lib_abs_path_list(lib_path)
        if single_donor_lib_abs_path_list:
            single_donor_lib_sequences = load_sequences(single_donor_lib_abs_path_list, len_limit, flag_filter_n)
            if single_donor_lib_sequences:
                all_donor_sequence_list.extend(single_donor_lib_sequences)
                all_donor_len_list.extend([len(seq) for seq in single_donor_lib_sequences])
                if weight:
                    len_single_donor_lib_sequences = len(single_donor_lib_sequences)
                    all_donor_weight_list.extend([weight / len_single_donor_lib_sequences] * len_single_donor_lib_sequences)

    # Sort by length
    all_donor_len_list, all_donor_sequence_list, all_donor_weight_list = sort_multiple_lists(
        all_donor_len_list, all_donor_sequence_list, all_donor_weight_list)

    # Uniform weights by default: If no weight provided, set all weights to 1
    if not weight_list:
        all_donor_weight_list = [1] * len(all_donor_sequence_list)

    return all_donor_sequence_list, all_donor_len_list, all_donor_weight_list

# ======================================================================================================================
# Main Sequence Generator Class

class SeqGenerator:
    def __init__(self, input_file: str, insertion: Union[str, int], batch: int, processors: int, output_dir_path: str,
                 donor_lib: Optional[List[str]] = None, donor_lib_weight: Optional[List[float]] = None,
                 donor_len_limit: Optional[int] = None, flag_filter_n: bool = False, flag_track: bool = False,
                 tsd_length: Optional[int] = None, flag_visual: bool = False, flag_recursive: bool = False, 
                 iteration: int = 1):
        """
        Initialize the sequence generator.

        Args:
            input_file: Path to the input sequence file
            insertion: Number of insertions per sequence
            batch: Number of independent result files to generate
            processors: Number of processors to use
            output_dir_path: Directory to save output files
            donor_lib: List of donor library file paths
            donor_lib_weight: Weights for donor libraries
            donor_len_limit: Maximum length limit for donor sequences to load
            flag_filter_n: Whether to filter out sequences containing N
            flag_track: Whether to track donor sequences used
            tsd_length: Length of Target Site Duplication (TSD) to generate
            flag_visual: Whether to generate graphviz visualization of the sequence tree
            iteration: Number of insertion iterations to perform on each sequence
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
        self.iteration: int = iteration

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

    def __print_header(self):
        """
        Print the program header with basic information.
        """
        print(f"=== RandSeqInsert {VERSION} ===")
        print(f"Processing input file: {self.input_file}")
        print(f"Number of input sequences: {len(self.input)}")
        print(f"Insertion settings: {self.insertion} insertions per sequence")
        if self.iteration > 1:
            print(f"Performing {self.iteration} iterations of insertion")
        if self.tsd_length:
            print(f"TSD settings: Generating TSD of length {self.tsd_length} at insertion sites")
        print(f"Donor library: {len(self.donor_sequences)} sequences loaded")
        print(f"Generating {self.batch} independent result file(s)")
        if self.flag_recursive:
            print(f"Using recursive insertion method")
        if self.flag_visual:
            print(f"Tree and Graph visualization enabled")

    def __pre_check(self):
        """
        Perform pre-execution checks.
        """
        if not os.path.exists(self.output_dir_path):
            os.makedirs(self.output_dir_path)

    def __process_batch_multiprocessing(self) -> Tuple[List[SeqRecord], 
                                                      Optional[List[SeqRecord]],
                                                      Optional[List[SeqRecord]]]:
        """
        Process a batch of sequences using multiprocessing.

        Returns:
            Tuple containing:
                - All processed sequences
                - Donor records (if tracking enabled)
                - Reconstructed donor records (if tracking enabled)
        """
        if not self.input:
            return [], None, None

        total_sequences = len(self.input)
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
        Process a single sequence with random insertions, potentially over multiple iterations.

        Args:
            idx_or_record: Either an index to look up in self.input, or a SeqRecord to process directly

        Returns:
            Tuple containing:
                - The processed sequence
                - Used donors (if tracking enabled)
                - Reconstructed donors (if tracking enabled)
        """
        random.seed(os.getpid() + int(time.time()))
        # If idx_or_record is an integer, retrieve the sequence record
        # If it's already a SeqRecord, use it directly
        if isinstance(idx_or_record, int):
            seq_record = self.input[idx_or_record]
        else:
            seq_record = idx_or_record

        # Check if there are donor sequences to insert
        if not self.insertion or not self.donor_sequences:
            # If no insertions requested or no donor sequences available,
            # return the original sequence with modified ID
            new_id = f"{seq_record.id}_ins0"
            return create_sequence_record(str(seq_record.seq), new_id), None, None

        # Create a new sequence tree for this sequence
        seq_tree = SequenceTree(str(seq_record.seq), 0)
        # Track all used donors across iterations if tracking is enabled
        all_used_donors = [] if self.flag_track else None
        all_reconstructed_donors = [] if self.flag_track else None
        
        # Perform multiple iterations of insertion
        for iteration in range(1, self.iteration + 1):
            # Get the current total sequence length (updated after each iteration)
            current_seq = str(seq_tree)
            total_length = len(current_seq)

            # Generate insertion positions and donor sequences
            # Use 1-based positions (1 to total_length+1)
            insert_positions = random.choices(range(1, total_length + 2), k=self.insertion)
            selected_donors = random.choices(self.donor_sequences, weights=self.donor_weights, k=self.insertion)

            # Sort positions in descending order to maintain position integrity during insertion
            insert_positions, selected_donors = sort_multiple_lists(insert_positions, selected_donors, reverse=True)

            # Insert donors into the tree
            for pos, donor_seq in zip(insert_positions, selected_donors):
                seq_tree.insert(pos, donor_seq, None, self.tsd_length, self.flag_recursive)

        # Collect donor sequences once after all iterations if tracking is enabled
        if self.flag_track:
            used_donors, reconstructed_donors = seq_tree.donors(seq_record.id)
            all_used_donors = used_donors if used_donors else []
            all_reconstructed_donors = reconstructed_donors if reconstructed_donors else []

        # Generate final sequence
        new_id = f"{seq_record.id}_ins{self.insertion}"
        if self.iteration > 1:
            new_id += f"_iter{self.iteration}"
        if self.tsd_length:
            new_id += f"_tsd{self.tsd_length}"
        new_seq_record = create_sequence_record(str(seq_tree), new_id)

        # Generate visualization if requested
        if self.flag_visual:
            # Generate tree visualization
            graphviz_str = seq_tree.to_graphviz_dot(node_id_prefix=seq_record.id)
            tree_visual_dir_path = os.path.join(self.output_dir_path, VISUAL_DIR_NAME)
            os.makedirs(tree_visual_dir_path, exist_ok=True)
            with open(os.path.join(tree_visual_dir_path, f"{seq_record.id}_tree_visual.dot"), "w") as f:
                f.write(graphviz_str)
            
            # Generate nesting graph visualization
            graph_visual_dir_path = os.path.join(self.output_dir_path, VISUAL_DIR_NAME)
            os.makedirs(graph_visual_dir_path, exist_ok=True)
            nesting_graphviz_str = seq_tree.nesting_graph.to_graphviz_dot(prefix=seq_record.id)
            with open(os.path.join(graph_visual_dir_path, f"{seq_record.id}_graph_visual.dot"), "w") as f:
                f.write(nesting_graphviz_str)
            
            print(f"Generated tree and graph visualizations for {seq_record.id}")

        return new_seq_record, all_used_donors, all_reconstructed_donors

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

    def __print_summary(self, total_elapsed_time: float):
        """
        Print a summary of the processing results.

        Args:
            total_elapsed_time: Total time taken for all batches
        """
        print(f"\nAll batches completed in {total_elapsed_time:.2g} seconds")
        print(f"Results saved to \"{os.path.abspath(self.output_dir_path)}\"")

    def __process_single_batch(self, batch_num: int) -> float:
        """
        Process a single batch of sequences.

        Args:
            batch_num: Batch number (starting from 1)

        Returns:
            float: Time taken to process this batch (seconds)
        """
        batch_start_time = time.time()

        # Set random seed for this batch
        random.seed(int(time.time()) + batch_num)

        print(f"\nProcessing batch {batch_num}/{self.batch}")

        # Choose processing method
        all_sequences, all_donors, all_reconstructed_donors = self.__process_batch_multiprocessing()

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
            output_dir: Directory to save the results
            sequences: List of sequence records to save
            donors: List of donor records to save if tracking is enabled
            reconstructed_donors: List of reconstructed donor records
            batch_num: The batch number (1-based index)
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
                full_recs = [d for d in reconstructed_donors if "reconstruction_type" in d.annotations and d.annotations["reconstruction_type"] == "full"]
                clean_recs = [d for d in reconstructed_donors if "reconstruction_type" in d.annotations and d.annotations["reconstruction_type"] == "clean"]

                if full_recs:
                    print(f"Saving {len(full_recs)} full reconstructed donor records")
                    output_dict[f"reconstructed_donors{suffix}.fa"] = full_recs

                if clean_recs:
                    print(f"Saving {len(clean_recs)} clean reconstructed donor records")
                    output_dict[f"clean_reconstructed_donors{suffix}.fa"] = clean_recs

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
    core_group.add_argument("-it", "--iteration", type=int, default=1, metavar="INT",
                            help="Number of insertion iterations to perform on each sequence. Each iteration will use the sequence from the previous iteration as input. Default: 1")

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
        flag_recursive=parsed_args.recursive,
        iteration=parsed_args.iteration
    )
    generator.execute()

if __name__ == "__main__":
    main()
