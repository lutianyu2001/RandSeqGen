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
        实现对树的中序遍历，以左-根-右的顺序迭代所有节点
        Yields:
            SequenceNode: 树中的每个节点
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

    def insert(self, abs_position: int, donor_seq: str, donor_id: str = None, tsd_length: int = 0) -> None:
        """
        Insert a donor sequence at the specified position.

        Args:
            abs_position (int): Absolute position for insertion (1-based)
            donor_seq (str): Donor sequence to insert
            donor_id (str): Identifier for the donor sequence
            tsd_length (int): Length of Target Site Duplication (TSD) to generate
        """
        # Skip insertion if donor sequence is empty
        if not donor_seq:
            return

        # Convert 1-based position to 0-based for internal processing
        zero_based_position = abs_position - 1
        self.root = self._insert_iterative(self.root, zero_based_position, donor_seq, donor_id, tsd_length)

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

        return current if not parent_stack else node.balance()

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
            node_id_prefix (str): Prefix for node IDs 
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
        
        # 检查是否被切割
        is_fragment = self.nesting_graph.is_fragment(node.uid)
        if is_fragment:
            orig_uid = self.nesting_graph.get_original_donor(node.uid)
            frag_info = f"Fragment of {orig_uid}"
            fill_color = "lightpink"  # 被切割的donor显示为粉色
        else:
            frag_info = ""
            
        # 检查是否有嵌套关系
        nested_donors = self.nesting_graph.get_nested_donors(node.uid)
        nested_in = self.nesting_graph.get_containers(node.uid)
        
        if nested_donors and nested_in:
            # 既被嵌套又有嵌套的节点
            fill_color = "plum"
        elif nested_donors:
            # 含有嵌套donor的节点
            fill_color = "yellow"
        elif nested_in:
            # 被嵌套在其他donor中的节点
            fill_color = "orange"

        # Create node label with position information
        label = "".join([node_type, " | ", str(node.uid), "\\n",
                         str(start_pos_1based), "-", str(end_pos_1based), "\\n",
                         "Length: ", str(node.length), "\\n"])
                         
        if frag_info:
            label += frag_info + "\\n"
            
        if nested_donors:
            label += f"Contains: {len(nested_donors)} donors\\n"
            
        if nested_in:
            label += f"Nested in: {len(nested_in)} donors\\n"

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
        reconstructed = []
        excluded = set()
        
        # 处理所有被切割的donor
        for cut_uid, cut_info_list in self.cut_by.items():
            for cutter_uid, left_uid, right_uid in cut_info_list:
                # 确保节点存在
                if not all(uid in self.nodes for uid in [cut_uid, cutter_uid, left_uid, right_uid]):
                    continue
                
                # 如果已经处理过，跳过
                if cut_uid in excluded or cutter_uid in excluded:
                    continue
                
                # 获取序列
                left_seq = self.nodes[left_uid]['sequence']
                right_seq = self.nodes[right_uid]['sequence']
                cutter_seq = self.nodes[cutter_uid]['sequence']
                
                # 重建两种形式的donor：带切割donor和不带切割donor
                
                # 1. 包含切割donor的完整重建
                full_seq = left_seq + cutter_seq + right_seq
                full_id = f"{seq_id}_reconstructed_{cut_uid}_{cutter_uid}"
                full_rec = create_sequence_record(full_seq, full_id)
                full_rec.annotations["reconstruction_type"] = "full"
                full_rec.annotations["original_uid"] = cut_uid
                full_rec.annotations["cutter_uid"] = cutter_uid
                full_rec.annotations["left_uid"] = left_uid
                full_rec.annotations["right_uid"] = right_uid
                reconstructed.append(full_rec)
                
                # 2. 不包含切割donor的清洁重建
                clean_seq = left_seq + right_seq
                clean_id = f"{seq_id}_clean_reconstructed_{cut_uid}"
                clean_rec = create_sequence_record(clean_seq, clean_id)
                clean_rec.annotations["reconstruction_type"] = "clean"
                clean_rec.annotations["original_uid"] = cut_uid
                clean_rec.annotations["left_uid"] = left_uid
                clean_rec.annotations["right_uid"] = right_uid
                reconstructed.append(clean_rec)
                
                # 添加到排除集
                excluded.update([cut_uid, cutter_uid, left_uid, right_uid])
                
        # 递归处理可能的嵌套关系（处理嵌套在切割片段中的donor）
        # 这部分逻辑可以根据需要进一步扩展
                
        return reconstructed, excluded
    
    def to_graphviz(self, prefix: str = "donor_graph") -> str:
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
                 
        # 添加所有节点
        for uid, info in self.nodes.items():
            node_type = "Donor"
            fragment_info = ""
            if self.is_fragment(uid):
                orig_uid, pos, is_left = self.fragments[uid]
                side = "Left" if is_left else "Right"
                fragment_info = f"\\nFragment of {orig_uid} ({side})"
                node_type = "Fragment"
                
            label = f"{node_type}\\nUID: {uid}\\nLen: {info['length']}{fragment_info}"
            
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

class SeqGenerator:
    def __init__(self, input_file: str, insertion: Union[str, int], batch: int, processors: int, output_dir_path: str,
                 donor_lib: Optional[List[str]] = None, donor_lib_weight: Optional[List[float]] = None,
                 donor_len_limit: Optional[int] = None, flag_filter_n: bool = False, flag_track: bool = False,
                 tsd_length: Optional[int] = None, flag_visual: bool = False, iteration: int = 1):
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
            iteration (int): Number of insertion iterations to perform on each sequence
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
        self.iteration: int = iteration

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
        if self.flag_visual:
            print(f"Tree visualization enabled: DOT files will be generated for each sequence")

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
        iteration=parsed_args.iteration
    )
    generator.execute()

if __name__ == "__main__":
    main()
