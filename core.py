#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Tianyu (Sky) Lu (tianyu@lu.fm)

import random
from typing import Tuple, List, Dict, Optional
from Bio.SeqRecord import SeqRecord

from utils import create_sequence_record, generate_TSD


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

    def get_balance_factor(self):
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
        balance = self.get_balance_factor()

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
        return self.left.get_balance_factor() if self.left else 0

    def __get_right_balance(self):
        """Get balance factor of right child"""
        return self.right.get_balance_factor() if self.right else 0


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

        if uid not in self.available_uids:
            self.available_uids.append(uid)

    def _create_node(self, data: str, is_donor: bool = False, donor_id: str = None, uid: int = None) -> SequenceNode:
        """
        Create a new SequenceNode with a unique UID.

        Args:
            data (str): Sequence data
            is_donor (bool): Whether this node contains a donor sequence
            donor_id (str): Donor ID for tracking and visualization
            uid (int): 预设的UID，如果提供则使用此UID

        Returns:
            SequenceNode: The newly created node
        """
        uid = uid or self._get_next_uid()
        node = SequenceNode(data, is_donor, donor_id, uid)
        self.node_dict[uid] = node
        return node

    def insert(self, abs_position: int, donor_seq: str, donor_id: str = None, tsd_length: int = 0,
               recursive: bool = False, debug: bool = False) -> None:
        """
        Insert a donor sequence at the specified position.

        Args:
            debug:
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
            self.root = self._insert_iterative(self.root, zero_based_position, donor_seq, donor_id, tsd_length, debug)
        else:
            self.root = self._insert_recursive(self.root, zero_based_position, donor_seq, donor_id, tsd_length, debug)

    def _insert_recursive(self, node: SequenceNode, abs_position: int, donor_seq: str,
                          donor_id: str = None, tsd_length: int = 0, debug: bool = False) -> SequenceNode:
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

        # 创建donor节点UID
        donor_node_uid = self._get_next_uid()
        # 添加到嵌套关系图
        self.nesting_graph.add_node(donor_node_uid, donor_seq, donor_id, abs_position)

        # Calculate positions in tree
        node_start = node.left.total_length if node.left else 0
        node_end = node_start + node.length

        # Case 1: Position is in the current node's left subtree
        if abs_position <= node_start:
            if node.left:
                node.left = self._insert_recursive(node.left, abs_position, donor_seq, donor_id, tsd_length)
            else:
                # Insert as left child - 使用预先创建的UID
                node.left = self._create_node(donor_seq, True, donor_id, donor_node_uid)

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

            # 保存当前节点的信息
            old_node_uid = node.uid
            is_current_donor = node.is_donor
            current_donor_id = node.donor_id

            # 创建左、右片段节点的UID
            left_node_uid = self._get_next_uid()
            right_node_uid = self._get_next_uid()

            # 添加左右片段节点到关系图
            self.nesting_graph.add_node(left_node_uid, left_data, current_donor_id, node_start)
            self.nesting_graph.add_node(right_node_uid, right_data, current_donor_id, abs_position + len(donor_seq))

            # 使用预设的UID创建新节点
            new_left = self._create_node(left_data, is_current_donor, current_donor_id, left_node_uid)
            if node.left:
                new_left.left = node.left
                new_left.update()
                # Balance left subtree
                new_left = new_left.balance()

            # 使用预设的UID创建右节点
            new_right = self._create_node(right_data, is_current_donor, current_donor_id, right_node_uid)
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
                    old_node_uid,    # 被切割donor的UID
                    left_node_uid,   # 切割后左片段的UID
                    right_node_uid   # 切割后右片段的UID
                )

            self._release_uid(old_node_uid)

            # 更新节点信息
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
                # Insert as right child - 使用预先创建的UID
                node.right = self._create_node(donor_seq, True, donor_id, donor_node_uid)

            node.update()
            return node.balance()

        raise RuntimeError("[ERROR] Should not reach here")

    def _insert_iterative(self, node: SequenceNode, abs_position: int, donor_seq: str,
                          donor_id: str = None, tsd_length: int = 0, debug: bool = False) -> SequenceNode:
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

        # 创建donor节点UID
        donor_node_uid = self._get_next_uid()
        # 添加到嵌套关系图
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
                    # Create new left child node - 使用预先创建的UID
                    current.left = self._create_node(donor_seq, True, donor_id, donor_node_uid)

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

                # 保存当前节点信息
                old_node_uid = current.uid
                is_current_donor = current.is_donor
                current_donor_id = current.donor_id

                # 创建左、右片段节点UID
                left_node_uid = self._get_next_uid()
                right_node_uid = self._get_next_uid()

                # 添加左右片段节点到关系图
                self.nesting_graph.add_node(left_node_uid, left_data, current_donor_id, node_start)
                self.nesting_graph.add_node(right_node_uid, right_data, current_donor_id, abs_position + len(donor_seq))

                # 使用预设的UID创建左节点
                new_left = self._create_node(left_data, is_current_donor, current_donor_id, left_node_uid)
                if current.left:
                    new_left.left = current.left
                    new_left.update()

                # 使用预设的UID创建右节点
                new_right = self._create_node(right_data, is_current_donor, current_donor_id, right_node_uid)
                if current.right:
                    new_right.right = current.right
                    new_right.update()

                # 如果当前节点是donor，记录切割关系
                if is_current_donor:
                    # 添加切割关系到嵌套关系图
                    self.nesting_graph.add_cut_relation(
                        donor_node_uid,  # 切割者donor的UID
                        old_node_uid,    # 被切割donor的UID
                        left_node_uid,   # 切割后左片段的UID
                        right_node_uid   # 切割后右片段的UID
                    )

                self._release_uid(old_node_uid)

                # 更新节点信息
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
                    # Create new right child node - 使用预先创建的UID
                    current.right = self._create_node(donor_seq, True, donor_id, donor_node_uid)

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
        nodes, edges = self.__build_graphviz_dot_nodes_edges(self.root, node_id_prefix)

        # Add all nodes and edges to the DOT string
        for node in nodes:
            dot_str.append(f"  {node}")
        for edge in edges:
            dot_str.append(f"  {edge}")

        dot_str.append('}')
        return '\n'.join(dot_str)

    def __build_graphviz_dot_nodes_edges(self, node: SequenceNode, node_id_prefix: str, abs_pos: int = 0) -> tuple:
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
        cut_half = ""

        # 重新设计颜色赋值逻辑，使其更加明确
        # 1. 先设置基础颜色（donor为蓝色，acceptor为绿色）
        # 2. 如果是切割片段，改为粉色
        # 3. 如果有嵌套关系，改为黄色
        # 4. 如果同时满足2和3，则为紫色

        # Check if this is a fragment of a cut donor
        if self.nesting_graph.is_fragment(node.uid):
            # Get fragment info (original_uid, position, is_left)
            fragment_info = self.nesting_graph.fragments.get(node.uid)
            if fragment_info:
                orig_uid, _, is_left = fragment_info
                half_type = "L" if is_left else "R"
                cut_half = f"Cut: {half_type}\\n"
                fill_color = "lightpink"  # Cut fragments shown in pink

        # Create node label with position information
        label = "".join([node_type, " | ", str(node.uid), "\\n",
                         str(start_pos_1based), "\\l",
                         str(end_pos_1based), "\\l",
                         "Length: ", str(node.length), "\\n",
                         cut_half])

        # Add the node to the nodes list
        nodes.append(f'{node_id} [label="{label}", fillcolor="{fill_color}"];')

        # Process left child if exists
        if node.left:
            # Left child should start at the same absolute position as its parent
            left_abs_pos = abs_pos
            left_nodes, left_edges = self.__build_graphviz_dot_nodes_edges(
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
            right_nodes, right_edges = self.__build_graphviz_dot_nodes_edges(
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
    表示donor序列和切割关系的图数据结构。
    用于高效追踪和重建被切割的donor序列。
    """
    def __init__(self):
        # 节点信息：UID -> {sequence, length, donor_id, positions}
        self.nodes = {}

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

    def __str__(self) -> str:
        """
        打印图的所有属性信息，用于调试

        Returns:
            str: 格式化的图信息字符串
        """
        lines = []

        # 标题
        lines.append("=" * 32)
        lines.append("DonorNestingGraph 信息")
        lines.append("=" * 32)

        # 节点信息
        lines.append(f"\n节点数量: {len(self.nodes)}")
        for uid, info in self.nodes.items():
            seq = info.get('sequence', '')
            seq_display = seq if len(seq) <= 20 else f"{seq[:10]}...{seq[-10:]}"
            lines.append(f"  UID {uid}: 长度={info.get('length', 0)}, "
                         f"donor_id={info.get('donor_id', 'None')}, "
                         f"位置={info.get('position', 'None')}, "
                         f"序列={seq_display}")

        # 切割关系
        cut_count = sum(len(cuts) for cuts in self.cuts.values())
        lines.append(f"\n切割关系数量: {cut_count}")
        for cutter, cut_list in self.cuts.items():
            if cut_list:
                for cut_info in cut_list:
                    cut_uid, left_uid, right_uid = cut_info
                    lines.append(f"  切割者 {cutter} 切割了 {cut_uid} 产生片段: 左={left_uid}, 右={right_uid}")

        # 反向切割关系
        cut_by_count = sum(len(cuts) for cuts in self.cut_by.values())
        lines.append(f"\n反向切割关系数量: {cut_by_count}")
        for cut, cutter_list in self.cut_by.items():
            if cutter_list:
                for cut_info in cutter_list:
                    cutter_uid, left_uid, right_uid = cut_info
                    lines.append(f"  被切割者 {cut} 被 {cutter_uid} 切割产生片段: 左={left_uid}, 右={right_uid}")

        # 片段关系
        lines.append(f"\n片段关系数量: {len(self.fragments)}")
        for frag_uid, frag_info in self.fragments.items():
            orig_uid, pos, is_left = frag_info
            side = "左" if is_left else "右"
            lines.append(f"  片段 {frag_uid} 来自 {orig_uid} 的{side}侧, 位置={pos}")

        # 片段-原始映射
        lines.append(f"\n片段-原始映射数量: {len(self.fragment_to_original)}")
        for frag_uid, orig_uid in self.fragment_to_original.items():
            lines.append(f"  片段 {frag_uid} -> 原始 {orig_uid}")

        return "\n".join(lines)

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

    def get_cut_by(self, cut_uid: int) -> list:
        """获取切割了指定donor的所有donor信息"""
        return self.cut_by.get(cut_uid, [])

    def get_cuts(self, cutter_uid: int) -> list:
        """获取被指定donor切割的所有donor信息"""
        return self.cuts.get(cutter_uid, [])

    def reconstruct_donors(self, seq_id: str) -> tuple:
        """
        改进的donor重建算法，处理多重切割情况。使用迭代方法代替递归，以更好地处理循环依赖。

        Args:
            seq_id: 序列标识符，用于构建输出记录ID

        Returns:
            tuple: (重建的donor记录列表, 排除的UID集合)
        """
        reconstructed = []
        excluded = set()

        # 收集所有切割关系
        cut_relations = {}  # donor_uid -> [(cutter_uid, left_uid, right_uid, left_seq, right_seq, cutter_seq)]

        # 预处理所有切割关系
        for donor_uid, cut_info_list in self.cut_by.items():
            cut_relations[donor_uid] = []
            for cutter_uid, left_uid, right_uid in cut_info_list:
                # 验证节点存在
                if not all(uid in self.nodes for uid in [donor_uid, cutter_uid, left_uid, right_uid]):
                    continue

                try:
                    left_seq = self.nodes[left_uid]['sequence']
                    right_seq = self.nodes[right_uid]['sequence']
                    cutter_seq = self.nodes[cutter_uid]['sequence']
                    cut_relations[donor_uid].append((cutter_uid, left_uid, right_uid, left_seq, right_seq, cutter_seq))
                except KeyError:
                    continue

        # 待处理donor队列
        to_process = list(cut_relations.keys())

        # 优先处理被切割次数多的donor
        to_process.sort(key=lambda uid: len(cut_relations.get(uid, [])), reverse=True)

        while to_process:
            # 获取下一个待处理的donor
            donor_uid = to_process.pop(0)

            # 如果已处理过，跳过
            if donor_uid in excluded:
                continue

            # 获取有效的切割关系
            relations = cut_relations.get(donor_uid, [])
            if not relations:
                continue

            # 筛选出cutter尚未被处理的切割关系
            valid_relations = [r for r in relations if r[0] not in excluded]

            # 如果没有有效切割关系，视为处理过
            if not valid_relations:
                excluded.add(donor_uid)
                continue

            # 判断是否存在多重切割
            is_multiple_cuts = len(valid_relations) > 1

            # 处理每个有效切割关系，生成完整重建
            for relation in valid_relations:
                cutter_uid, left_uid, right_uid, left_seq, right_seq, cutter_seq = relation

                # 创建完整重建
                full_seq = left_seq + cutter_seq + right_seq
                full_id = f"{seq_id}_reconstructed_{donor_uid}_{cutter_uid}"
                full_rec = create_sequence_record(full_seq, full_id)
                full_rec.annotations["reconstruction_type"] = "full"
                full_rec.annotations["original_uid"] = donor_uid
                full_rec.annotations["cutter_uid"] = cutter_uid
                full_rec.annotations["left_uid"] = left_uid
                full_rec.annotations["right_uid"] = right_uid

                # 将完整重建添加到结果列表
                reconstructed.append(full_rec)

            # 使用第一个切割关系创建清洁重建
            first_relation = valid_relations[0]
            cutter_uid, left_uid, right_uid, left_seq, right_seq, cutter_seq = first_relation

            clean_seq = left_seq + right_seq
            clean_id = f"{seq_id}_clean_reconstructed_{donor_uid}"
            clean_rec = create_sequence_record(clean_seq, clean_id)
            clean_rec.annotations["reconstruction_type"] = "clean"
            clean_rec.annotations["original_uid"] = donor_uid
            clean_rec.annotations["left_uid"] = left_uid
            clean_rec.annotations["right_uid"] = right_uid

            # 多重切割标记
            if is_multiple_cuts:
                clean_rec.annotations["multiple_cuts"] = True
                clean_rec.annotations["cut_count"] = len(valid_relations)
                clean_rec.annotations["all_cutters"] = [r[0] for r in valid_relations]

            # 将清洁重建添加到结果列表
            reconstructed.append(clean_rec)

            # 标记当前donor为已处理
            excluded.add(donor_uid)

            # 如果还存在循环依赖，需要特殊处理
            # 先找出待处理列表中互相切割的donor对
            cyclic_pairs = []
            for i, uid1 in enumerate(to_process):
                for uid2 in to_process[i+1:]:
                    # 检查uid1是否切割uid2
                    cuts_forward = any(r[0] == uid1 for r in cut_relations.get(uid2, []))
                    # 检查uid2是否切割uid1
                    cuts_backward = any(r[0] == uid2 for r in cut_relations.get(uid1, []))

                    if cuts_forward and cuts_backward:
                        cyclic_pairs.append((uid1, uid2))

            # 对每个循环依赖对，强制处理其切割关系
            for uid1, uid2 in cyclic_pairs:
                # 如果已处理过，跳过
                if uid1 in excluded or uid2 in excluded:
                    continue

                # 处理 uid1 切割 uid2
                for relation in cut_relations.get(uid2, []):
                    if relation[0] != uid1:
                        continue

                    cutter_uid, left_uid, right_uid, left_seq, right_seq, cutter_seq = relation

                    # 创建完整重建 (uid1切割uid2)
                    full_seq = left_seq + cutter_seq + right_seq
                    full_id = f"{seq_id}_reconstructed_{uid2}_{uid1}_CYCLIC"
                    full_rec = create_sequence_record(full_seq, full_id)
                    full_rec.annotations["reconstruction_type"] = "full"
                    full_rec.annotations["original_uid"] = uid2
                    full_rec.annotations["cutter_uid"] = uid1
                    full_rec.annotations["left_uid"] = left_uid
                    full_rec.annotations["right_uid"] = right_uid
                    full_rec.annotations["cyclic"] = True

                    # 将完整重建添加到结果列表
                    reconstructed.append(full_rec)

                # 处理 uid2 切割 uid1
                for relation in cut_relations.get(uid1, []):
                    if relation[0] != uid2:
                        continue

                    cutter_uid, left_uid, right_uid, left_seq, right_seq, cutter_seq = relation

                    # 创建完整重建 (uid2切割uid1)
                    full_seq = left_seq + cutter_seq + right_seq
                    full_id = f"{seq_id}_reconstructed_{uid1}_{uid2}_CYCLIC"
                    full_rec = create_sequence_record(full_seq, full_id)
                    full_rec.annotations["reconstruction_type"] = "full"
                    full_rec.annotations["original_uid"] = uid1
                    full_rec.annotations["cutter_uid"] = uid2
                    full_rec.annotations["left_uid"] = left_uid
                    full_rec.annotations["right_uid"] = right_uid
                    full_rec.annotations["cyclic"] = True

                    # 将完整重建添加到结果列表
                    reconstructed.append(full_rec)

        return reconstructed, excluded

    def to_graphviz_dot(self) -> str:
        """
        生成Graphviz DOT格式可视化

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

            # 安全获取节点长度
            length = info.get('length', 0)
            donor_id = info.get('donor_id', '')

            # 创建标签，包含donor_id信息
            label = f"{node_type}\\nUID: {uid}\\nLen: {length}"
            if donor_id:
                label += f"\\nID: {donor_id}"
            label += fragment_info

            # 设置颜色
            # 默认donor为蓝色，其他为绿色
            # 切割片段为粉色
            color = "#AAAAFF" if node_type == "Donor" else "#AAFFAA"

            # 检查是否为切割片段
            if fragment_info:
                color = "#FFAAAA"  # 切割片段为粉色

            lines.append(f'  node_{uid} [label="{label}", fillcolor="{color}"];')

        # 添加切割关系边
        for cutter_uid, cut_list in self.cuts.items():
            for cut_uid, left_uid, right_uid in cut_list:
                lines.append(f'  node_{cutter_uid} -> node_{cut_uid} [label="cuts", color="red", penwidth=2.0];')
                lines.append(f'  node_{cut_uid} -> node_{left_uid} [label="left_frag", style="dashed", color="blue"];')
                lines.append(f'  node_{cut_uid} -> node_{right_uid} [label="right_frag", style="dashed", color="blue"];')

        lines.append("}")
        return '\n'.join(lines)
