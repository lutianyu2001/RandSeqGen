#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from typing import Dict, List, Tuple, Set, Optional, Any, NamedTuple
from Bio.SeqRecord import SeqRecord
from utils import create_sequence_record


# Define InsertionEvent as a namedtuple for simplicity and efficiency
class InsertionEvent(NamedTuple):
    """Represents a sequence insertion event in the event sourcing architecture"""
    event_id: int
    donor_uid: int
    target_uid: int
    left_uid: int
    right_uid: int
    tsd_info: Dict[str, Any] = None  # Optional TSD information
    
    def __str__(self) -> str:
        """String representation for debugging"""
        tsd_str = f", TSD({self.tsd_info['length']}bp)" if self.tsd_info else ""
        return f"Event#{self.event_id}: Donor({self.donor_uid}) → Target({self.target_uid}) → [L({self.left_uid}), R({self.right_uid})]{tsd_str}"


class SequenceEventJournal:
    """
    Event sourcing architecture for sequence insertion events.
    Records all sequence insertion events and supports sequence reconstruction from event history.
    """
    def __init__(self, tree_ref):
        """
        Initialize the event log.

        Args:
            tree_ref: Reference to the host SequenceTree
        """
        self.tree_ref = tree_ref
        self.events = []  # Event log list
        self.next_event_id = 0  # Event ID counter
        self.node_events = {}  # Node event index (uid -> [event_ids])

    def __str__(self) -> str:
        """String representation for debugging"""
        lines = ["=" * 32,
                 "SequenceEventJournal Information",
                 "=" * 32,
                 f"\nEvent count: {len(self.events)}"]

        # Add information for each event
        for event in self.events:
            lines.append(f"  {event}")

        # Add node event index information
        lines.append(f"\nNode event index count: {len(self.node_events)}")
        for uid, event_ids in self.node_events.items():
            lines.append(f"  Node {uid} events: {event_ids}")

        return "\n".join(lines)

    def record_insertion(self, donor_uid: int, target_uid: int, 
                         left_uid: int, right_uid: int, 
                         tsd_info: Dict[str, Any] = None) -> InsertionEvent:
        """
        Record an insertion event.

        Args:
            donor_uid: Inserted donor node UID
            target_uid: Target node UID being inserted into
            left_uid: Left fragment UID after splitting
            right_uid: Right fragment UID after splitting
            tsd_info: TSD information (optional)

        Returns:
            InsertionEvent: The recorded event
        """
        # Create new event
        event = InsertionEvent(
            event_id=self.next_event_id,
            donor_uid=donor_uid,
            target_uid=target_uid,
            left_uid=left_uid,
            right_uid=right_uid,
            tsd_info=tsd_info or {}
        )

        # Add to event list
        self.events.append(event)

        # Update node event index
        for uid in [donor_uid, target_uid, left_uid, right_uid]:
            if uid not in self.node_events:
                self.node_events[uid] = []
            self.node_events[uid].append(self.next_event_id)

        # Increment event ID counter
        self.next_event_id += 1

        return event

    def is_donor(self, uid: int) -> bool:
        """
        Determine if a node is a donor (appears as donor_uid in events).

        Args:
            uid: Node UID

        Returns:
            bool: True if the node is a donor
        """
        # First check if node has the donor flag in the tree
        node = self.tree_ref.node_dict.get(uid)
        if node and node.is_donor:
            return True
            
        # Also check if it appears as a donor in events
        for event in self.events:
            if event.donor_uid == uid:
                return True
                
        return False

    def is_fragment(self, uid: int) -> bool:
        """
        Determine if a node is a fragment (appears as left_uid or right_uid in events).

        Args:
            uid: Node UID

        Returns:
            bool: True if the node is a fragment
        """
        # Skip uid 0 which is a sentinel value
        if uid == 0:
            return False
            
        for event in self.events:
            if event.left_uid == uid or event.right_uid == uid:
                return True
        return False

    def get_fragment_info(self, uid: int) -> Optional[Tuple[int, bool, int]]:
        """
        Get complete information about a fragment node.

        Args:
            uid: Fragment node UID

        Returns:
            Optional[Tuple[int, bool, int]]: (original node UID, is left fragment, cutter UID)
                                             or None if not found
        """
        for event in self.events:
            if event.left_uid == uid:
                return (event.target_uid, True, event.donor_uid)
            if event.right_uid == uid:
                return (event.target_uid, False, event.donor_uid)
        return None

    def find_node_events(self, uid: int) -> List[InsertionEvent]:
        """
        Find all events a node participates in.

        Args:
            uid: Node UID

        Returns:
            List[InsertionEvent]: List of events
        """
        event_ids = self.node_events.get(uid, [])
        return [self.events[event_id] for event_id in event_ids]

    def find_fragment_creation_event(self, uid: int) -> Optional[InsertionEvent]:
        """
        Find the event that created a specific fragment.

        Args:
            uid: Fragment node UID

        Returns:
            Optional[InsertionEvent]: Creation event, or None if not found
        """
        # Skip uid 0 which is a sentinel value for direct insertion
        if uid == 0:
            return None
            
        for event in self.events:
            if event.left_uid == uid or event.right_uid == uid:
                return event
        return None

    def find_target_events(self, uid: int) -> List[InsertionEvent]:
        """
        Find all insertion events targeting a node.

        Args:
            uid: Target node UID

        Returns:
            List[InsertionEvent]: List of events, sorted by event ID
        """
        target_events = []
        for event in self.events:
            if event.target_uid == uid:
                target_events.append(event)
        
        # Sort by event ID (chronological order)
        target_events.sort(key=lambda e: e.event_id)
        
        return target_events

    def reconstruct(self, uid: int, seq_id: str) -> Dict[str, Any]:
        """
        Comprehensive reconstruction algorithm, generating three types of reconstructions in one pass.

        Args:
            uid: Node UID to reconstruct
            seq_id: Sequence identifier for output record IDs

        Returns:
            Dict[str, Any]: Dictionary with three reconstruction types:
                - 'full': Complete sequence with all nested structures
                - 'clean': Original sequence without nested structures
                - 'events': Sequence state after each event
        """
        # Get node sequence
        node_seq = self._get_node_sequence(uid)
        if node_seq is None:
            return None

        # Initialize reconstruction results
        result = {
            'full': node_seq,           # Complete reconstruction with nesting
            'clean': node_seq,          # Clean reconstruction without nesting
            'events': [(None, node_seq)]  # Event history [(event_id, sequence), ...]
        }

        # Determine node type and handle reconstruction accordingly
        if self.is_fragment(uid):
            # Fragment node: recursively process the original node
            fragment_result = self._reconstruct_fragment(uid, seq_id)
            if fragment_result:
                result = fragment_result
        else:
            # Check if it's a target node (has been inserted into)
            target_events = self.find_target_events(uid)
            if target_events:
                # Target node: process all insertion events
                result = self._reconstruct_target(uid, target_events, seq_id)

        return result

    def _get_node_sequence(self, uid: int) -> Optional[str]:
        """Get the original sequence of a node"""
        node = self.tree_ref.node_dict.get(uid)
        return node.data if node else None

    def _reconstruct_fragment(self, fragment_uid: int, seq_id: str) -> Dict[str, Any]:
        """
        Reconstruct a fragment node.

        Args:
            fragment_uid: Fragment node UID
            seq_id: Sequence identifier

        Returns:
            Dict[str, Any]: Reconstruction results
        """
        # Find the event that created this fragment
        event = self.find_fragment_creation_event(fragment_uid)
        if not event:
            return None

        # Get original node information
        original_uid = event.target_uid
        is_left_fragment = (event.left_uid == fragment_uid)

        # Recursively reconstruct the original node
        original_result = self.reconstruct(original_uid, seq_id)
        if not original_result:
            return None

        # Get fragment sequence from the original node's reconstruction
        original_seq = original_result['full']
        original_clean_seq = original_result['clean']
        
        # Get the current node sequence
        node_seq = self._get_node_sequence(fragment_uid)
        
        # Extract the appropriate part from the original node's reconstruction
        if is_left_fragment:
            # Left fragment: take the first part of the original sequence
            fragment_pos = len(node_seq)
            fragment_full = original_seq[:fragment_pos]
            fragment_clean = original_clean_seq[:fragment_pos]
        else:
            # Right fragment: take the latter part of the original sequence
            donor_seq = self._get_node_sequence(event.donor_uid)
            left_seq = self._get_node_sequence(event.left_uid)
            fragment_pos = len(left_seq) + len(donor_seq)
            fragment_full = original_seq[fragment_pos:]
            fragment_clean = original_clean_seq[fragment_pos:]

        # Build event history reconstruction
        fragment_events = []
        for event_id, seq in original_result['events']:
            if event_id is None:
                # Initial state
                if is_left_fragment:
                    fragment_events.append((None, seq[:len(node_seq)]))
                else:
                    fragment_events.append((None, seq[-len(node_seq):]))
            else:
                # State after event
                if is_left_fragment:
                    fragment_len = min(len(node_seq), len(seq))
                    fragment_events.append((event_id, seq[:fragment_len]))
                else:
                    # Right fragment may need more complex calculation
                    e = self.events[event_id]
                    if e.target_uid == original_uid:
                        # If the event affects the original node, calculate right fragment position
                        donor_len = len(self._get_node_sequence(e.donor_uid))
                        left_len = len(self._get_node_sequence(e.left_uid))
                        right_pos = left_len + donor_len
                        fragment_events.append((event_id, seq[right_pos:]))
                    else:
                        # Otherwise, keep the same right fragment
                        fragment_events.append((event_id, fragment_events[-1][1]))

        return {
            'full': fragment_full,
            'clean': fragment_clean,
            'events': fragment_events
        }

    def _reconstruct_target(self, target_uid: int, target_events: List[InsertionEvent], seq_id: str) -> Dict[str, Any]:
        """
        Reconstruct a target node (node that was inserted into).

        Args:
            target_uid: Target node UID
            target_events: List of events targeting this node
            seq_id: Sequence identifier

        Returns:
            Dict[str, Any]: Reconstruction results
        """
        # Initial state is the node's original sequence
        original_seq = self._get_node_sequence(target_uid)
        if original_seq is None:
            return None

        # Initialize results
        full_seq = original_seq
        clean_seq = original_seq
        event_seqs = [(None, original_seq)]

        # Process each event in chronological order
        for event in target_events:
            # Recursively reconstruct the donor node
            donor_result = self.reconstruct(event.donor_uid, seq_id)
            donor_full_seq = donor_result['full'] if donor_result else self._get_node_sequence(event.donor_uid)
            
            # Get left and right fragment sequences
            left_seq = self._get_node_sequence(event.left_uid)
            right_seq = self._get_node_sequence(event.right_uid)
            
            # Apply TSD if present
            if event.tsd_info:
                tsd_5 = event.tsd_info.get('tsd_5', '')
                tsd_3 = event.tsd_info.get('tsd_3', '')
                donor_with_tsd = tsd_5 + donor_full_seq + tsd_3
            else:
                donor_with_tsd = donor_full_seq
            
            # Update full reconstruction (with nesting)
            full_seq = left_seq + donor_with_tsd + right_seq
            
            # Update clean reconstruction (without nesting)
            clean_seq = left_seq + right_seq
            
            # Update event history
            event_seqs.append((event.event_id, full_seq))

        return {
            'full': full_seq,
            'clean': clean_seq,
            'events': event_seqs
        }

    def get_reconstructed_donor_uids(self) -> Set[int]:
        """
        Get all donor UIDs that need reconstruction.
        These are donors that have been split by other donors.

        Returns:
            Set[int]: Set of UIDs of donors to reconstruct
        """
        reconstructed_uids = set()
        
        for event in self.events:
            # If a node is a target and it's a donor (not the root node)
            if event.target_uid != self.tree_ref.root.uid and self.is_donor(event.target_uid):
                reconstructed_uids.add(event.target_uid)
        
        return reconstructed_uids

    def reconstruct_donors_to_records(self, seq_id: str) -> List[SeqRecord]:
        """
        Reconstruct nested sequences and create SeqRecord objects.

        Args:
            seq_id: Sequence identifier

        Returns:
            List[SeqRecord]: List of reconstructed records
        """
        reconstructed_records = []
        
        # Get all donors that need reconstruction
        reconstructed_uids = self.get_reconstructed_donor_uids()
        
        for uid in reconstructed_uids:
            # Reconstruct sequence
            reconstruct_result = self.reconstruct(uid, seq_id)
            if not reconstruct_result:
                continue
            
            # Create full reconstruction record
            full_seq = reconstruct_result['full']
            full_id = f"{seq_id}_reconstructed_{uid}"
            full_rec = create_sequence_record(full_seq, full_id)
            full_rec.annotations["reconstruction_type"] = "full"
            full_rec.annotations["original_uid"] = uid
            
            # Create clean reconstruction record
            clean_seq = reconstruct_result['clean']
            clean_id = f"{seq_id}_clean_reconstructed_{uid}"
            clean_rec = create_sequence_record(clean_seq, clean_id)
            clean_rec.annotations["reconstruction_type"] = "clean"
            clean_rec.annotations["original_uid"] = uid
            
            # Add to result list
            reconstructed_records.append(full_rec)
            reconstructed_records.append(clean_rec)
        
        return reconstructed_records

    def collect_donor_records(self, seq_id: str) -> List[SeqRecord]:
        """
        Collect donor records from the tree.

        Args:
            seq_id: Sequence identifier

        Returns:
            List[SeqRecord]: Donor record list
        """
        donor_records = []
        abs_position_map = self._calculate_absolute_positions()
        
        # Iterate through all nodes, find donors
        for uid, node in self.tree_ref.node_dict.items():
            if node.is_donor:
                # Get absolute position information
                start_pos = abs_position_map.get(uid, 0)
                end_pos = start_pos + node.length
                
                # Convert to 1-based index for output
                start_pos_1based = start_pos + 1
                end_pos_1based = end_pos
                
                # Create donor ID
                donor_id = f"{seq_id}_{start_pos_1based}_{end_pos_1based}-+-{node.length}"
                if node.donor_id:
                    donor_id += f"-{node.donor_id}"
                
                # Create record
                donor_record = create_sequence_record(node.data, donor_id)
                donor_record.annotations["uid"] = uid
                donor_record.annotations["position"] = start_pos_1based
                donor_record.annotations["length"] = node.length
                
                # Add to result list
                donor_records.append(donor_record)
        
        return donor_records

    def to_graphviz_dot(self) -> str:
        """
        Generate GraphViz DOT format graph of events and node relationships.

        Returns:
            str: Graphviz DOT format string
        """
        if not self.tree_ref.root:
            return "digraph SequenceEventJournal { }"
        
        # Initialize DOT header
        dot_str = ["digraph SequenceEventJournal {",
                  "  bgcolor=\"#FFFFFF\";",
                  "  node [fontcolor=\"#000000\", shape=box, style=filled];",
                  "  edge [fontcolor=\"#000000\", penwidth=2.0];",
                  "  rankdir=LR;"]  # Left-to-right layout, better for event flow
        
        # Add root node if it exists
        root = self.tree_ref.root
        if root:
            root_label = f"Root\\nUID: {root.uid}\\nLen: {len(root.data)}"
            dot_str.append(f'  node_{root.uid} [label="{root_label}", fillcolor="lightgreen"];')
        
        # If there are no events, we still need at least one node to pass the test
        if not self.events and root:
            dot_str.append('}')
            return '\n'.join(dot_str)
        
        # Collect nodes used in events
        nodes = set()
        for event in self.events:
            nodes.add(event.donor_uid)
            nodes.add(event.target_uid)
            nodes.add(event.left_uid)
            nodes.add(event.right_uid)
        
        # Generate representation for each node
        for uid in nodes:
            node = self.tree_ref.node_dict.get(uid)
            if not node:
                continue
            
            # Determine node type and color
            if self.is_donor(uid):
                node_type = "Donor"
                fill_color = "lightblue"
            elif self.is_fragment(uid):
                fragment_info = self.get_fragment_info(uid)
                if fragment_info:
                    original_uid, is_left, cutter_uid = fragment_info
                    if is_left:
                        node_type = "Left Fragment"
                        fill_color = "lightpink"
                    else:
                        node_type = "Right Fragment"
                        fill_color = "salmon"
                else:
                    node_type = "Fragment"
                    fill_color = "lightpink"
            else:
                node_type = "Node"
                fill_color = "lightgreen"
            
            # Create node label
            label = f"{node_type}\\nUID: {uid}\\nLen: {len(node.data)}"
            
            # Add node to graph
            dot_str.append(f'  node_{uid} [label="{label}", fillcolor="{fill_color}"];')
        
        # Add event nodes and relationships
        for event in self.events:
            # Add event node
            event_label = f"Event #{event.event_id}"
            if event.tsd_info:
                event_label += f"\\nTSD: {event.tsd_info.get('length', 0)}bp"
            
            dot_str.append(f'  event_{event.event_id} [label="{event_label}", shape=ellipse, fillcolor="lightyellow"];')
            
            # Add relationship edges
            # Donor to event
            dot_str.append(f'  node_{event.donor_uid} -> event_{event.event_id} [label="donor", color="blue"];')
            # Event to target
            dot_str.append(f'  event_{event.event_id} -> node_{event.target_uid} [label="target", color="red"];')
            # Target to left/right fragments
            dot_str.append(f'  node_{event.target_uid} -> node_{event.left_uid} [label="left", style="dashed", color="darkgreen"];')
            dot_str.append(f'  node_{event.target_uid} -> node_{event.right_uid} [label="right", style="dashed", color="darkgreen"];')
        
        dot_str.append('}')
        return '\n'.join(dot_str)

    def _calculate_absolute_positions(self) -> Dict[int, int]:
        """
        Calculate absolute positions of all nodes in the tree.

        Returns:
            Dict[int, int]: Mapping from node UID to absolute position
        """
        positions = {}

        def _traverse(node, current_pos=0):
            if not node:
                return current_pos

            # Process left subtree
            left_end_pos = _traverse(node.left, current_pos)

            # Calculate current node position
            node_pos = left_end_pos
            positions[node.uid] = node_pos

            # Process right subtree
            right_end_pos = _traverse(node.right, node_pos + node.length)

            return right_end_pos

        _traverse(self.tree_ref.root, 0)
        return positions 