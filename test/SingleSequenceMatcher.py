#!/usr/bin/env python3

from Bio.Align import PairwiseAligner
from Bio import SeqIO
from functools import lru_cache

class SingleSequenceMatcher:
    """
    A class for matching a single query sequence against a reference sequence database.
    """

    def __init__(self, min_identity=90.0, min_length=50, kmer_size=7, kmer_threshold=0.1):
        """
        Initialize the SingleSequenceMatcher with matching parameters.

        Args:
            min_identity (float): Minimum sequence identity percentage (default: 90.0)
            min_length (int): Minimum alignment length (default: 50)
            kmer_size (int): Size of k-mers for pre-filtering (default: 7)
            kmer_threshold (float): Minimum k-mer similarity threshold for pre-filtering (default: 0.1)
        """
        self.min_identity = min_identity
        self.min_length = min_length
        self.kmer_size = kmer_size
        self.kmer_threshold = kmer_threshold

        # Validate inputs
        if self.min_identity <= 0 or self.min_identity > 100:
            raise ValueError("min_identity must be between 0 and 100")
        if self.min_length <= 0:
            raise ValueError("min_length must be greater than 0")
        if self.kmer_size <= 0:
            raise ValueError("kmer_size must be greater than 0")
        if self.kmer_threshold < 0 or self.kmer_threshold > 1:
            raise ValueError("kmer_threshold must be between 0 and 1")

        # Configure aligner for speed
        self.aligner = PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
        self.aligner.open_gap_score = -2
        self.aligner.extend_gap_score = -0.5

    def _generate_kmers(self, sequence):
        """
        Generate k-mers from a sequence for fast pre-filtering.

        Args:
            sequence (str): Input sequence

        Returns:
            set: Set of k-mers in the sequence
        """
        return {sequence[i:i+self.kmer_size] for i in range(len(sequence) - self.kmer_size + 1)}

    def _kmer_similarity(self, kmers1, kmers2):
        """
        Calculate Jaccard similarity between two sets of k-mers.

        Args:
            kmers1 (set): First set of k-mers
            kmers2 (set): Second set of k-mers

        Returns:
            float: Jaccard similarity score
        """
        intersection = len(kmers1.intersection(kmers2))
        union = len(kmers1.union(kmers2))
        return intersection / union if union > 0 else 0

    @lru_cache(maxsize=1024)
    def _calculate_sequence_identity(self, seq1, seq2):
        """
        Calculate the identity percentage between two sequences.
        Uses LRU cache to avoid recalculating the same alignments.

        Args:
            seq1 (str): First sequence
            seq2 (str): Second sequence

        Returns:
            tuple: (identity_percentage, alignment_length, alignment_str)
        """
        # Convert sequences to strings to ensure hashability for lru_cache
        seq1, seq2 = str(seq1), str(seq2)

        # Early termination if length difference is too great
        len_diff_ratio = abs(len(seq1) - len(seq2)) / max(len(seq1), len(seq2))
        if len_diff_ratio > 0.3:  # If length difference > 30%, no match is likely
            return 0, 0, ""

        try:
            # Use the faster Bio.Align.PairwiseAligner
            alignment = self.aligner.align(seq1, seq2)[0]

            # Get the full alignment string for visualization
            alignment_str = str(alignment)

            # Extract aligned sequences
            lines = alignment_str.split('\n')

            if len(lines) >= 3:
                seq1_aligned = lines[0]
                seq2_aligned = lines[2]

                # Calculate matching positions
                matches = 0
                align_length = 0

                for a, b in zip(seq1_aligned, seq2_aligned):
                    if a != '-' or b != '-':
                        align_length += 1
                        if a == b and a != '-' and b != '-':
                            matches += 1

                # Calculate identity percentage
                if align_length > 0:
                    identity = (matches / align_length) * 100
                else:
                    identity = 0

                return identity, align_length, alignment_str
            else:
                return 0, 0, ""
        except Exception as e:
            # Handle any alignment errors gracefully
            print(f"Alignment error: {e}")
            return 0, 0, ""

    def find_best_match(self, query_sequence, reference_file):
        """
        Find the best match for a query sequence in a reference sequence database.

        Args:
            query_sequence (str): Query sequence to match
            reference_file (str): Path to reference sequences file (FASTA format)

        Returns:
            tuple: (best_ref_id, best_identity, best_length, alignment_str)
                  If no match is found, returns (None, 0, 0, "")
        """
        # Generate k-mers for query sequence
        query_kmers = self._generate_kmers(query_sequence)

        # Track best match
        best_ref_id = None
        best_identity = 0
        best_length = 0
        best_alignment_str = ""

        # Iterate through reference sequences
        for record in SeqIO.parse(reference_file, "fasta"):
            ref_id = record.id
            ref_sequence = str(record.seq)

            # Generate k-mers for reference sequence
            ref_kmers = self._generate_kmers(ref_sequence)

            # Calculate k-mer similarity for pre-filtering
            kmer_sim = self._kmer_similarity(query_kmers, ref_kmers)
            if kmer_sim < self.kmer_threshold:
                continue

            # Calculate sequence identity and alignment
            identity, align_length, alignment_str = self._calculate_sequence_identity(query_sequence, ref_sequence)

            # Check if it meets threshold requirements
            if identity >= self.min_identity and align_length >= self.min_length:
                # Update best match
                if identity > best_identity or (identity == best_identity and align_length > best_length):
                    best_ref_id = ref_id
                    best_identity = identity
                    best_length = align_length
                    best_alignment_str = alignment_str

        return (best_ref_id, best_identity, best_length, best_alignment_str)


# Example usage
def main():
    """
    Main function for command-line execution.
    """
    import argparse

    parser = argparse.ArgumentParser(description='Find the best match for a query sequence in a reference database')
    parser.add_argument('--query', required=True, help='Query sequence')
    parser.add_argument('--ref-file', required=True, help='Reference sequence file (FASTA format)')
    parser.add_argument('--min-identity', type=float, default=90.0, help='Minimum sequence identity percentage (default: 90.0)')
    parser.add_argument('--min-length', type=int, default=50, help='Minimum alignment length (default: 50)')
    parser.add_argument('--kmer-size', type=int, default=7, help='Size of k-mers for pre-filtering (default: 7)')
    parser.add_argument('--kmer-threshold', type=float, default=0.1, help='Minimum k-mer similarity threshold (default: 0.1)')

    args = parser.parse_args()

    # Create matcher instance
    matcher = SingleSequenceMatcher(
        min_identity=args.min_identity,
        min_length=args.min_length,
        kmer_size=args.kmer_size,
        kmer_threshold=args.kmer_threshold
    )

    # Execute matching
    best_ref_id, best_identity, best_length, alignment_str = matcher.find_best_match(
        query_sequence=args.query,
        reference_file=args.ref_file
    )

    if best_ref_id:
        print(f"\nBest match: {best_ref_id}")
        print(f"Sequence identity: {best_identity:.2f}%")
        print(f"Alignment length: {best_length}")
        print("\nSequence alignment visualization:")
        print(alignment_str)
    else:
        print("\nNo match found that meets the criteria")

    return 0


# For command-line execution
if __name__ == "__main__":
    exit(main())
