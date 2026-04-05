"""
dna-analyzer — A terminal-based DNA sequence analysis toolkit.
"""

__version__ = "1.0.0"
__author__ = "Your Name"
__license__ = "MIT"

from dna_analyzer.analyzer import (
    validate_sequence,
    clean_sequence,
    gc_content,
    nucleotide_frequency,
    reverse_complement,
    find_motif,
    find_orfs,
    melting_temperature,
    sequence_summary,
)

__all__ = [
    "validate_sequence",
    "clean_sequence",
    "gc_content",
    "nucleotide_frequency",
    "reverse_complement",
    "find_motif",
    "find_orfs",
    "melting_temperature",
    "sequence_summary",
]
