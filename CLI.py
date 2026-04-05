"""
CLI entry point for DNA Analyzer.
"""

import argparse
import sys

from dna_analyzer.analyzer import (
    validate_sequence,
    clean_sequence,
    find_orfs,
    find_motif,
    sequence_summary,
    reverse_complement,
)
from dna_analyzer.visualizer import (
    print_header,
    visualize_nucleotide_freq,
    visualize_gc_content,
    visualize_orfs,
    visualize_dinucleotide,
    visualize_sequence_map,
    visualize_summary,
    RESET, DIM, BOLD, CYAN, WHITE,
)


def read_sequence(args) -> str:
    """Read sequence from argument, file, or stdin."""
    if args.sequence:
        return args.sequence
    if args.file:
        try:
            with open(args.file, "r") as f:
                content = f.read()
            # Basic FASTA stripping
            lines = content.splitlines()
            seq_lines = [l for l in lines if not l.startswith(">")]
            return "".join(seq_lines)
        except FileNotFoundError:
            print(f"\n  ✖  File not found: {args.file}", file=sys.stderr)
            sys.exit(1)
    if not sys.stdin.isatty():
        return sys.stdin.read()
    return ""


def cmd_analyze(args):
    """Full analysis pipeline."""
    raw = read_sequence(args)
    if not raw:
        print("\n  ✖  No sequence provided. Use -s, -f, or pipe via stdin.\n", file=sys.stderr)
        sys.exit(1)

    valid, err = validate_sequence(raw)
    if not valid:
        print(f"\n  ✖  Invalid sequence: {err}\n", file=sys.stderr)
        sys.exit(1)

    seq = clean_sequence(raw)
    print_header("🧬  DNA SEQUENCE ANALYZER")

    summary = sequence_summary(seq)
    visualize_summary(summary)
    visualize_sequence_map(seq)
    visualize_nucleotide_freq(seq)
    visualize_gc_content(seq)

    orfs = find_orfs(seq, min_length=args.min_orf)
    visualize_orfs(orfs, len(seq))

    if not args.no_dinuc:
        visualize_dinucleotide(seq)

    # Print reverse complement
    print(f"\n  {DIM}Reverse Complement:{RESET}")
    rc = reverse_complement(seq)
    print(f"  {DIM}{rc[:80]}{'...' if len(rc) > 80 else ''}{RESET}")

    print()


def cmd_motif(args):
    """Find motif positions."""
    raw = read_sequence(args)
    valid, err = validate_sequence(raw)
    if not valid:
        print(f"\n  ✖  {err}\n", file=sys.stderr)
        sys.exit(1)

    motif_valid, motif_err = validate_sequence(args.motif)
    if not motif_valid:
        print(f"\n  ✖  Invalid motif: {motif_err}\n", file=sys.stderr)
        sys.exit(1)

    seq = clean_sequence(raw)
    positions = find_motif(seq, args.motif)

    print_header(f"🔍  MOTIF SEARCH: {args.motif.upper()}")
    if positions:
        print(f"  {BOLD}{WHITE}Found {len(positions)} occurrence(s):{RESET}\n")
        for pos in positions:
            snippet_start = max(0, pos - 4)
            snippet_end = min(len(seq), pos + len(args.motif) + 3)
            before = seq[snippet_start:pos - 1]
            match = seq[pos - 1:pos - 1 + len(args.motif)]
            after = seq[pos - 1 + len(args.motif):snippet_end]
            print(f"  Position {BOLD}{CYAN}{pos:>5}{RESET}  …{DIM}{before}{RESET}"
                  f"\033[93m\033[1m{match}{RESET}{DIM}{after}{RESET}…")
    else:
        print(f"  {DIM}No occurrences of '{args.motif.upper()}' found.{RESET}")
    print()


def cmd_orfs(args):
    """List ORFs only."""
    raw = read_sequence(args)
    valid, err = validate_sequence(raw)
    if not valid:
        print(f"\n  ✖  {err}\n", file=sys.stderr)
        sys.exit(1)

    seq = clean_sequence(raw)
    orfs = find_orfs(seq, min_length=args.min_orf)

    print_header("🧬  OPEN READING FRAMES")
    visualize_orfs(orfs, len(seq))

    if orfs and args.proteins:
        print(f"\n  {BOLD}Translated Proteins:{RESET}\n")
        for i, orf in enumerate(orfs[:5], 1):
            print(f"  {DIM}#{i} Frame {orf['frame']} ({orf['length']} bp):{RESET}")
            print(f"     {CYAN}{orf['protein'][:80]}{RESET}\n")
    print()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dna-analyzer",
        description="🧬 DNA Sequence Analyzer — bioinformatics in your terminal",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  dna-analyzer analyze -s "ATGCGATCGATCGAAATTTGCGATCG"
  dna-analyzer analyze -f genome.fasta --min-orf 60
  dna-analyzer motif -s "ATGCGATCGATG" -m "ATG"
  dna-analyzer orfs -f sequence.fasta --proteins
  cat sequence.txt | dna-analyzer analyze
        """,
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Shared sequence args
    def add_seq_args(p):
        group = p.add_mutually_exclusive_group()
        group.add_argument("-s", "--sequence", metavar="SEQ", help="DNA sequence string")
        group.add_argument("-f", "--file", metavar="FILE", help="FASTA or plain text file")

    # analyze
    p_analyze = subparsers.add_parser("analyze", help="Full sequence analysis")
    add_seq_args(p_analyze)
    p_analyze.add_argument("--min-orf", type=int, default=30, metavar="N",
                           help="Minimum ORF length in bp (default: 30)")
    p_analyze.add_argument("--no-dinuc", action="store_true",
                           help="Skip dinucleotide heatmap")
    p_analyze.set_defaults(func=cmd_analyze)

    # motif
    p_motif = subparsers.add_parser("motif", help="Find motif occurrences")
    add_seq_args(p_motif)
    p_motif.add_argument("-m", "--motif", required=True, metavar="MOTIF",
                         help="Motif/pattern to search for")
    p_motif.set_defaults(func=cmd_motif)

    # orfs
    p_orfs = subparsers.add_parser("orfs", help="Find open reading frames")
    add_seq_args(p_orfs)
    p_orfs.add_argument("--min-orf", type=int, default=30, metavar="N",
                        help="Minimum ORF length in bp (default: 30)")
    p_orfs.add_argument("--proteins", action="store_true",
                        help="Show translated amino acid sequences")
    p_orfs.set_defaults(func=cmd_orfs)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
