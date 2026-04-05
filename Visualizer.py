"""
Terminal visualization for DNA analysis results.
"""

from dna_analyzer.analyzer import nucleotide_frequency, gc_content, dinucleotide_frequency


NUCLEOTIDE_COLORS = {
    "A": "\033[92m",   # green
    "T": "\033[91m",   # red
    "C": "\033[94m",   # blue
    "G": "\033[93m",   # yellow
    "N": "\033[90m",   # grey
}
RESET = "\033[0m"
BOLD = "\033[1m"
DIM = "\033[2m"
CYAN = "\033[96m"
MAGENTA = "\033[95m"
WHITE = "\033[97m"


def _bar(value: float, max_value: float, width: int = 30, char: str = "█") -> str:
    filled = int((value / max_value) * width) if max_value > 0 else 0
    return char * filled + DIM + "░" * (width - filled) + RESET


def print_header(title: str) -> None:
    width = 60
    print()
    print(CYAN + BOLD + "╔" + "═" * (width - 2) + "╗" + RESET)
    padding = (width - 2 - len(title)) // 2
    print(CYAN + BOLD + "║" + " " * padding + WHITE + title + " " * (width - 2 - padding - len(title)) + CYAN + "║" + RESET)
    print(CYAN + BOLD + "╚" + "═" * (width - 2) + "╝" + RESET)
    print()


def print_section(title: str) -> None:
    print(MAGENTA + BOLD + f"\n  ▶ {title}" + RESET)
    print(DIM + "  " + "─" * 50 + RESET)


def visualize_nucleotide_freq(sequence: str) -> None:
    """Print a colored bar chart of nucleotide frequencies."""
    print_section("Nucleotide Frequency")
    freq = nucleotide_frequency(sequence)
    total = len(sequence.replace(" ", "").replace("\n", ""))
    max_count = max(freq.values()) if freq.values() else 1

    for nuc, count in freq.items():
        if nuc == "N" and count == 0:
            continue
        color = NUCLEOTIDE_COLORS.get(nuc, RESET)
        pct = (count / total * 100) if total > 0 else 0
        bar = _bar(count, max_count)
        print(f"  {color}{BOLD}{nuc}{RESET}  {bar}  {color}{count:>5}{RESET}  {DIM}({pct:5.1f}%){RESET}")


def visualize_gc_content(sequence: str) -> None:
    """Display a GC content meter."""
    print_section("GC Content")
    gc = gc_content(sequence)
    at = 100 - gc

    gc_bar = _bar(gc, 100, width=40, char="▓")
    at_bar = _bar(at, 100, width=40, char="░")

    print(f"  {NUCLEOTIDE_COLORS['G']}G+C{RESET}  {gc_bar}  {BOLD}{gc:5.1f}%{RESET}")
    print(f"  {NUCLEOTIDE_COLORS['A']}A+T{RESET}  {at_bar}  {BOLD}{at:5.1f}%{RESET}")

    # Stability label
    print()
    if gc > 60:
        label = f"{NUCLEOTIDE_COLORS['G']}High GC — thermostable, resists denaturation{RESET}"
    elif gc < 40:
        label = f"{NUCLEOTIDE_COLORS['A']}Low GC — AT-rich, lower melting point{RESET}"
    else:
        label = f"{CYAN}Balanced GC content{RESET}"
    print(f"  {DIM}Stability:{RESET} {label}")


def visualize_orfs(orfs: list[dict], seq_length: int) -> None:
    """Print a map of discovered ORFs."""
    print_section(f"Open Reading Frames  ({len(orfs)} found)")
    if not orfs:
        print(f"  {DIM}No ORFs detected with current minimum length filter.{RESET}")
        return

    for i, orf in enumerate(orfs[:10], 1):  # show top 10
        frame_color = CYAN if orf["frame"].startswith("+") else MAGENTA
        bar_start = int((orf["start"] / seq_length) * 40)
        bar_end = int((orf["end"] / seq_length) * 40)
        bar_len = max(1, bar_end - bar_start)
        map_str = DIM + "·" * bar_start + RESET + frame_color + "█" * bar_len + RESET + DIM + "·" * (40 - bar_start - bar_len) + RESET

        print(f"  {BOLD}#{i:<2}{RESET}  [{map_str}]  "
              f"{frame_color}Frame {orf['frame']}{RESET}  "
              f"{WHITE}{orf['start']:>5}–{orf['end']:<5}{RESET}  "
              f"{DIM}{orf['length']} bp{RESET}")


def visualize_dinucleotide(sequence: str) -> None:
    """Print a 4x4 dinucleotide heatmap."""
    print_section("Dinucleotide Heatmap")
    freq = dinucleotide_frequency(sequence)
    max_val = max(freq.values()) if freq.values() else 1
    nucleotides = "ATCG"

    # Header row
    print("       " + "  ".join(f"{BOLD}{WHITE}{n}{RESET}" for n in nucleotides))
    for a in nucleotides:
        row = f"  {BOLD}{WHITE}{a}{RESET}   "
        for b in nucleotides:
            count = freq.get(f"{a}{b}", 0)
            intensity = count / max_val if max_val > 0 else 0
            if intensity > 0.75:
                color = "\033[97m\033[41m"  # white on red
            elif intensity > 0.5:
                color = "\033[91m"
            elif intensity > 0.25:
                color = "\033[93m"
            elif intensity > 0:
                color = "\033[92m"
            else:
                color = DIM
            row += f"{color}{count:>3}{RESET}  "
        print(row)


def visualize_sequence_map(sequence: str, width: int = 60) -> None:
    """Print a color-coded sequence map (first 120 bases)."""
    print_section("Sequence Map  (first 120 bp, color-coded)")
    seq = sequence.replace(" ", "").replace("\n", "").upper()[:120]

    print("  ", end="")
    for i, base in enumerate(seq):
        color = NUCLEOTIDE_COLORS.get(base, RESET)
        print(f"{color}{base}{RESET}", end="")
        if (i + 1) % width == 0:
            print(f"\n  ", end="")
    print()

    # Legend
    print()
    print("  " + DIM + "Legend: " + RESET + "  ".join(
        f"{NUCLEOTIDE_COLORS[n]}{BOLD}{n}{RESET}" for n in "ATCG"
    ))


def visualize_summary(summary: dict) -> None:
    """Print the full summary table."""
    print_section("Sequence Summary")
    rows = [
        ("Length", f"{summary['length']} bp"),
        ("GC Content", f"{summary['gc_content']}%"),
        ("AT Content", f"{summary['at_content']}%"),
        ("Melting Temp (Tm)", f"{summary['melting_temp']} °C"),
        ("ORFs Found", str(summary["num_orfs"])),
    ]
    for label, value in rows:
        print(f"  {DIM}{label:<22}{RESET} {BOLD}{WHITE}{value}{RESET}")

    counts = summary["nucleotide_counts"]
    print(f"\n  {DIM}{'Nucleotide counts':<22}{RESET}", end="")
    for nuc, cnt in counts.items():
        if nuc == "N" and cnt == 0:
            continue
        color = NUCLEOTIDE_COLORS.get(nuc, RESET)
        print(f" {color}{nuc}:{BOLD}{cnt}{RESET}", end="")
    print()
