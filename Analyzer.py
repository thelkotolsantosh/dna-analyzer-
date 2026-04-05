"""
Core DNA sequence analysis functions.
"""

from collections import Counter
from typing import Optional


COMPLEMENT_MAP = str.maketrans("ATCGatcg", "TAGCtagc")

CODON_TABLE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met (START)",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "STOP", "TAG": "STOP",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "STOP", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

STOP_CODONS = {"TAA", "TAG", "TGA"}
START_CODON = "ATG"


def validate_sequence(sequence: str) -> tuple[bool, str]:
    """
    Validate a DNA sequence string.

    Returns:
        (is_valid, error_message)
    """
    cleaned = sequence.upper().strip()
    if not cleaned:
        return False, "Sequence is empty."
    invalid = set(cleaned) - set("ATCGN")
    if invalid:
        return False, f"Invalid characters found: {', '.join(sorted(invalid))}"
    return True, ""


def clean_sequence(sequence: str) -> str:
    """Strip whitespace/newlines and uppercase."""
    return "".join(sequence.upper().split())


def gc_content(sequence: str) -> float:
    """Calculate GC content percentage."""
    seq = clean_sequence(sequence)
    if not seq:
        return 0.0
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100


def nucleotide_frequency(sequence: str) -> dict[str, int]:
    """Count occurrences of each nucleotide."""
    seq = clean_sequence(sequence)
    counts = Counter(seq)
    return {nuc: counts.get(nuc, 0) for nuc in "ATCGN"}


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    seq = clean_sequence(sequence)
    return seq.translate(COMPLEMENT_MAP)[::-1]


def find_motif(sequence: str, motif: str) -> list[int]:
    """
    Find all 1-based positions of a motif in the sequence.

    Returns:
        List of 1-based start positions.
    """
    seq = clean_sequence(sequence)
    motif = clean_sequence(motif)
    positions = []
    start = 0
    while True:
        pos = seq.find(motif, start)
        if pos == -1:
            break
        positions.append(pos + 1)  # 1-based
        start = pos + 1
    return positions


def find_orfs(sequence: str, min_length: int = 30) -> list[dict]:
    """
    Find all Open Reading Frames (ORFs) in all 6 reading frames.

    Args:
        sequence: DNA sequence string.
        min_length: Minimum ORF length in nucleotides (default 30).

    Returns:
        List of dicts with keys: frame, start, end, length, protein.
    """
    seq = clean_sequence(sequence)
    orfs = []

    strands = [
        (seq, "+"),
        (reverse_complement(seq), "-"),
    ]

    for strand_seq, strand in strands:
        for frame in range(3):
            i = frame
            while i < len(strand_seq) - 2:
                codon = strand_seq[i:i + 3]
                if codon == START_CODON:
                    # Found start — scan until STOP
                    protein_codons = []
                    j = i
                    while j < len(strand_seq) - 2:
                        c = strand_seq[j:j + 3]
                        if len(c) < 3:
                            break
                        if c in STOP_CODONS:
                            orf_len = j + 3 - i
                            if orf_len >= min_length:
                                if strand == "+":
                                    start_pos = i + 1
                                    end_pos = j + 3
                                else:
                                    end_pos = len(seq) - i
                                    start_pos = len(seq) - (j + 3) + 1
                                orfs.append({
                                    "frame": f"{strand}{frame + 1}",
                                    "start": start_pos,
                                    "end": end_pos,
                                    "length": orf_len,
                                    "protein": "-".join(
                                        CODON_TABLE.get(protein_codons[k], "???")
                                        for k in range(len(protein_codons))
                                    ),
                                })
                            break
                        protein_codons.append(c)
                        j += 3
                i += 1

    orfs.sort(key=lambda x: x["length"], reverse=True)
    return orfs


def melting_temperature(sequence: str) -> float:
    """
    Estimate melting temperature (Tm) of a DNA sequence.
    Uses the Wallace rule for short sequences (<14 bp),
    and the GC% formula for longer ones.
    """
    seq = clean_sequence(sequence)
    n = len(seq)
    a = seq.count("A")
    t = seq.count("T")
    g = seq.count("G")
    c = seq.count("C")

    if n < 14:
        return 2 * (a + t) + 4 * (g + c)
    else:
        gc = gc_content(seq)
        return 81.5 + 16.6 * (0) + 0.41 * gc - 675 / n


def dinucleotide_frequency(sequence: str) -> dict[str, int]:
    """Count all dinucleotide pairs in the sequence."""
    seq = clean_sequence(sequence)
    pairs = [seq[i:i + 2] for i in range(len(seq) - 1)]
    counts = Counter(pairs)
    nucleotides = "ATCG"
    return {
        f"{a}{b}": counts.get(f"{a}{b}", 0)
        for a in nucleotides for b in nucleotides
    }


def sequence_summary(sequence: str) -> dict:
    """Return a comprehensive summary of the DNA sequence."""
    seq = clean_sequence(sequence)
    freq = nucleotide_frequency(seq)
    return {
        "length": len(seq),
        "gc_content": round(gc_content(seq), 2),
        "at_content": round(100 - gc_content(seq), 2),
        "nucleotide_counts": freq,
        "melting_temp": round(melting_temperature(seq), 2),
        "reverse_complement": reverse_complement(seq),
        "num_orfs": len(find_orfs(seq)),
    }
