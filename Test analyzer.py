"""
Unit tests for dna_analyzer.analyzer module.
"""

import pytest
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


# ---------------------------------------------------------------------------
# validate_sequence
# ---------------------------------------------------------------------------

class TestValidateSequence:
    def test_valid_simple(self):
        ok, err = validate_sequence("ATCG")
        assert ok is True
        assert err == ""

    def test_valid_with_n(self):
        ok, err = validate_sequence("ATCGN")
        assert ok is True

    def test_valid_lowercase(self):
        ok, err = validate_sequence("atcg")
        assert ok is True

    def test_empty(self):
        ok, err = validate_sequence("")
        assert ok is False
        assert "empty" in err.lower()

    def test_invalid_chars(self):
        ok, err = validate_sequence("ATXCG")
        assert ok is False
        assert "X" in err

    def test_whitespace_only(self):
        ok, err = validate_sequence("   ")
        assert ok is False


# ---------------------------------------------------------------------------
# clean_sequence
# ---------------------------------------------------------------------------

class TestCleanSequence:
    def test_strips_whitespace(self):
        assert clean_sequence("  AT CG  ") == "ATCG"

    def test_uppercases(self):
        assert clean_sequence("atcg") == "ATCG"

    def test_removes_newlines(self):
        assert clean_sequence("AT\nCG\nTA") == "ATCGTA"


# ---------------------------------------------------------------------------
# gc_content
# ---------------------------------------------------------------------------

class TestGCContent:
    def test_all_gc(self):
        assert gc_content("GCGCGC") == pytest.approx(100.0)

    def test_all_at(self):
        assert gc_content("ATATAT") == pytest.approx(0.0)

    def test_mixed(self):
        assert gc_content("ATCG") == pytest.approx(50.0)

    def test_empty(self):
        assert gc_content("") == pytest.approx(0.0)

    def test_rounding(self):
        result = gc_content("ATCGATCG")
        assert 0 <= result <= 100


# ---------------------------------------------------------------------------
# nucleotide_frequency
# ---------------------------------------------------------------------------

class TestNucleotideFrequency:
    def test_counts(self):
        freq = nucleotide_frequency("AAATTTCCCGGG")
        assert freq["A"] == 3
        assert freq["T"] == 3
        assert freq["C"] == 3
        assert freq["G"] == 3
        assert freq["N"] == 0

    def test_with_n(self):
        freq = nucleotide_frequency("ATCGN")
        assert freq["N"] == 1

    def test_empty(self):
        freq = nucleotide_frequency("")
        assert all(v == 0 for v in freq.values())


# ---------------------------------------------------------------------------
# reverse_complement
# ---------------------------------------------------------------------------

class TestReverseComplement:
    def test_simple(self):
        assert reverse_complement("ATCG") == "CGAT"

    def test_palindrome(self):
        # AATT complement is TTAA, reversed is AATT
        assert reverse_complement("AATT") == "AATT"

    def test_self_complement(self):
        rc = reverse_complement("GCGC")
        assert rc == "GCGC"

    def test_double_reverse(self):
        seq = "ATCGATCG"
        assert reverse_complement(reverse_complement(seq)) == seq


# ---------------------------------------------------------------------------
# find_motif
# ---------------------------------------------------------------------------

class TestFindMotif:
    def test_single_match(self):
        positions = find_motif("ATCGATCG", "ATC")
        assert 1 in positions

    def test_multiple_matches(self):
        positions = find_motif("ATGATGATG", "ATG")
        assert len(positions) == 3
        assert positions == [1, 4, 7]

    def test_no_match(self):
        positions = find_motif("AAAA", "TTT")
        assert positions == []

    def test_overlapping(self):
        positions = find_motif("AAAA", "AA")
        assert len(positions) == 3

    def test_case_insensitive(self):
        positions = find_motif("atcg", "ATC")
        assert 1 in positions


# ---------------------------------------------------------------------------
# find_orfs
# ---------------------------------------------------------------------------

class TestFindOrfs:
    # ATG + several codons + stop codon
    SIMPLE_ORF = "AAATGCCCAAATAAGGG"  # ATG CCC AAA TAA

    def test_finds_orf(self):
        seq = "ATGCCCAAATAA"  # ATG CCC AAA TAA — 12 bp ORF
        orfs = find_orfs(seq, min_length=9)
        assert len(orfs) >= 1

    def test_min_length_filter(self):
        seq = "ATGCCCAAATAA"
        long_orfs = find_orfs(seq, min_length=100)
        assert long_orfs == []

    def test_returns_frame_info(self):
        seq = "ATGCCCAAATAAGGGG"
        orfs = find_orfs(seq, min_length=9)
        if orfs:
            assert "frame" in orfs[0]
            assert "start" in orfs[0]
            assert "end" in orfs[0]
            assert "length" in orfs[0]

    def test_sorted_by_length(self):
        seq = "ATGCCCAAATAAATGAAATAAATGCCCGGGAAATAA"
        orfs = find_orfs(seq, min_length=9)
        if len(orfs) > 1:
            lengths = [o["length"] for o in orfs]
            assert lengths == sorted(lengths, reverse=True)


# ---------------------------------------------------------------------------
# melting_temperature
# ---------------------------------------------------------------------------

class TestMeltingTemperature:
    def test_short_sequence(self):
        # Short seq: 2*(A+T) + 4*(G+C)
        # "AAAA" → 2*4 + 4*0 = 8
        tm = melting_temperature("AAAA")
        assert tm == pytest.approx(8.0)

    def test_long_sequence(self):
        tm = melting_temperature("ATCGATCGATCGATCG")  # 16 bp
        assert 50 <= tm <= 100  # reasonable range for a 16 bp sequence

    def test_gc_rich_higher_tm(self):
        gc_rich = melting_temperature("G" * 20)
        at_rich = melting_temperature("A" * 20)
        assert gc_rich > at_rich


# ---------------------------------------------------------------------------
# sequence_summary
# ---------------------------------------------------------------------------

class TestSequenceSummary:
    def test_returns_all_keys(self):
        summary = sequence_summary("ATCGATCGATCG")
        expected_keys = {"length", "gc_content", "at_content",
                         "nucleotide_counts", "melting_temp",
                         "reverse_complement", "num_orfs"}
        assert expected_keys.issubset(summary.keys())

    def test_gc_at_sum_to_100(self):
        summary = sequence_summary("ATCGATCG")
        assert summary["gc_content"] + summary["at_content"] == pytest.approx(100.0, abs=0.1)

    def test_length_correct(self):
        summary = sequence_summary("ATCG")
        assert summary["length"] == 4
