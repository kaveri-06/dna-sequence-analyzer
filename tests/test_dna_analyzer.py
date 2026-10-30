"""
Comprehensive Automated Tests for DNA Sequence Analyzer
=========================================================
Tests core bioinformatics algorithms, edge cases, formulas,
validation, FASTA parsing, and file exports.
"""

import json
import os
import tempfile
import unittest

from dna_analyzer import (
    analyze_dna,
    calculate_composition,
    calculate_gc_content,
    export_results_csv,
    export_results_json,
    find_motif,
    parse_fasta,
    sliding_window_gc,
    validate_sequence,
)


class TestDNAValidation(unittest.TestCase):
    """Tests for DNA sequence normalization and validation."""

    def test_valid_sequence(self):
        seq = "ATGCGATCGATC"
        self.assertEqual(validate_sequence(seq), "ATGCGATCGATC")

    def test_lowercase_conversion(self):
        seq = "atgcgatcgatc"
        self.assertEqual(validate_sequence(seq), "ATGCGATCGATC")

    def test_whitespace_and_newline_stripping(self):
        seq = "  ATGC \n GATC \t GATC \r\n "
        self.assertEqual(validate_sequence(seq), "ATGCGATCGATC")

    def test_ambiguous_n_allowed(self):
        seq = "ATGCNNNGATC"
        self.assertEqual(validate_sequence(seq), "ATGCNNNGATC")

    def test_empty_sequence_raises(self):
        with self.assertRaises(ValueError) as ctx:
            validate_sequence("")
        self.assertIn("cannot be empty", str(ctx.exception).lower())

    def test_whitespace_only_raises(self):
        with self.assertRaises(ValueError) as ctx:
            validate_sequence("   \n\t  ")
        self.assertIn("empty", str(ctx.exception).lower())

    def test_invalid_characters_raise(self):
        with self.assertRaises(ValueError) as ctx:
            validate_sequence("ATGCXGATC")
        self.assertIn("Invalid nucleotide", str(ctx.exception))
        self.assertIn("'X'", str(ctx.exception))

    def test_fasta_formatted_input_in_validation(self):
        raw_fasta = ">test_seq\nATGCGATC\nGATCGATC"
        self.assertEqual(validate_sequence(raw_fasta), "ATGCGATCGATCGATC")


class TestFASTAParser(unittest.TestCase):
    """Tests for FASTA header and sequence parsing."""

    def test_single_record_fasta(self):
        fasta = ">Homo_sapiens_gene\nATGC\nGATC\n"
        header, seq = parse_fasta(fasta)
        self.assertEqual(header, "Homo_sapiens_gene")
        self.assertEqual(seq, "ATGCGATC")

    def test_raw_sequence_without_header(self):
        raw = "ATGCGATC"
        header, seq = parse_fasta(raw)
        self.assertEqual(header, "")
        self.assertEqual(seq, "ATGCGATC")

    def test_fasta_with_spaces_and_numbers(self):
        fasta = ">seq1\n 10 ATGC 20 GATC\n"
        header, seq = parse_fasta(fasta)
        self.assertEqual(header, "seq1")
        self.assertEqual(seq, "ATGCGATC")


class TestCompositionCalculation(unittest.TestCase):
    """Tests for nucleotide base counting and percentage calculation."""

    def test_composition_counts(self):
        seq = "AAATTTCCGGN"
        counts, pcts, length = calculate_composition(seq)
        self.assertEqual(length, 11)
        self.assertEqual(counts["A"], 3)
        self.assertEqual(counts["T"], 3)
        self.assertEqual(counts["C"], 2)
        self.assertEqual(counts["G"], 2)
        self.assertEqual(counts["N"], 1)

    def test_composition_percentages(self):
        seq = "ATCG"
        counts, pcts, length = calculate_composition(seq)
        self.assertEqual(length, 4)
        for base in ["A", "T", "C", "G"]:
            self.assertEqual(pcts[base], 25.0)
        self.assertEqual(pcts["N"], 0.0)

    def test_empty_composition(self):
        counts, pcts, length = calculate_composition("")
        self.assertEqual(length, 0)
        self.assertEqual(counts["A"], 0)
        self.assertEqual(pcts["A"], 0.0)


class TestGCContentCalculation(unittest.TestCase):
    """Tests for GC content using canonical bases (A, T, C, G) denominator."""

    def test_balanced_gc(self):
        # 2 GC out of 4 canonical -> 50.0%
        self.assertEqual(calculate_gc_content("ATGC"), 50.0)

    def test_zero_gc(self):
        self.assertEqual(calculate_gc_content("AATT"), 0.0)

    def test_100_percent_gc(self):
        self.assertEqual(calculate_gc_content("GGCC"), 100.0)

    def test_canonical_denominator_with_ambiguous_bases(self):
        # 4 GC bases (GGCC) + 4 ambiguous bases (NNNN)
        # Canonical bases = 4, GC = 4. Old bug divided by 8 (giving 50%).
        # Correct biological requirement: denominator is canonical bases (4) -> 100.0%
        self.assertEqual(calculate_gc_content("GGCCNNNN"), 100.0)

    def test_no_canonical_bases_returns_none(self):
        # Only ambiguous bases
        self.assertIsNone(calculate_gc_content("NNNN"))

    def test_empty_sequence_returns_none(self):
        self.assertIsNone(calculate_gc_content(""))


class TestMotifFinding(unittest.TestCase):
    """Tests for motif detection, overlapping matches, and 1-based indexing."""

    def test_ecori_site_detection(self):
        seq = "ATGCGAATTCATGC"
        matches = find_motif(seq, "GAATTC")
        # 'GAATTC' starts at 0-based index 4 -> 1-based index 5
        self.assertEqual(matches, [5])

    def test_motif_at_start(self):
        seq = "GAATTCATGC"
        matches = find_motif(seq, "GAATTC")
        self.assertEqual(matches, [1])

    def test_overlapping_motifs(self):
        # 'AAA' in 'AAAA' occurs at 1-based index 1 and 2
        seq = "AAAA"
        matches = find_motif(seq, "AAA")
        self.assertEqual(matches, [1, 2])

    def test_multiple_occurrences(self):
        seq = "GAATTCGGGAATTC"
        matches = find_motif(seq, "GAATTC")
        self.assertEqual(matches, [1, 9])

    def test_motif_not_found(self):
        seq = "ATCGATCG"
        matches = find_motif(seq, "GAATTC")
        self.assertEqual(matches, [])

    def test_case_insensitive_motif(self):
        seq = "ATGCGAATTCATGC"
        matches = find_motif(seq, "gaattc")
        self.assertEqual(matches, [5])

    def test_invalid_motif_characters_raises(self):
        seq = "ATGCGAATTC"
        with self.assertRaises(ValueError) as ctx:
            find_motif(seq, "GAATTCN")  # 'N' is not canonical A, T, C, G
        self.assertIn("canonical nucleotides", str(ctx.exception))

    def test_empty_motif_raises(self):
        seq = "ATGCGAATTC"
        with self.assertRaises(ValueError):
            find_motif(seq, "")


class TestSlidingWindowGC(unittest.TestCase):
    """Tests for sliding-window GC calculations and short sequence safety."""

    def test_sliding_window_normal(self):
        seq = "GGCC" * 10  # 40 bp, all GC
        windows = sliding_window_gc(seq, window_size=20, step_size=5)
        self.assertGreater(len(windows), 0)
        for win in windows:
            self.assertEqual(win["gc_content"], 100.0)
            self.assertEqual(win["end"] - win["start"] + 1, 20)

    def test_sequence_shorter_than_window_does_not_crash(self):
        # Sequence of length 10 with window 20: 3 A, 3 T, 2 G, 2 C -> 4/10 = 40.0%
        seq = "ATGCATGCAT"
        windows = sliding_window_gc(seq, window_size=20)
        self.assertEqual(len(windows), 1)
        self.assertEqual(windows[0]["start"], 1)
        self.assertEqual(windows[0]["end"], 10)
        self.assertEqual(windows[0]["gc_content"], 40.0)

    def test_empty_sequence_sliding_window(self):
        self.assertEqual(sliding_window_gc("", window_size=20), [])

    def test_invalid_window_size_raises(self):
        with self.assertRaises(ValueError):
            sliding_window_gc("ATGC", window_size=0)
        with self.assertRaises(ValueError):
            sliding_window_gc("ATGC", window_size=10, step_size=0)


class TestAnalysisCoordinatorAndExport(unittest.TestCase):
    """Tests for analyze_dna() coordinator and CSV/JSON exports."""

    def test_full_analysis(self):
        seq = "ATGCGAATTC"
        res = analyze_dna(seq, motif="GAATTC", window_size=5)
        self.assertEqual(res["sequence_length"], 10)
        self.assertEqual(res["motif_match_count"], 1)
        self.assertEqual(res["motif_matches_1based"], [5])
        self.assertIsNotNone(res["gc_content"])
        self.assertIn("Balanced", res["gc_interpretation"])

    def test_csv_export(self):
        seq = "ATGCGAATTCATGC"
        res = analyze_dna(seq, motif="GAATTC")
        with tempfile.NamedTemporaryFile(suffix=".csv", delete=False) as f:
            temp_path = f.name
            
        try:
            export_results_csv(res, temp_path)
            self.assertTrue(os.path.exists(temp_path))
            with open(temp_path, "r", encoding="utf-8") as f:
                content = f.read()
            self.assertIn("DNA Sequence Analysis Report", content)
            self.assertIn("Total Sequence Length", content)
            self.assertIn("GAATTC", content)
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)

    def test_json_export(self):
        seq = "ATGCGAATTCATGC"
        res = analyze_dna(seq, motif="GAATTC")
        with tempfile.NamedTemporaryFile(suffix=".json", delete=False) as f:
            temp_path = f.name
            
        try:
            export_results_json(res, temp_path)
            self.assertTrue(os.path.exists(temp_path))
            with open(temp_path, "r", encoding="utf-8") as f:
                data = json.load(f)
            self.assertEqual(data["tool"], "DNA Sequence Analyzer")
            self.assertEqual(data["sequence_length_bp"], 14)
            self.assertEqual(data["motif_search"]["match_count"], 1)
            self.assertEqual(data["motif_search"]["positions_1based"], [5])
        finally:
            if os.path.exists(temp_path):
                os.remove(temp_path)


if __name__ == "__main__":
    unittest.main()
