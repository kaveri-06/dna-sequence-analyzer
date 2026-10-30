"""
DNA Sequence Analyzer - Core Bioinformatics Analysis Engine
===========================================================
A robust, beginner-friendly computational biology module for analyzing DNA
sequences, calculating nucleotide composition, computing GC content, detecting
sequence motifs (e.g. EcoRI restriction sites), and generating sliding-window
genomic metrics.

All calculations run 100% locally with zero external network requests.
"""

import csv
import json
import os
import re
from typing import Any, Dict, List, Optional, Tuple


def parse_fasta(content: str) -> Tuple[str, str]:
    """
    Parse a DNA sequence from raw text or FASTA formatted content.
    
    Extracts the FASTA header (if present) and concatenates all sequence
    lines, stripping headers, numbering, and whitespace.
    
    Args:
        content: Raw text or FASTA file content.
        
    Returns:
        Tuple of (header, sequence_string). If no FASTA header, header is empty.
    """
    lines = content.strip().splitlines()
    header_parts = []
    seq_parts = []
    
    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith(">"):
            header_parts.append(stripped[1:].strip())
        elif stripped.startswith(";"):
            # Semi-colon comment line in older FASTA formats
            continue
        else:
            # Remove line numbers or spacing sometimes found in GenBank / FASTA dumps
            cleaned_line = re.sub(r"[\s\d]+", "", stripped)
            seq_parts.append(cleaned_line)
            
    header = " | ".join(header_parts) if header_parts else ""
    sequence = "".join(seq_parts)
    return header, sequence


def validate_sequence(dna_seq: str) -> str:
    """
    Validate and normalize a DNA sequence.
    
    Strips whitespace and newlines, handles FASTA headers if passed directly,
    converts characters to uppercase, and ensures that only valid nucleotide
    characters (A, T, C, G, N) are present.
    
    Args:
        dna_seq: Input DNA sequence string.
        
    Returns:
        Cleaned, uppercase DNA sequence string containing only A, T, C, G, N.
        
    Raises:
        ValueError: If sequence is empty or contains invalid characters.
    """
    if not dna_seq or not isinstance(dna_seq, str):
        raise ValueError("DNA sequence cannot be empty. Please enter or paste a valid sequence.")
    
    # If a full FASTA string was passed, parse it
    if dna_seq.strip().startswith(">"):
        _, dna_seq = parse_fasta(dna_seq)
        
    # Remove all whitespace, line breaks, carriage returns, and tabs
    cleaned = re.sub(r"\s+", "", dna_seq).upper()
    
    if not cleaned:
        raise ValueError("DNA sequence is empty after removing whitespace. Please enter valid nucleotides.")
        
    valid_nucleotides = set("ATCGN")
    invalid_chars = sorted(list(set(c for c in cleaned if c not in valid_nucleotides)))
    
    if invalid_chars:
        sample_invalid = ", ".join(repr(c) for c in invalid_chars[:5])
        if len(invalid_chars) > 5:
            sample_invalid += "..."
        raise ValueError(
            f"Invalid nucleotide character(s) found: {sample_invalid}. "
            "DNA sequences must only contain A, T, C, G, or ambiguous base N."
        )
        
    return cleaned


def calculate_composition(dna_seq: str) -> Tuple[Dict[str, int], Dict[str, float], int]:
    """
    Calculate nucleotide counts and percentages for A, T, C, G, and N.
    
    Args:
        dna_seq: Normalized DNA sequence string.
        
    Returns:
        Tuple of (counts_dict, percentages_dict, sequence_length).
    """
    length = len(dna_seq)
    composition = {
        "A": dna_seq.count("A"),
        "T": dna_seq.count("T"),
        "C": dna_seq.count("C"),
        "G": dna_seq.count("G"),
        "N": dna_seq.count("N"),
    }
    
    if length == 0:
        percentages = {base: 0.0 for base in composition}
    else:
        percentages = {
            base: round((count / length) * 100, 2)
            for base, count in composition.items()
        }
        
    return composition, percentages, length


def calculate_gc_content(dna_seq: str) -> Optional[float]:
    """
    Calculate GC content percentage using canonical bases (A, T, C, G) as the denominator.
    
    If no canonical bases exist (e.g. sequence is empty or only consists of 'N' bases),
    the value is unavailable and returns None.
    
    Args:
        dna_seq: Normalized DNA sequence string.
        
    Returns:
        GC percentage rounded to 2 decimal places, or None if no canonical bases exist.
    """
    canonical_count = sum(dna_seq.count(base) for base in "ATCG")
    if canonical_count == 0:
        return None
        
    gc_count = dna_seq.count("G") + dna_seq.count("C")
    return round((gc_count / canonical_count) * 100, 2)


def find_motif(dna_seq: str, motif: str) -> List[int]:
    """
    Find all occurrence locations of a DNA motif, including overlapping matches.
    
    Positions are returned using standard biological 1-based indexing.
    
    Args:
        dna_seq: Normalized DNA sequence string.
        motif: Target nucleotide sequence to locate (e.g. EcoRI 'GAATTC').
        
    Returns:
        List of 1-based start positions where the motif occurs.
        
    Raises:
        ValueError: If motif is empty or contains characters other than A, T, C, G.
    """
    if not motif or not isinstance(motif, str):
        raise ValueError("Motif search pattern cannot be empty.")
        
    cleaned_motif = re.sub(r"\s+", "", motif).upper()
    if not cleaned_motif:
        raise ValueError("Motif search pattern cannot be blank.")
        
    if not all(base in "ATCG" for base in cleaned_motif):
        invalid = [b for b in cleaned_motif if b not in "ATCG"]
        raise ValueError(
            f"Motif must contain only canonical nucleotides (A, T, C, G). "
            f"Found invalid character(s): {', '.join(set(invalid))}."
        )
        
    # Use positive lookahead (?=...) to capture overlapping matches
    pattern = f"(?={re.escape(cleaned_motif)})"
    matches = [match.start() + 1 for match in re.finditer(pattern, dna_seq)]
    return matches


def sliding_window_gc(
    dna_seq: str, window_size: int = 20, step_size: int = 1
) -> List[Dict[str, Any]]:
    """
    Calculate GC content across fixed-size sliding sequence windows.
    
    Positions are 1-based. Uses canonical bases (A, T, C, G) as denominator
    for each window. If the sequence is shorter than window_size, handles
    gracefully without crashing by returning a single window spanning the
    available sequence.
    
    Args:
        dna_seq: Normalized DNA sequence string.
        window_size: Length of each window in base pairs (default 20).
        step_size: Step increment to advance the window (default 1).
        
    Returns:
        List of dicts:
        [{'start': int, 'end': int, 'midpoint': float, 'gc_content': Optional[float]}]
        
    Raises:
        ValueError: If window_size < 1 or step_size < 1.
    """
    if window_size < 1 or step_size < 1:
        raise ValueError("Window size and step size must be positive integers (>= 1).")
        
    seq_len = len(dna_seq)
    if seq_len == 0:
        return []
        
    # Graceful handling when sequence is shorter than requested window size
    if seq_len < window_size:
        canonical_count = sum(dna_seq.count(base) for base in "ATCG")
        gc = None
        if canonical_count > 0:
            gc = round(((dna_seq.count("G") + dna_seq.count("C")) / canonical_count) * 100, 2)
        return [{
            "start": 1,
            "end": seq_len,
            "midpoint": (1 + seq_len) / 2.0,
            "gc_content": gc,
            "note": "Window size exceeded sequence length; calculated across full sequence."
        }]
        
    results = []
    for start in range(0, seq_len - window_size + 1, step_size):
        window = dna_seq[start : start + window_size]
        canonical_count = sum(window.count(base) for base in "ATCG")
        gc = None
        if canonical_count > 0:
            gc_count = window.count("G") + window.count("C")
            gc = round((gc_count / canonical_count) * 100, 2)
            
        start_1based = start + 1
        end_1based = start + window_size
        midpoint = (start_1based + end_1based) / 2.0
        
        results.append({
            "start": start_1based,
            "end": end_1based,
            "midpoint": midpoint,
            "gc_content": gc,
        })
        
    return results


def analyze_dna(
    dna_seq: str,
    motif: str = "GAATTC",
    window_size: int = 20,
    step_size: int = 1,
    header: str = ""
) -> Dict[str, Any]:
    """
    Perform a complete biological sequence analysis on a DNA sequence.
    
    Args:
        dna_seq: Raw or pre-cleaned DNA sequence string.
        motif: Motif to search for (default EcoRI 'GAATTC').
        window_size: Sliding window size for GC profiling.
        step_size: Step size for sliding window.
        header: Optional sequence identifier or filename.
        
    Returns:
        Structured dictionary containing all analysis metrics.
    """
    cleaned_seq = validate_sequence(dna_seq)
    counts, percentages, length = calculate_composition(cleaned_seq)
    gc_content = calculate_gc_content(cleaned_seq)
    
    canonical_count = counts["A"] + counts["T"] + counts["C"] + counts["G"]
    ambiguous_count = counts["N"]
    
    # Motif search (if motif provided)
    motif_matches = []
    motif_error = None
    if motif:
        try:
            motif_matches = find_motif(cleaned_seq, motif)
        except ValueError as err:
            motif_error = str(err)
            
    # Sliding window GC profile
    gc_profile = sliding_window_gc(cleaned_seq, window_size=window_size, step_size=step_size)
    
    # Plain English summary explanation
    if gc_content is not None:
        if gc_content > 60:
            gc_interpretation = "High GC content (elevated thermal stability, strong G-C triple hydrogen bonding)"
        elif gc_content < 40:
            gc_interpretation = "Low GC content / AT-rich (lower thermal stability, easier strand separation)"
        else:
            gc_interpretation = "Moderate / Balanced GC content (typical for many standard genomic regions)"
    else:
        gc_interpretation = "Unavailable (sequence has no canonical A, T, C, G nucleotides)"
        
    return {
        "header": header,
        "sequence_length": length,
        "canonical_bases_count": canonical_count,
        "ambiguous_bases_count": ambiguous_count,
        "ambiguous_percentage": percentages.get("N", 0.0),
        "composition_counts": counts,
        "composition_percentages": percentages,
        "gc_content": gc_content,
        "gc_interpretation": gc_interpretation,
        "motif_searched": motif.upper() if motif else "",
        "motif_matches_1based": motif_matches,
        "motif_match_count": len(motif_matches),
        "motif_error": motif_error,
        "window_size": window_size,
        "step_size": step_size,
        "sliding_window_gc": gc_profile,
        "cleaned_sequence": cleaned_seq,
    }


def export_results_csv(analysis: Dict[str, Any], filepath: str) -> None:
    """
    Export sequence analysis summary and nucleotide composition to a CSV file.
    
    Args:
        analysis: Dictionary returned by analyze_dna().
        filepath: Destination file path.
    """
    with open(filepath, mode="w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["DNA Sequence Analysis Report"])
        writer.writerow(["Generated by", "DNA Sequence Analyzer"])
        if analysis.get("header"):
            writer.writerow(["Sequence Header", analysis["header"]])
        writer.writerow([])
        
        # Summary metrics
        writer.writerow(["Metric", "Value", "Unit / Note"])
        writer.writerow(["Total Sequence Length", analysis["sequence_length"], "base pairs (bp)"])
        writer.writerow(["Canonical Bases (A,T,C,G)", analysis["canonical_bases_count"], "bp"])
        writer.writerow(["Ambiguous Bases (N)", analysis["ambiguous_bases_count"], f"{analysis['ambiguous_percentage']}%"])
        writer.writerow([
            "Overall GC Content",
            f"{analysis['gc_content']}%" if analysis["gc_content"] is not None else "Unavailable",
            analysis["gc_interpretation"]
        ])
        writer.writerow([])
        
        # Composition breakdown
        writer.writerow(["Nucleotide", "Count (bp)", "Percentage (%)"])
        for base in ["A", "T", "C", "G", "N"]:
            writer.writerow([
                base,
                analysis["composition_counts"].get(base, 0),
                f"{analysis['composition_percentages'].get(base, 0.0):.2f}"
            ])
        writer.writerow([])
        
        # Motif search
        writer.writerow(["Motif Search Analysis"])
        writer.writerow(["Target Motif", analysis.get("motif_searched", "")])
        writer.writerow(["Total Matches", analysis.get("motif_match_count", 0)])
        positions_str = ", ".join(map(str, analysis.get("motif_matches_1based", [])))
        writer.writerow(["Positions (1-based)", positions_str if positions_str else "None"])
        writer.writerow([])
        
        # Sliding window GC profile (sample rows)
        writer.writerow(["Sliding Window GC Profile (Window Size: {} bp)".format(analysis.get("window_size", 20))])
        writer.writerow(["Start (1-based)", "End (1-based)", "Midpoint (bp)", "GC Content (%)"])
        for win in analysis.get("sliding_window_gc", []):
            gc_str = f"{win['gc_content']:.2f}" if win["gc_content"] is not None else "N/A"
            writer.writerow([win["start"], win["end"], win["midpoint"], gc_str])


def export_results_json(analysis: Dict[str, Any], filepath: str) -> None:
    """
    Export sequence analysis data to a structured JSON file.
    
    Args:
        analysis: Dictionary returned by analyze_dna().
        filepath: Destination file path.
    """
    # Create exportable copy
    export_data = {
        "tool": "DNA Sequence Analyzer",
        "header": analysis.get("header", ""),
        "sequence_length_bp": analysis.get("sequence_length", 0),
        "canonical_bases_count": analysis.get("canonical_bases_count", 0),
        "ambiguous_bases_count": analysis.get("ambiguous_bases_count", 0),
        "ambiguous_percentage": analysis.get("ambiguous_percentage", 0.0),
        "gc_content_percent": analysis.get("gc_content"),
        "gc_interpretation": analysis.get("gc_interpretation", ""),
        "composition_counts": analysis.get("composition_counts", {}),
        "composition_percentages": analysis.get("composition_percentages", {}),
        "motif_search": {
            "motif": analysis.get("motif_searched", ""),
            "match_count": analysis.get("motif_match_count", 0),
            "positions_1based": analysis.get("motif_matches_1based", []),
            "error": analysis.get("motif_error"),
        },
        "sliding_window": {
            "window_size_bp": analysis.get("window_size", 20),
            "step_size_bp": analysis.get("step_size", 1),
            "windows": analysis.get("sliding_window_gc", []),
        }
    }
    with open(filepath, mode="w", encoding="utf-8") as f:
        json.dump(export_data, f, indent=2)


def get_default_sample() -> str:
    """Returns a realistic sample DNA sequence containing an EcoRI site (GAATTC)."""
    return (
        "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        "GAATTCCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        "TACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    )


def main():
    """Application entry point: launches GUI by default or CLI if --cli is passed."""
    import sys
    
    if "--cli" in sys.argv:
        print("=" * 60)
        print("           DNA Sequence Analyzer (Console Mode)")
        print("=" * 60)
        sample_dna = get_default_sample()
        print("Running analysis on default sample sequence (contains EcoRI GAATTC site)...")
        
        try:
            results = analyze_dna(sample_dna, motif="GAATTC", window_size=20)
            print(f"\nSequence Length: {results['sequence_length']} bp")
            print(f"Overall GC Content: {results['gc_content']}% ({results['gc_interpretation']})")
            print("\nNucleotide Breakdown:")
            for base, count in results["composition_counts"].items():
                pct = results["composition_percentages"][base]
                print(f"  {base}: {count:4d} bp ({pct:5.2f}%)")
                
            print(f"\nMotif Search ('{results['motif_searched']}'):")
            if results["motif_matches_1based"]:
                print(f"  Matches found at 1-based position(s): {results['motif_matches_1based']}")
            else:
                print("  No matches found.")
                
            print(f"\nSliding window data calculated ({len(results['sliding_window_gc'])} windows).")
        except Exception as e:
            print(f"[ERROR] {e}")
    else:
        try:
            from gui import launch
            launch()
        except Exception as e:
            print(f"[NOTE] GUI could not be opened: {e}")
            print("Falling back to CLI mode...")
            sys.argv.append("--cli")
            main()


if __name__ == "__main__":
    main()