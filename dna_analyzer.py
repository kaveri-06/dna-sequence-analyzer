"""
Biological Sequence Analysis Tool
Author: B.Pharm / Bioinformatics Student
Description: A beginner-friendly computational biology script to process DNA sequences, 
calculate composition, find motifs, and analyze GC content distribution.
"""

import re
import matplotlib.pyplot as plt

def validate_sequence(dna_seq):
    """Ensures the DNA sequence contains only valid nucleotides (A, T, C, G)."""
    cleaned = dna_seq.upper().strip()
    cleaned = re.sub(r'\s+', '', cleaned)
    valid_nucleotides = set("ATCGN")
    if not all(n in valid_nucleotides for n in cleaned):
        raise ValueError("Invalid nucleotide found in sequence. Use only A, T, C, G, N.")
    return cleaned

def calculate_composition(dna_seq):
    """Calculates counts and percentages of each nucleotide."""
    length = len(dna_seq)
    composition = {
        'A': dna_seq.count('A'),
        'T': dna_seq.count('T'),
        'C': dna_seq.count('C'),
        'G': dna_seq.count('G'),
        'N': dna_seq.count('N')
    }
    percentages = {n: (count / length) * 100 for n, count in composition.items()}
    return composition, percentages, length

def calculate_gc_content(dna_seq):
    """Calculates overall GC content percentage (important for DNA stability/melting temp)."""
    g_count = dna_seq.count('G')
    c_count = dna_seq.count('C')
    gc_content = ((g_count + c_count) / len(dna_seq)) * 100
    return round(gc_content, 2)

def find_motif(dna_seq, motif):
    """Finds position indices of a specific pattern (e.g., restriction site GAATTC)."""
    motif = motif.upper()
    matches = [m.start() for m in re.finditer(f'(?={motif})', dna_seq)]
    return matches

def sliding_window_gc(dna_seq, window_size=30):
    """Calculates GC content across a sliding window for regional sequence analysis."""
    gc_values = []
    for i in range(len(dna_seq) - window_size + 1):
        window = dna_seq[i:i + window_size]
        g = window.count('G')
        c = window.count('C')
        gc_val = ((g + c) / window_size) * 100
        gc_values.append(gc_val)
    return gc_values

def main():
    print("=== DNA Sequence Analysis Tool ===")
    
    # Sample DNA sequence containing an EcoRI restriction site (GAATTC)
    sample_dna = (
        "ATGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        "GAATTCCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
        "TACGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG"
    )
    
    try:
        dna = validate_sequence(sample_dna)
        print(f"\n[INFO] DNA sequence loaded successfully.")
        
        # 1. Base Composition & Length
        counts, percentages, length = calculate_composition(dna)
        print(f"Sequence Length: {length} base pairs (bp)")
        print("\n--- Nucleotide Composition ---")
        for nuc in ['A', 'T', 'C', 'G']:
            print(f"  {nuc}: {counts[nuc]} ({percentages[nuc]:.2f}%)")
            
        # 2. GC Content Ratio
        gc = calculate_gc_content(dna)
        print(f"\nOverall GC Content: {gc}%")
        
        # 3. Restriction Site / Motif Search
        target_motif = "GAATTC"  # EcoRI Recognition Site
        matches = find_motif(dna, target_motif)
        print(f"\nSearching for EcoRI site ('{target_motif}'):")
        if matches:
            print(f"  Found site at index position: {matches}")
        else:
            print("  Site not found.")
            
        # 4. Graphical Output
        print("\n[INFO] Generating GC content plot...")
        window_size = 20
        gc_profile = sliding_window_gc(dna, window_size=window_size)
      
        plt.figure(figsize=(10, 4))
        plt.plot(range(len(gc_profile)), gc_profile, color='teal', linewidth=2)
        plt.title(f"Sliding Window GC Content Distribution (Window = {window_size} bp)")
        plt.xlabel("Sequence Position (bp)")
        plt.ylabel("GC Content (%)")
        plt.grid(True, linestyle='--', alpha=0.6)
        plt.ylim(0, 100)
        plt.tight_layout()
        plt.show()
        
    except ValueError as e:
        print(f"[ERROR] {e}")

if __name__ == "__main__":
    main()