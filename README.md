# 🧬 DNA Sequence Analyzer

A modern, fast, and beginner-friendly desktop application for analyzing DNA sequences, measuring nucleotide composition, detecting restriction motifs (like EcoRI `GAATTC`), and visualizing regional GC content distribution.

Designed specifically for **extreme simplicity**—a non-technical user can run the application with **a single double-click**, without knowing Python, typing commands, or installing software.

---

## 🚀 Quick Start (Running the Application)

### Option 1: Standalone Application (Recommended - No Python Needed)
Simply double-click:
```
DNA Sequence Analyzer.exe
```
* **Zero prerequisites**: Runs immediately on Windows without installing Python, Git, or libraries.
* **100% Local & Private**: All sequence computations run entirely on your computer. No data is ever sent to external servers or the cloud.

---

### Option 2: One-Click Launcher (`run.bat`)
If you prefer running from source or want to modify the Python code:
1. Double-click `run.bat`.
2. The launcher will automatically:
   * Detect Python on your computer.
   * Set up a clean, isolated virtual environment (`.venv`).
   * Automatically install required packages (`matplotlib`).
   * Open the graphical application interface automatically.

---

## ✨ Features

- **🛡️ Strict DNA Validation**: Handles raw sequences and FASTA formats. Accurately verifies valid nucleotides (`A`, `T`, `C`, `G`, `N`), strips formatting and line numbers, and presents friendly error messages for any invalid characters.
- **📊 Precise Nucleotide Breakdown**: Calculates exact counts and percentage frequencies for Adenine (A), Thymine (T), Cytosine (C), Guanine (G), and Ambiguous bases (N).
- **🔬 Canonical GC Content Ratio**: Uses the biologically standard canonical bases (`A`, `T`, `C`, `G`) as the denominator. If a sequence contains ambiguous bases (e.g., `N`), they do not distort the GC ratio. If no canonical bases exist, it safely displays `Unavailable`.
- **🎯 Restriction Site & Motif Search**:
  - Pre-configured to detect the widely used **EcoRI restriction enzyme site (`GAATTC`)**.
  - Supports searching for **any custom motif** (e.g., `HindIII`, `BamHI`, or custom oligomers).
  - Identifies **overlapping matches** (e.g., `AAA` in `AAAA`).
  - Reports positions using biological standard **1-based indexing**.
- **📈 Regional GC Content Profile**: Displays a sliding-window GC distribution curve using an embedded Matplotlib graph. Includes a reference line for sequence-wide average GC and markers for motif match locations. Handles sequences shorter than the selected window gracefully without crashing.
- **💾 One-Click Data & Graph Export**:
  - Save publication-ready high-resolution graphs as **PNG** images.
  - Export structured numerical reports as **CSV** spreadsheets.
  - Export complete machine-readable profiles as structured **JSON**.
- **🧪 Built-in Sample Sequences**: Includes a `⚡ Load Sample` button and sample files in the `sample_data/` folder to test features instantly.

---

## 🖥️ Using the Interface

1. **Provide a DNA Sequence**:
   - **Type or Paste**: Paste sequence text directly into the input box.
   - **Upload File**: Click `📂 Upload File` and choose a `.fasta`, `.fa`, `.fna`, or `.txt` file.
   - **Try Sample**: Click `⚡ Load Sample` to load a realistic plasmid sequence with an EcoRI site.
2. **Configure Parameters (Optional)**:
   - **Search Motif**: Default is `GAATTC` (EcoRI). Enter any sequence consisting of `A`, `T`, `C`, `G`.
   - **GC Window Size**: Adjust the sliding window length (default: `20 bp`).
3. **Analyze**:
   - Click the prominent **`🔬 Analyze DNA Sequence`** button.
4. **Review & Save**:
   - Inspect the KPI summary cards (Length, Overall GC %, Ambiguous Bases).
   - View the nucleotide frequency breakdown and motif positions.
   - Analyze the sliding-window GC curve.
   - Click `💾 Export CSV`, `💾 Export JSON`, or `📈 Save Graph (PNG)` to save your results.

---

## 📁 Project Structure

```
DNA-Analyzer/
├── DNA Sequence Analyzer.exe   # Standalone Windows executable (1-click launch)
├── run.bat                     # Automatic Python environment runner & launcher
├── build.bat                   # 1-click script to rebuild the standalone .exe
├── dna_analyzer.py             # Core bioinformatics calculation & analysis engine
├── gui.py                      # Modern scientific desktop interface (Tkinter + Matplotlib)
├── app.py                      # Application launch entry point
├── requirements.txt            # Python dependencies (matplotlib, pyinstaller)
├── .gitignore                  # Git exclusion rules
├── sample_data/                # Realistic sample sequences for testing
│   ├── ecori_plasmid_fragment.fasta
│   ├── sample_sequence.txt
│   └── ambiguous_sample.fasta
└── tests/                      # Automated unit test suite (35 tests)
    └── test_dna_analyzer.py
```

---

## 🛠️ Building the Standalone Executable (`build.bat`)

To rebuild `DNA Sequence Analyzer.exe` at any time:
1. Ensure Python 3.9+ is installed with `Add Python to PATH` enabled.
2. Double-click `build.bat`.
3. The script will automatically:
   - Verify Python and Pip.
   - Install build tools (`PyInstaller`, `matplotlib`).
   - Run the 35 automated unit tests to ensure correctness.
   - Compile and package everything into `DNA Sequence Analyzer.exe`.

---

## 🧪 Automated Testing

Automated tests cover nucleotide counting, GC calculations with canonical denominators, ambiguous base handling, motif lookaheads, sliding window boundaries, FASTA parsing, and file exports.

To run the tests manually:
```bash
python -m unittest discover tests
```
*Current test suite status: **35 passing tests** (100% pass rate).*

---

## 📚 Biological Background & Terminology

- **GC Content**: The percentage of nitrogenous bases in a DNA molecule that are either Guanine (G) or Cytosine (C). Because G-C pairs are linked by three hydrogen bonds (compared to two for A-T pairs), DNA with higher GC content is more thermally stable and has a higher melting temperature ($T_m$).
- **Canonical vs. Ambiguous Bases**:
  - Canonical bases: `A` (Adenine), `T` (Thymine), `C` (Cytosine), `G` (Guanine).
  - Ambiguous base `N`: Represents an unknown or uncalled base from sequencing equipment. The application excludes `N` from the GC denominator so that incomplete sequencing reads do not artificially skew GC stability calculations.
- **EcoRI Restriction Site (`GAATTC`)**: A palindromic 6-base recognition site cleaved by the restriction endonuclease enzyme *EcoRI* derived from *Escherichia coli*. Widely used in molecular cloning, plasmid construction, and restriction fragment length polymorphism (RFLP) assays.
