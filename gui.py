"""
DNA Sequence Analyzer - Desktop Graphical Interface
====================================================
A modern, scientific, beginner-friendly desktop application for analyzing DNA
sequences, calculating nucleotide composition, locating restriction motifs,
and visualizing regional GC content distribution.

Designed for high reliability, standalone execution, and clean usability.
"""

import os
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Any, Dict, Optional

# Enable High-DPI scaling on Windows displays for razor-sharp typography
if sys.platform.startswith("win"):
    try:
        from ctypes import windll
        windll.shcore.SetProcessDpiAwareness(1)
    except Exception:
        pass

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

from dna_analyzer import (
    analyze_dna,
    export_results_csv,
    export_results_json,
    get_default_sample,
    parse_fasta,
    validate_sequence,
)


# Color Palette (Scientific Modern Slate & Teal Theme)
BG_MAIN = "#f8fafc"         # Slate-50
BG_CARD = "#ffffff"         # Pure White
BG_HEADER = "#0f172a"       # Slate-900
TEXT_MAIN = "#0f172a"       # Slate-900
TEXT_MUTED = "#64748b"      # Slate-500
BORDER_COLOR = "#e2e8f0"    # Slate-200
ACCENT_TEAL = "#0d9488"     # Teal-600
ACCENT_TEAL_HOVER = "#0f766e"
ACCENT_BLUE = "#2563eb"     # Blue-600
ACCENT_BLUE_HOVER = "#1d4ed8"
BTN_SECONDARY_BG = "#f1f5f9"
BTN_SECONDARY_HOVER = "#e2e8f0"


class DNAAnalyzerGUI:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("DNA Sequence Analyzer")
        self.root.geometry("1180x820")
        self.root.minsize(980, 700)
        self.root.configure(bg=BG_MAIN)
        
        # State
        self.last_analysis: Optional[Dict[str, Any]] = None
        self.current_filename: str = ""
        
        # Style configuration
        self._setup_styles()
        
        # Build UI layout
        self._build_header()
        self._build_main_content()
        self._build_status_bar()
        
        # Load sample by default so the user immediately sees a working state
        self.load_sample_sequence()

    def _setup_styles(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
            
        style.configure("TFrame", background=BG_MAIN)
        style.configure("Card.TFrame", background=BG_CARD, relief="solid", borderwidth=1)
        style.configure("TLabel", background=BG_MAIN, foreground=TEXT_MAIN, font=("Segoe UI", 10))
        style.configure("Card.TLabel", background=BG_CARD, foreground=TEXT_MAIN, font=("Segoe UI", 10))
        style.configure("Muted.TLabel", background=BG_CARD, foreground=TEXT_MUTED, font=("Segoe UI", 9))
        style.configure("Heading.TLabel", background=BG_CARD, foreground=TEXT_MAIN, font=("Segoe UI", 11, "bold"))
        
        # Custom button styles
        style.configure("Primary.TButton",
                        font=("Segoe UI", 11, "bold"),
                        background=ACCENT_TEAL,
                        foreground="#ffffff",
                        padding=8)
        style.map("Primary.TButton",
                  background=[("active", ACCENT_TEAL_HOVER), ("pressed", "#115e59")])
                  
        style.configure("Action.TButton",
                        font=("Segoe UI", 9, "bold"),
                        background=ACCENT_BLUE,
                        foreground="#ffffff",
                        padding=5)
        style.map("Action.TButton",
                  background=[("active", ACCENT_BLUE_HOVER), ("pressed", "#1e40af")])

        style.configure("Secondary.TButton",
                        font=("Segoe UI", 9),
                        background=BTN_SECONDARY_BG,
                        foreground=TEXT_MAIN,
                        padding=5)
        style.map("Secondary.TButton",
                  background=[("active", BTN_SECONDARY_HOVER)])

    def _build_header(self):
        header_frame = tk.Frame(self.root, bg=BG_HEADER, height=75)
        header_frame.pack(fill=tk.X, side=tk.TOP)
        header_frame.pack_propagate(False)
        
        content = tk.Frame(header_frame, bg=BG_HEADER)
        content.pack(fill=tk.BOTH, expand=True, padx=24, pady=12)
        
        title_label = tk.Label(
            content,
            text="🧬 DNA Sequence Analyzer",
            font=("Segoe UI", 18, "bold"),
            bg=BG_HEADER,
            fg="#ffffff"
        )
        title_label.pack(side=tk.LEFT, anchor=tk.W)
        
        subtitle_label = tk.Label(
            content,
            text="Local Genomic Metrics • Motif Detection • GC Sliding Window",
            font=("Segoe UI", 9),
            bg=BG_HEADER,
            fg="#94a3b8"
        )
        subtitle_label.pack(side=tk.LEFT, anchor=tk.W, padx=16, pady=(4, 0))
        
        # Top right sample quick button
        top_btn_frame = tk.Frame(content, bg=BG_HEADER)
        top_btn_frame.pack(side=tk.RIGHT)
        
        load_sample_top = tk.Button(
            top_btn_frame,
            text="⚡ Load Sample",
            font=("Segoe UI", 9, "bold"),
            bg="#334155",
            fg="#ffffff",
            activebackground="#475569",
            activeforeground="#ffffff",
            relief="flat",
            padx=12,
            pady=4,
            cursor="hand2",
            command=self.load_sample_sequence
        )
        load_sample_top.pack(side=tk.RIGHT, padx=4)

    def _build_main_content(self):
        # Two-column layout: Left (Input & Controls), Right (Results & Visualization)
        main_container = tk.Frame(self.root, bg=BG_MAIN)
        main_container.pack(fill=tk.BOTH, expand=True, padx=20, pady=12)
        
        # Left Panel (Width ~420)
        left_panel = tk.Frame(main_container, bg=BG_MAIN, width=440)
        left_panel.pack(side=tk.LEFT, fill=tk.BOTH, padx=(0, 10))
        left_panel.pack_propagate(False)
        
        self._build_input_card(left_panel)
        self._build_motif_card(left_panel)
        
        # Right Panel (Results & Plot)
        right_panel = tk.Frame(main_container, bg=BG_MAIN)
        right_panel.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=(10, 0))
        
        self._build_metrics_card(right_panel)
        self._build_plot_card(right_panel)

    def _build_input_card(self, parent: tk.Frame):
        card = tk.Frame(parent, bg=BG_CARD, highlightbackground=BORDER_COLOR, highlightthickness=1)
        card.pack(fill=tk.BOTH, expand=True, pady=(0, 10))
        
        header = tk.Frame(card, bg=BG_CARD)
        header.pack(fill=tk.X, padx=14, pady=(12, 6))
        
        tk.Label(
            header,
            text="DNA Sequence Input",
            font=("Segoe UI", 12, "bold"),
            bg=BG_CARD,
            fg=TEXT_MAIN
        ).pack(side=tk.LEFT)
        
        # File buttons
        btn_bar = tk.Frame(header, bg=BG_CARD)
        btn_bar.pack(side=tk.RIGHT)
        
        upload_btn = tk.Button(
            btn_bar,
            text="📂 Upload File",
            font=("Segoe UI", 9),
            bg=BTN_SECONDARY_BG,
            fg=TEXT_MAIN,
            activebackground=BTN_SECONDARY_HOVER,
            relief="flat",
            padx=8,
            pady=2,
            cursor="hand2",
            command=self.upload_file
        )
        upload_btn.pack(side=tk.LEFT, padx=3)
        
        clear_btn = tk.Button(
            btn_bar,
            text="Clear",
            font=("Segoe UI", 9),
            bg=BTN_SECONDARY_BG,
            fg=TEXT_MAIN,
            activebackground=BTN_SECONDARY_HOVER,
            relief="flat",
            padx=8,
            pady=2,
            cursor="hand2",
            command=self.clear_input
        )
        clear_btn.pack(side=tk.LEFT, padx=3)
        
        # Sequence Text Area
        txt_frame = tk.Frame(card, bg=BG_CARD)
        txt_frame.pack(fill=tk.BOTH, expand=True, padx=14, pady=6)
        
        scroll = tk.Scrollbar(txt_frame)
        scroll.pack(side=tk.RIGHT, fill=tk.Y)
        
        self.seq_text = tk.Text(
            txt_frame,
            wrap=tk.CHAR,
            font=("Consolas", 10),
            bg="#fdfefe",
            fg="#1e293b",
            relief="solid",
            borderwidth=1,
            highlightthickness=0,
            yscrollcommand=scroll.set,
            height=12
        )
        self.seq_text.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scroll.config(command=self.seq_text.yview)
        
        # Sequence info bar under text box
        self.seq_info_label = tk.Label(
            card,
            text="Accepted bases: A, T, C, G, N. Formats: Raw sequence, FASTA (.fasta, .fa, .fna, .txt)",
            font=("Segoe UI", 8),
            bg=BG_CARD,
            fg=TEXT_MUTED,
            anchor="w",
            justify=tk.LEFT
        )
        self.seq_info_label.pack(fill=tk.X, padx=14, pady=(2, 10))

    def _build_motif_card(self, parent: tk.Frame):
        card = tk.Frame(parent, bg=BG_CARD, highlightbackground=BORDER_COLOR, highlightthickness=1)
        card.pack(fill=tk.X, pady=(0, 0))
        
        inner = tk.Frame(card, bg=BG_CARD)
        inner.pack(fill=tk.X, padx=14, pady=12)
        
        tk.Label(
            inner,
            text="Motif & Window Parameters",
            font=("Segoe UI", 11, "bold"),
            bg=BG_CARD,
            fg=TEXT_MAIN
        ).pack(anchor="w", pady=(0, 8))
        
        # Motif Search Entry
        motif_row = tk.Frame(inner, bg=BG_CARD)
        motif_row.pack(fill=tk.X, pady=4)
        
        tk.Label(
            motif_row,
            text="Search Motif:",
            font=("Segoe UI", 9, "bold"),
            bg=BG_CARD,
            fg=TEXT_MAIN,
            width=14,
            anchor="w"
        ).pack(side=tk.LEFT)
        
        self.motif_var = tk.StringVar(value="GAATTC")
        self.motif_entry = tk.Entry(
            motif_row,
            textvariable=self.motif_var,
            font=("Consolas", 10, "bold"),
            relief="solid",
            borderwidth=1,
            width=14
        )
        self.motif_entry.pack(side=tk.LEFT, padx=(0, 8))
        
        tk.Label(
            motif_row,
            text="(Default: EcoRI 'GAATTC')",
            font=("Segoe UI", 8),
            bg=BG_CARD,
            fg=TEXT_MUTED
        ).pack(side=tk.LEFT)
        
        # Window size row
        win_row = tk.Frame(inner, bg=BG_CARD)
        win_row.pack(fill=tk.X, pady=4)
        
        tk.Label(
            win_row,
            text="GC Window Size:",
            font=("Segoe UI", 9, "bold"),
            bg=BG_CARD,
            fg=TEXT_MAIN,
            width=14,
            anchor="w"
        ).pack(side=tk.LEFT)
        
        self.window_var = tk.IntVar(value=20)
        self.window_spin = tk.Spinbox(
            win_row,
            from_=5,
            to=500,
            increment=5,
            textvariable=self.window_var,
            font=("Segoe UI", 9),
            relief="solid",
            borderwidth=1,
            width=6
        )
        self.window_spin.pack(side=tk.LEFT, padx=(0, 8))
        
        tk.Label(
            win_row,
            text="bp (sliding window resolution)",
            font=("Segoe UI", 8),
            bg=BG_CARD,
            fg=TEXT_MUTED
        ).pack(side=tk.LEFT)
        
        # Prominent Analyze DNA Button
        self.analyze_btn = tk.Button(
            inner,
            text="🔬 Analyze DNA Sequence",
            font=("Segoe UI", 12, "bold"),
            bg=ACCENT_TEAL,
            fg="#ffffff",
            activebackground=ACCENT_TEAL_HOVER,
            activeforeground="#ffffff",
            relief="flat",
            padx=16,
            pady=8,
            cursor="hand2",
            command=self.run_analysis
        )
        self.analyze_btn.pack(fill=tk.X, pady=(12, 4))

    def _build_metrics_card(self, parent: tk.Frame):
        card = tk.Frame(parent, bg=BG_CARD, highlightbackground=BORDER_COLOR, highlightthickness=1)
        card.pack(fill=tk.X, pady=(0, 10))
        
        inner = tk.Frame(card, bg=BG_CARD)
        inner.pack(fill=tk.X, padx=14, pady=12)
        
        # Header with export buttons
        top_row = tk.Frame(inner, bg=BG_CARD)
        top_row.pack(fill=tk.X, pady=(0, 8))
        
        tk.Label(
            top_row,
            text="Sequence Metrics & Composition",
            font=("Segoe UI", 12, "bold"),
            bg=BG_CARD,
            fg=TEXT_MAIN
        ).pack(side=tk.LEFT)
        
        export_box = tk.Frame(top_row, bg=BG_CARD)
        export_box.pack(side=tk.RIGHT)
        
        self.csv_btn = tk.Button(
            export_box,
            text="💾 Export CSV",
            font=("Segoe UI", 9),
            bg=BTN_SECONDARY_BG,
            fg=TEXT_MAIN,
            activebackground=BTN_SECONDARY_HOVER,
            relief="flat",
            padx=8,
            pady=2,
            cursor="hand2",
            command=self.export_csv
        )
        self.csv_btn.pack(side=tk.LEFT, padx=3)
        
        self.json_btn = tk.Button(
            export_box,
            text="💾 Export JSON",
            font=("Segoe UI", 9),
            bg=BTN_SECONDARY_BG,
            fg=TEXT_MAIN,
            activebackground=BTN_SECONDARY_HOVER,
            relief="flat",
            padx=8,
            pady=2,
            cursor="hand2",
            command=self.export_json
        )
        self.json_btn.pack(side=tk.LEFT, padx=3)
        
        # Quick Stat KPI tiles
        kpi_frame = tk.Frame(inner, bg=BG_CARD)
        kpi_frame.pack(fill=tk.X, pady=6)
        
        # 3 KPI boxes: Length, GC Content, Ambiguous Bases
        self.kpi_len = self._create_kpi_box(kpi_frame, "Sequence Length", "-- bp", "Total Nucleotides")
        self.kpi_len.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(0, 6))
        
        self.kpi_gc = self._create_kpi_box(kpi_frame, "Overall GC Content", "-- %", "Canonical Bases (A,T,C,G)")
        self.kpi_gc.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=3)
        
        self.kpi_amb = self._create_kpi_box(kpi_frame, "Ambiguous Bases (N)", "-- bp", "Uncalled Bases")
        self.kpi_amb.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(6, 0))
        
        # Nucleotide Composition Table
        comp_frame = tk.Frame(inner, bg=BG_CARD)
        comp_frame.pack(fill=tk.X, pady=(10, 4))
        
        tk.Label(
            comp_frame,
            text="Nucleotide Counts & Frequencies:",
            font=("Segoe UI", 9, "bold"),
            bg=BG_CARD,
            fg=TEXT_MAIN
        ).pack(anchor="w", pady=(0, 4))
        
        self.base_labels: Dict[str, tk.Label] = {}
        bases_row = tk.Frame(comp_frame, bg=BG_CARD)
        bases_row.pack(fill=tk.X)
        
        # Colors for bases: A, T, C, G, N
        colors = {
            "A": "#2563eb",  # Blue
            "T": "#dc2626",  # Red
            "C": "#d97706",  # Amber/Orange
            "G": "#059669",  # Green
            "N": "#64748b",  # Grey
        }
        
        for base in ["A", "T", "C", "G", "N"]:
            pill = tk.Frame(bases_row, bg="#f8fafc", highlightbackground=BORDER_COLOR, highlightthickness=1)
            pill.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=2)
            
            top_bar = tk.Frame(pill, bg=colors[base], height=3)
            top_bar.pack(fill=tk.X)
            
            b_name = "Ambiguous (N)" if base == "N" else f"Base {base}"
            tk.Label(pill, text=b_name, font=("Segoe UI", 8, "bold"), bg="#f8fafc", fg=colors[base]).pack(pady=(3, 0))
            
            lbl = tk.Label(pill, text="0 (0.00%)", font=("Segoe UI", 9), bg="#f8fafc", fg=TEXT_MAIN)
            lbl.pack(pady=(0, 4))
            self.base_labels[base] = lbl
            
        # Motif Result Summary Row
        self.motif_result_frame = tk.Frame(inner, bg="#f1f5f9", highlightbackground=BORDER_COLOR, highlightthickness=1)
        self.motif_result_frame.pack(fill=tk.X, pady=(8, 0))
        
        self.motif_result_label = tk.Label(
            self.motif_result_frame,
            text="Motif search: Ready. Click 'Analyze DNA Sequence' to search.",
            font=("Segoe UI", 9),
            bg="#f1f5f9",
            fg=TEXT_MAIN,
            anchor="w",
            padx=10,
            pady=6,
            justify=tk.LEFT
        )
        self.motif_result_label.pack(fill=tk.X)

    def _create_kpi_box(self, parent: tk.Frame, title: str, value: str, subtitle: str) -> tk.Frame:
        box = tk.Frame(parent, bg="#f8fafc", highlightbackground=BORDER_COLOR, highlightthickness=1)
        
        tk.Label(box, text=title, font=("Segoe UI", 9), bg="#f8fafc", fg=TEXT_MUTED).pack(anchor="w", padx=10, pady=(6, 0))
        val_lbl = tk.Label(box, text=value, font=("Segoe UI", 14, "bold"), bg="#f8fafc", fg=TEXT_MAIN)
        val_lbl.pack(anchor="w", padx=10, pady=(2, 0))
        sub_lbl = tk.Label(box, text=subtitle, font=("Segoe UI", 8), bg="#f8fafc", fg=TEXT_MUTED)
        sub_lbl.pack(anchor="w", padx=10, pady=(0, 6))
        
        # Store label references
        box.val_label = val_lbl
        box.sub_label = sub_lbl
        return box

    def _build_plot_card(self, parent: tk.Frame):
        card = tk.Frame(parent, bg=BG_CARD, highlightbackground=BORDER_COLOR, highlightthickness=1)
        card.pack(fill=tk.BOTH, expand=True)
        
        header = tk.Frame(card, bg=BG_CARD)
        header.pack(fill=tk.X, padx=14, pady=(10, 4))
        
        tk.Label(
            header,
            text="GC Content Distribution (Sliding Window)",
            font=("Segoe UI", 11, "bold"),
            bg=BG_CARD,
            fg=TEXT_MAIN
        ).pack(side=tk.LEFT)
        
        save_plot_btn = tk.Button(
            header,
            text="📈 Save Graph (PNG)",
            font=("Segoe UI", 9),
            bg=BTN_SECONDARY_BG,
            fg=TEXT_MAIN,
            activebackground=BTN_SECONDARY_HOVER,
            relief="flat",
            padx=8,
            pady=2,
            cursor="hand2",
            command=self.save_plot_image
        )
        save_plot_btn.pack(side=tk.RIGHT)
        
        # Matplotlib embedded canvas
        self.fig, self.ax = plt.subplots(figsize=(7.0, 3.2), dpi=100)
        self.fig.patch.set_facecolor("#ffffff")
        self._init_empty_plot()
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=card)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True, padx=14, pady=(4, 10))

    def _init_empty_plot(self):
        self.ax.clear()
        self.ax.set_facecolor("#fafafa")
        self.ax.text(
            0.5, 0.5,
            "Enter sequence and click 'Analyze DNA Sequence' to plot GC distribution",
            ha="center", va="center", transform=self.ax.transAxes,
            color=TEXT_MUTED, fontsize=10, style="italic"
        )
        self.ax.set_xlim(0, 100)
        self.ax.set_ylim(0, 100)
        self.ax.set_xlabel("Sequence Position (bp)", fontsize=9, color=TEXT_MUTED)
        self.ax.set_ylabel("GC Content (%)", fontsize=9, color=TEXT_MUTED)
        self.ax.tick_params(colors=TEXT_MUTED, labelsize=8)
        for spine in self.ax.spines.values():
            spine.set_color(BORDER_COLOR)
        self.fig.tight_layout()

    def _build_status_bar(self):
        self.status_bar = tk.Frame(self.root, bg="#f1f5f9", height=28, highlightbackground=BORDER_COLOR, highlightthickness=1)
        self.status_bar.pack(fill=tk.X, side=tk.BOTTOM)
        self.status_bar.pack_propagate(False)
        
        self.status_label = tk.Label(
            self.status_bar,
            text="Ready • Standalone local analysis mode",
            font=("Segoe UI", 9),
            bg="#f1f5f9",
            fg=TEXT_MUTED,
            anchor="w"
        )
        self.status_label.pack(side=tk.LEFT, padx=16, fill=tk.X, expand=True)

    # ---------------- User Interactions ---------------- #

    def set_status(self, text: str, is_error: bool = False):
        color = "#dc2626" if is_error else TEXT_MUTED
        self.status_label.config(text=text, fg=color)

    def load_sample_sequence(self):
        sample = get_default_sample()
        self.seq_text.delete("1.0", tk.END)
        self.seq_text.insert(tk.END, sample)
        self.motif_var.set("GAATTC")
        self.window_var.set(20)
        self.current_filename = "sample_pBR322_fragment.dna"
        self.seq_info_label.config(
            text=f"Loaded realistic sample sequence ({len(sample)} bp, contains EcoRI 'GAATTC' site)."
        )
        self.run_analysis()

    def upload_file(self):
        filetypes = [
            ("All Supported Formats", "*.fasta;*.fa;*.fna;*.txt"),
            ("FASTA Files", "*.fasta;*.fa;*.fna"),
            ("Text Files", "*.txt"),
            ("All Files", "*.*")
        ]
        path = filedialog.askopenfilename(title="Select DNA Sequence File", filetypes=filetypes)
        if not path:
            return
            
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                content = f.read()
                
            header, sequence = parse_fasta(content)
            if not sequence and not content.strip():
                messagebox.showwarning("Empty File", "The selected file is empty. Please choose a valid DNA sequence file.")
                return
                
            self.seq_text.delete("1.0", tk.END)
            self.seq_text.insert(tk.END, sequence if sequence else content)
            self.current_filename = os.path.basename(path)
            
            header_display = f" | Header: {header}" if header else ""
            self.seq_info_label.config(
                text=f"Loaded file: {self.current_filename}{header_display}"
            )
            self.set_status(f"File loaded successfully: {self.current_filename}")
            
            # Automatically analyze upon loading
            self.run_analysis()
            
        except Exception as err:
            messagebox.showerror(
                "Error Reading File",
                f"Could not read the selected file:\n{err}\n\nPlease verify that the file is an accessible text or FASTA file."
            )
            self.set_status(f"Error opening file: {err}", is_error=True)

    def clear_input(self):
        self.seq_text.delete("1.0", tk.END)
        self.current_filename = ""
        self.seq_info_label.config(
            text="Accepted bases: A, T, C, G, N. Formats: Raw sequence, FASTA (.fasta, .fa, .fna, .txt)"
        )
        self.set_status("Sequence cleared.")

    def run_analysis(self):
        raw_text = self.seq_text.get("1.0", tk.END).strip()
        if not raw_text:
            messagebox.showwarning(
                "Input Needed",
                "Please enter or paste a DNA sequence, or click '⚡ Load Sample' to try the application."
            )
            self.set_status("Analysis halted: No sequence provided.", is_error=True)
            return
            
        motif_query = self.motif_var.get().strip()
        try:
            win_size = int(self.window_var.get())
            if win_size < 1:
                win_size = 20
                self.window_var.set(20)
        except ValueError:
            win_size = 20
            self.window_var.set(20)
            
        try:
            analysis = analyze_dna(
                raw_text,
                motif=motif_query,
                window_size=win_size,
                header=self.current_filename
            )
        except ValueError as val_err:
            messagebox.showerror("Sequence Validation Error", str(val_err))
            self.set_status(f"Validation Error: {val_err}", is_error=True)
            return
        except Exception as err:
            messagebox.showerror("Analysis Error", f"An unexpected error occurred during analysis:\n{err}")
            self.set_status(f"Analysis error: {err}", is_error=True)
            return
            
        self.last_analysis = analysis
        self._update_display(analysis)
        self.set_status(
            f"Analysis complete: {analysis['sequence_length']} bp analyzed. "
            f"GC: {analysis['gc_content']}% | Motifs found: {analysis['motif_match_count']}"
        )

    def _update_display(self, data: Dict[str, Any]):
        # 1. Update KPI tiles
        length = data["sequence_length"]
        self.kpi_len.val_label.config(text=f"{length:,} bp")
        self.kpi_len.sub_label.config(text=f"{data['canonical_bases_count']:,} canonical bases")
        
        gc = data["gc_content"]
        if gc is not None:
            self.kpi_gc.val_label.config(text=f"{gc:.2f}%")
            self.kpi_gc.sub_label.config(text=data["gc_interpretation"][:32] + "...")
        else:
            self.kpi_gc.val_label.config(text="Unavailable")
            self.kpi_gc.sub_label.config(text="No canonical bases (A, T, C, G)")
            
        amb_count = data["ambiguous_bases_count"]
        amb_pct = data["ambiguous_percentage"]
        self.kpi_amb.val_label.config(text=f"{amb_count:,} bp")
        self.kpi_amb.sub_label.config(text=f"{amb_pct:.2f}% of sequence")
        
        # 2. Update Composition Table
        for base in ["A", "T", "C", "G", "N"]:
            count = data["composition_counts"].get(base, 0)
            pct = data["composition_percentages"].get(base, 0.0)
            self.base_labels[base].config(text=f"{count:,} ({pct:.2f}%)")
            
        # 3. Update Motif Search Results
        motif = data.get("motif_searched", "")
        matches = data.get("motif_matches_1based", [])
        motif_err = data.get("motif_error")
        
        if motif_err:
            self.motif_result_label.config(
                text=f"⚠️ Motif search error for '{motif}': {motif_err}",
                fg="#dc2626"
            )
        elif not motif:
            self.motif_result_label.config(
                text="Motif search: No motif specified.",
                fg=TEXT_MUTED
            )
        elif len(matches) == 0:
            self.motif_result_label.config(
                text=f"Motif '{motif}': No occurrences found in this sequence.",
                fg=TEXT_MAIN
            )
        else:
            # Format 1-based positions cleanly
            pos_preview = ", ".join(map(str, matches[:10]))
            if len(matches) > 10:
                pos_preview += f" ... (+{len(matches) - 10} more)"
            self.motif_result_label.config(
                text=f"🎯 Motif '{motif}' found {len(matches)} time(s) at 1-based position(s): {pos_preview}",
                fg="#047857"
            )
            
        # 4. Render GC Sliding Window Plot
        self._render_plot(data)

    def _render_plot(self, data: Dict[str, Any]):
        self.ax.clear()
        self.ax.set_facecolor("#ffffff")
        
        windows = data.get("sliding_window_gc", [])
        if not windows:
            self._init_empty_plot()
            self.canvas.draw()
            return
            
        x_positions = [w["midpoint"] for w in windows]
        y_gc = [w["gc_content"] if w["gc_content"] is not None else float("nan") for w in windows]
        
        # Plot styled line and filled area
        line = self.ax.plot(
            x_positions,
            y_gc,
            color=ACCENT_TEAL,
            linewidth=2.0,
            label=f"Sliding GC (Window: {data['window_size']} bp)"
        )
        self.ax.fill_between(x_positions, y_gc, color=ACCENT_TEAL, alpha=0.12)
        
        # Horizontal reference line for mean GC content
        mean_gc = data.get("gc_content")
        if mean_gc is not None:
            self.ax.axhline(
                mean_gc,
                color="#f59e0b",
                linestyle="--",
                linewidth=1.4,
                label=f"Average GC ({mean_gc:.2f}%)"
            )
            
        # Highlight motif matches on the plot as vertical markers
        matches = data.get("motif_matches_1based", [])
        if matches and len(matches) <= 25:
            for idx, pos in enumerate(matches):
                lbl = "Motif Match" if idx == 0 else ""
                self.ax.axvline(pos, color="#ef4444", linestyle=":", linewidth=1.2, alpha=0.7, label=lbl)
                
        self.ax.set_title(
            f"Sliding-Window GC Content Profile (Window Size = {data['window_size']} bp)",
            fontsize=10,
            fontweight="bold",
            color=TEXT_MAIN,
            pad=8
        )
        self.ax.set_xlabel("Sequence Position (bp)", fontsize=9, color=TEXT_MUTED)
        self.ax.set_ylabel("GC Content (%)", fontsize=9, color=TEXT_MUTED)
        self.ax.set_ylim(0, 100)
        self.ax.grid(True, linestyle="--", alpha=0.5, color=BORDER_COLOR)
        self.ax.legend(loc="upper right", fontsize=8, framealpha=0.85)
        
        # Polish axis spines
        for spine in self.ax.spines.values():
            spine.set_color(BORDER_COLOR)
            spine.set_linewidth(1.0)
            
        self.fig.tight_layout()
        self.canvas.draw()

    def export_csv(self):
        if not self.last_analysis:
            messagebox.showinfo("No Analysis", "Please run sequence analysis before exporting.")
            return
            
        default_name = f"dna_analysis_{self.last_analysis['sequence_length']}bp.csv"
        path = filedialog.asksaveasfilename(
            title="Export Analysis as CSV",
            defaultextension=".csv",
            initialfile=default_name,
            filetypes=[("CSV Files", "*.csv"), ("All Files", "*.*")]
        )
        if not path:
            return
            
        try:
            export_results_csv(self.last_analysis, path)
            messagebox.showinfo("Export Successful", f"Analysis results saved successfully to:\n{path}")
            self.set_status(f"Exported CSV: {os.path.basename(path)}")
        except Exception as err:
            messagebox.showerror("Export Failed", f"Could not save CSV file:\n{err}")
            self.set_status(f"CSV export error: {err}", is_error=True)

    def export_json(self):
        if not self.last_analysis:
            messagebox.showinfo("No Analysis", "Please run sequence analysis before exporting.")
            return
            
        default_name = f"dna_analysis_{self.last_analysis['sequence_length']}bp.json"
        path = filedialog.asksaveasfilename(
            title="Export Analysis as JSON",
            defaultextension=".json",
            initialfile=default_name,
            filetypes=[("JSON Files", "*.json"), ("All Files", "*.*")]
        )
        if not path:
            return
            
        try:
            export_results_json(self.last_analysis, path)
            messagebox.showinfo("Export Successful", f"Analysis data saved successfully to:\n{path}")
            self.set_status(f"Exported JSON: {os.path.basename(path)}")
        except Exception as err:
            messagebox.showerror("Export Failed", f"Could not save JSON file:\n{err}")
            self.set_status(f"JSON export error: {err}", is_error=True)

    def save_plot_image(self):
        if not self.last_analysis:
            messagebox.showinfo("No Analysis", "Please run sequence analysis before saving the graph.")
            return
            
        default_name = f"gc_distribution_{self.last_analysis['sequence_length']}bp.png"
        path = filedialog.asksaveasfilename(
            title="Save GC Content Graph",
            defaultextension=".png",
            initialfile=default_name,
            filetypes=[("PNG Image", "*.png"), ("PDF Document", "*.pdf"), ("All Files", "*.*")]
        )
        if not path:
            return
            
        try:
            self.fig.savefig(path, dpi=300, bbox_inches="tight")
            messagebox.showinfo("Graph Saved", f"GC Content graph saved successfully to:\n{path}")
            self.set_status(f"Saved graph: {os.path.basename(path)}")
        except Exception as err:
            messagebox.showerror("Save Failed", f"Could not save graph image:\n{err}")
            self.set_status(f"Graph save error: {err}", is_error=True)


def launch():
    root = tk.Tk()
    app = DNAAnalyzerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    launch()
