# 🧬 dna-analyzer > A terminal-based DNA sequence analysis toolkit — bioinformatics in your terminal, no external dependencies.

---

## ✨ Features

| Feature | Description |
|---|---|
| 📊 **Nucleotide Frequency** | Color-coded bar chart of A / T / C / G counts |
| 🌡️ **GC Content & Melting Temp** | GC% meter + Wallace-rule Tm estimation |
| 🔬 **ORF Detection** | Finds Open Reading Frames in all 6 reading frames |
| 🔍 **Motif Search** | Locate any subsequence with highlighted context |
| 🗺️ **Sequence Map** | Color-coded base-by-base view (first 120 bp) |
| 🔥 **Dinucleotide Heatmap** | 4×4 frequency matrix of all dinucleotide pairs |
| 🔄 **Reverse Complement** | Instant complement & reverse |
| 📄 **FASTA Support** | Read `.fasta` files or plain text files |
| 🔗 **Pipe-friendly** | Works with Unix pipes: `cat seq.txt \| dna-analyzer analyze` |
| ⚡ **Zero Dependencies** | Pure Python stdlib — nothing to install beyond the package |

---

## 🚀 Installation

### From PyPI *(coming soon)*

```bash
pip install dna-analyzer
```

### From source

```bash
git clone https://github.com/thelkotolsantosh/dna-analyzer.git
cd dna-analyzer
pip install -e .
```

After installation, the `dna-analyzer` command is available globally.

---

## 📖 Usage

### `analyze` — Full Analysis Pipeline

```bash
# From a sequence string
dna-analyzer analyze -s "ATGCGATCGATCGAAATTTGCGATCGATCG"

# From a FASTA file
dna-analyzer analyze -f examples/sample_sequences.fasta

# From stdin
cat my_sequence.txt | dna-analyzer analyze

# Custom minimum ORF length (default: 30 bp)
dna-analyzer analyze -f genome.fasta --min-orf 60

# Skip dinucleotide heatmap for faster output
dna-analyzer analyze -s "ATCGATCG..." --no-dinuc
```

### `motif` — Find Motif Occurrences

```bash
dna-analyzer motif -s "ATGCGATCGATGATCGATG" -m "ATG"
# Output shows each position with surrounding context highlighted
```

### `orfs` — List Open Reading Frames

```bash
dna-analyzer orfs -f sequence.fasta --min-orf 90

# Show translated amino acid sequences for top 5 ORFs
dna-analyzer orfs -f sequence.fasta --proteins
```

---

## 🖥️ Example Output

```
╔══════════════════════════════════════════════════════════╗
║                 🧬  DNA SEQUENCE ANALYZER                 ║
╚══════════════════════════════════════════════════════════╝

  ▶ Sequence Summary
  ──────────────────────────────────────────────────
  Length                 142 bp
  GC Content             48.59%
  AT Content             51.41%
  Melting Temp (Tm)      57.73 °C
  ORFs Found             3
  Nucleotide counts  A:37  T:36  C:34  G:35

  ▶ Nucleotide Frequency
  ──────────────────────────────────────────────────
  A  ████████████████░░░░░░░░░░░░░░     37  (26.1%)
  T  ███████████████░░░░░░░░░░░░░░░     36  (25.4%)
  C  ██████████████░░░░░░░░░░░░░░░░     34  (23.9%)
  G  ███████████████░░░░░░░░░░░░░░░     35  (24.6%)

  ▶ GC Content
  ──────────────────────────────────────────────────
  G+C  ████████████████████░░░░░░░░░░░░░░░░░░░░   48.6%
  A+T  █████████████████████░░░░░░░░░░░░░░░░░░░   51.4%

  Stability: Balanced GC content

  ▶ Open Reading Frames  (3 found)
  ──────────────────────────────────────────────────
  #1   [·····████████████████·················]  Frame +1    7–54    48 bp
  #2   [·············██████···················]  Frame -2   89–112   24 bp
  #3   [··············████·····················]  Frame +3  102–117  16 bp
```

---

## 🗂️ Project Structure

```
dna-analyzer/
├── dna_analyzer/
│   ├── __init__.py         # Public API exports
│   ├── analyzer.py         # Core analysis functions
│   ├── visualizer.py       # Terminal color rendering
│   └── cli.py              # argparse CLI entry point
├── tests/
│   ├── __init__.py
│   └── test_analyzer.py    # pytest unit tests
├── examples/
│   └── sample_sequences.fasta
├── .github/
│   └── workflows/
│       └── ci.yml          # GitHub Actions CI
├── .gitignore
├── CHANGELOG.md
├── CONTRIBUTING.md
├── LICENSE
├── README.md
├── pyproject.toml
├── requirements.txt
└── requirements-dev.txt
```

---

## 🧪 Running Tests

```bash
# Install dev dependencies
pip install -r requirements-dev.txt

# Run all tests
pytest

# With coverage report
pytest --cov=dna_analyzer --cov-report=term-missing

# Run linter
ruff check dna_analyzer/ tests/
```

---

## 🔬 Science Behind the Features

### GC Content
GC base pairs form **3 hydrogen bonds** (vs 2 for AT), making GC-rich sequences thermally stable. High GC > 60% indicates thermostable DNA; AT-rich regions often contain regulatory elements.

### Melting Temperature
- **Short sequences (<14 bp):** Wallace Rule → `Tm = 2(A+T) + 4(G+C)`
- **Longer sequences:** `Tm = 81.5 + 0.41 × GC% − 675/n`

### Open Reading Frames (ORFs)
Scans all **6 reading frames** (3 on each strand) for ATG → STOP codon spans. The reverse complement strand is analyzed using Watson-Crick complement rules.

### Codon Table
Uses the **standard genetic code** (NCBI translation table 1).

---

## 🤝 Contributing

Contributions are welcome! See [CONTRIBUTING.md](CONTRIBUTING.md) for guidelines on setting up your dev environment, running tests, and submitting pull requests.

---

## 📜 License

This project is licensed under the [MIT License](LICENSE).

---

## 💡 Roadmap

- [ ] Export results to JSON / CSV
- [ ] Multiple codon tables (mitochondrial, bacterial)
- [ ] Local alignment (Smith-Waterman)
- [ ] GFF / GenBank file parsing
- [ ] Web interface (FastAPI + htmx)
- [ ] Publish to PyPI
