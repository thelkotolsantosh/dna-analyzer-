# Changelog

All notable changes to **dna-analyzer** are documented here.

---

## [1.0.0] — 2024-01-01

### Added
- `analyze` command: full pipeline — summary, nucleotide frequency, GC content, ORFs, dinucleotide heatmap
- `motif` command: find all occurrences of a motif with context display
- `orfs` command: list ORFs in all 6 reading frames with optional protein translation
- Color-coded sequence map for first 120 bp
- Dinucleotide frequency heatmap
- Melting temperature estimation (Wallace rule + GC% formula)
- Reverse complement calculation
- FASTA file input support
- stdin pipe support
- Full unit test suite with pytest
- GitHub Actions CI across Python 3.10 / 3.11 / 3.12
- Zero external runtime dependencies

---

[1.0.0]: https://github.com/thelkotolsantosh/dna-analyzer/releases/tag/v1.0.0
