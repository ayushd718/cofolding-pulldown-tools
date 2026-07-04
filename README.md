# cofolding-pulldown-tools

A command-line toolkit designed to facilitate high-throughput protein sequence preparation for co-folding and PPI pulldown workflows, such as [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown).

## Key Features

- **Automated FASTA Retrieval**: Fetch sequences directly from UniProt via query strings, accession lists, or gene-symbol lists.
- **Header Standardization**: Convert complex UniProt headers into simplified, machine-readable formats.
- **Sequence Slicing**: Automatically slice long sequences into fixed-length overlapping windows to suit modeling constraints.
- **Domain Splitting**: Split proteins into structure-based domain fragments from the AlphaFold DB PAE (coherent folding units rather than fixed windows).
- **Complex Generation**: Generate bait-prey pair files for high-throughput interaction screens.
- **Boltz-2 Inputs**: Generate protein × ligand cofolding YAMLs with affinity prediction.
- **GPU MSAs**: Batch MSA generation via GPU-accelerated MMseqs2 (`scripts/run_msa_mmseqs_gpu.sh`).
- **Result Analysis**: Analyze cofolding outputs to generate summary CSVs including ipSAE, pDockQ, and iPTM scores.

## Requirements

- Python >= 3.10
- [uv](https://github.com/astral-sh/uv) (recommended for installation and dependency management)

## Installation & Usage

### For Users
You can install the tool globally using **uv tool**, which runs it in an isolated environment and adds the `cpt` command to your path:

```bash
uv tool install git+https://github.com/ayushd718/cofolding-pulldown-tools.git
```

Alternatively, use `pip`:

```bash
pip install git+https://github.com/ayushd718/cofolding-pulldown-tools.git
```

Use the `cpt` command for operations:
```bash
cpt --help
```

### For Developers
Clone the repository and sync dependencies using `uv`:

```bash
git clone https://github.com/ayushd718/cofolding-pulldown-tools
cd cofolding-pulldown-tools
uv sync
```

Use `uv run` to execute the CLI within the development environment:
```bash
uv run cpt --help
```

## CLI Reference

| Subcommand | Action |
| :--- | :--- |
| `cpt fetch query` | Fetch FASTA by query string/taxonomy/reviewed status |
| `cpt fetch acc` | Fetch FASTA by a file containing accession IDs |
| `cpt fetch genes` | Fetch one canonical sequence per gene symbol |
| `cpt fasta clean` | Simplify FASTA headers and format sequence lines |
| `cpt fasta slice` | Slice large sequences into overlapping windows |
| `cpt fasta complex` | Generate a bait-prey pairing file |
| `cpt fasta domains` | Split proteins into domain fragments using AlphaFold DB PAE |
| `cpt boltz inputs` | Generate Boltz-2 cofolding input YAMLs (protein × ligand) |
| `cpt analyze` | Analyze cofolding results and generate a summary CSV |

### Fetching from UniProt

#### By Query
```bash
cpt fetch query \
    --query p53 \
    --taxonomy_id 9606 \
    --reviewed true
```

#### By Accession List
```bash
cpt fetch acc --file accessions.txt
```

#### By Gene Symbol
For a list of gene symbols (one per line), fetch the reviewed (Swiss-Prot)
canonical human sequence for each. Symbols with no reviewed entry fall back to
the best unreviewed (TrEMBL) entry unless `--reviewed_only` is set.
```bash
cpt fetch genes --file genes.txt --organism_id 9606
```
Writes `genes.fasta` (raw UniProt headers, ready for `cpt fasta clean`),
`genes_report.tsv` (gene → accession, status, header), and `genes_failed.txt`
(unresolved symbols, which can be re-run through the same command).

### FASTA Processing

#### Cleaning
```bash
cpt fasta clean --file file.fasta
```
Converts `>sp|P04637|P53_HUMAN` to `>P04637_P53_HUMAN`. Ensures sequences are on a single line.

#### Slicing
```bash
cpt fasta slice \
    --file file.fasta \
    --max_slice 400 \
    --window 50
```
Slices sequences larger than `max_slice`, using a `window` for overlap.

#### Generating Complexes
```bash
cpt fasta complex \
    --file file.fasta \
    --bait 'protein1,protein2...' \
    --double_count
```
Generates a `file_complex.txt` file mapping bait proteins to all sequences in the input FASTA. This function validates that all provided bait proteins are present as headers in the input FASTA file.

#### Splitting by Domain (AlphaFold DB)
Splits each protein into structure-based domain fragments by clustering the
AlphaFold DB PAE matrix (residues AlphaFold places confidently relative to each
other become one domain). A better-founded alternative to fixed-width slicing —
fragments are coherent folding units rather than arbitrary windows.
```bash
cpt fasta domains --file file.fasta --ndr exclude
```
`--ndr` decides what happens to non-domain residues (disordered tails/linkers):
`exclude` drops them (default), `keep` emits each linker as its own fragment,
`pad` extends domains into flanking linker. `--resolution` tunes granularity
(0.1–0.5 stable; higher over-splits). Writes `<file>_domains.fasta` and a
`<file>_domains.tsv` map; proteins whose AFDB model is fragmented (very long
sequences) are logged to `_multifragment.txt` for separate handling.

> Requires the optional dependencies: `pip install "cofolding-pulldown-tools[domains]"`

### Generating Boltz-2 Inputs
Writes one Boltz-2 YAML per protein, each paired with a ligand (by SMILES), with
affinity prediction enabled by default. MSAs are referenced as
`<msa_dir>/<id>.a3m` — the naming produced by `scripts/run_msa_mmseqs_gpu.sh`.
```bash
cpt boltz inputs --file file.fasta --ligands ligands.smiles --ligand_name myligand
```
`--msa_mode` is `precomputed` (default), `empty` (single-sequence), or `server`
(fetch at predict time with `boltz predict --use_msa_server`). Proteins longer
than `--max_len` are skipped and recorded in `manifest.csv`.

### Generating MSAs (GPU MMseqs2)
`scripts/run_msa_mmseqs_gpu.sh` drives ColabFold's GPU MMseqs2 search to produce
one `<accession>.a3m` per sequence (the naming the Boltz inputs expect). Requires
an NVIDIA GPU, a GPU build of MMseqs2, `colabfold_search`, and the ColabFold
databases; see the header of the script for setup.
```bash
scripts/run_msa_mmseqs_gpu.sh sequences.fasta $HOME/colabfold_db msa
```

### Analyzing Results
```bash
cpt analyze --dir /path/to/structures_output
```
Processes all subdirectories in the specified directory to compute ipSAE and pDockQ scores for predicted structures, aggregating them with iPTM data into a `results_summary.csv` file.

---
*Note: All outputs are currently written to the working directory.*

## License
MIT License

## References
- Ahmad S, et al. "The UniProt website API: facilitating programmatic access to protein knowledge." Nucleic Acids Research, 2025.
- AlphaPulldown2: https://github.com/KosinskiLab/AlphaPulldown
- uv: https://github.com/astral-sh/uv
