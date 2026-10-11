# cofolding-pulldown-tools

A command-line toolkit designed to facilitate high-throughput protein sequence preparation for co-folding and PPI pulldown workflows, such as [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown).

## Key Features

- **Automated FASTA Retrieval**: Fetch sequences directly from UniProt via query strings, accession lists, or gene-symbol lists.
- **Header Standardization**: Convert complex UniProt headers into simplified, machine-readable formats.
- **Sequence Slicing**: Automatically slice long sequences into fixed-length overlapping windows to suit modeling constraints.
- **Domain Splitting**: Split proteins into structure-based domain fragments from the AlphaFold DB PAE (coherent folding units rather than fixed windows).
- **Annotation-Driven Domain Packing**: Pack whole UniProt-annotated domains into length-capped fragments for proteins AlphaFold DB has no model for (>2,700 aa), cutting only between domains so none is ever bisected.
- **Window Deduplication**: Collapse near-identical windowed fragments of tandem-repeat proteins (via `mmseqs easy-cluster`) so repeat domains aren't folded redundantly.
- **Complex Generation**: Generate bait-prey pair files for high-throughput interaction screens.
- **Boltz-2 Inputs**: Generate protein × ligand cofolding YAMLs with affinity prediction.
- **GPU MSAs**: Batch MSA generation via GPU-accelerated MMseqs2 (`scripts/run_msa_mmseqs_gpu.sh`).
- **Result Analysis**: Analyze cofolding outputs to generate summary CSVs including ipSAE, pDockQ, and iPTM scores.

## Requirements

- Python >= 3.10
- `pip`, or [uv](https://github.com/astral-sh/uv) (recommended for development)

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

Some subcommands need optional dependencies, declared as extras:

```bash
pip install "cofolding-pulldown-tools[domains] @ git+https://github.com/ayushd718/cofolding-pulldown-tools.git"
```

- `[domains]` adds numpy and networkx, needed by `cpt fasta domains`
- `[analysis]` adds ipsae, numpy and pandaprot, needed by `cpt analyze`
- `cpt fasta dedup` additionally needs the `mmseqs` binary, which is not a pip
  package: `conda install -c conda-forge -c bioconda mmseqs2`

Everything else, including `cpt fasta uniprot-domains`, needs only the base
install.

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

Run the test suite:
```bash
uv run pytest tests/ -q
```

### Publishing a Release
The project is a standard PEP 517 build, so a release needs no extra
configuration:

```bash
pip install build twine
python -m build            # writes dist/*.whl and dist/*.tar.gz
twine check dist/*         # validates the metadata PyPI will reject on
twine upload dist/*        # needs a PyPI API token
```

Bump `version` in `pyproject.toml` first. A version number can never be reused
on PyPI once uploaded, so `twine check` passing is worth confirming before the
upload step.

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
| `cpt fasta uniprot-domains` | Pack whole UniProt-annotated domains into length-capped fragments |
| `cpt fasta reformat-windows` | Rewrite `slice` output headers with real gene symbols |
| `cpt fasta dedup` | Cluster near-identical sequences and keep one representative per cluster |
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
`<file>_domains.tsv` map. Proteins whose AFDB model does not cover the full
sequence are logged to `_multifragment.txt`. In practice that means anything over
about 2,700 residues, which is the AlphaFold DB size ceiling: such entries have
**no model at all**, or only a short isoform model. Use
`cpt fasta uniprot-domains` for those.

> Requires the optional dependencies: `pip install "cofolding-pulldown-tools[domains]"`

#### Splitting by Domain (UniProt annotation)
Packs *whole* UniProt-annotated domains into fragments up to `--max_slice`,
cutting only between one domain's end and the next domain's start. No annotated
domain is ever bisected.

This exists because `cpt fasta domains` cannot help the proteins that most need
slicing. AlphaFold DB has no model above roughly 2,700 residues, so FBN1
(2,871 aa), HSPG2 (4,391 aa) and MUC16 (14,507 aa) return nothing at all, while
LAMA5, VWF and HMCN1 return only short isoform models covering residues 561,
351 and 555. UniProt annotation needs no structure model: all 37 proteins over
2,700 aa in a matrisome test set carry it, 1,174 domain and repeat features in
total.

```bash
cpt fasta uniprot-domains --file file.fasta --max_slice 1400
```

FBN1's 56 domains pack into **2 fragments** at `--max_slice 1400`, against 56
under naive per-domain splitting and 3 bisected windows under `cpt fasta slice`.

There is deliberately no attempt to find a safe residue inside a linker,
because there usually is no linker: the median inter-domain gap is 9 residues,
only 46.5% exceed 10, and inside repeat arrays it collapses to 1-5 residues,
with tenascin-C's EGF-like domains overlapping outright.

`--ndr` reuses the `exclude`/`keep`/`pad` vocabulary of `cpt fasta domains` and
adds `tile` (the default), which runs a tandem-period scan over each
unannotated stretch and emits a representative unit on a hit, falling back to
overlapping windows otherwise. Structural `Region` annotation such as
`Triple-helical region` is honoured and subtracted against what the domains
already cover, which lifts collagen coverage from 20% to 89%. `Region`
`Disordered` is deliberately ignored: it is predictor-derived and spans
COL1A1's rigid triple helix in full.

For mucin-heavy input, prefer `--ndr exclude`. MUC16's first 12,066 residues
are 83% of the protein, carry no annotation and have no recoverable
periodicity, so the default tiles them into 9 fragments of glycan-covered
sequence. Any unannotated stretch over `--big_ndr` is reported to stderr and to
`_large_unannotated.tsv` so the cost is visible. Across the 37-protein test set
`--ndr exclude` gives 107 fragments against 193.

`--collapse_repeats` groups tandem copies by family *and* qualifier, then
collapses only those above `--repeat_identity` (default 0.9), writing every
dropped copy and its representative to `<file>_repeat_clusters.tsv`. It is
**off by default, and the saving is smaller than it looks.** Measured pairwise
identity inside single-protein arrays: FBN1's 43 calcium-binding EGF-like
copies are a median 33% identical with zero pairs above 80%, and HMCN1's 44
Ig-like copies a median 23%. Nothing collapses at 0.9, and that is correct --
those are homologous folds with genuinely different binding surfaces. Only
mucin-type arrays are redundant (MUC16's SEA domains, median 72%). Across the
test set, collapsing removes 4,634 residues, 3.1% of the output, but raises the
fragment count from 193 to 235 because a dropped copy splits its run in two.
Do not lower the threshold below about 0.6; it starts discarding distinct
domain copies. The fragment-count win in this subcommand comes from packing,
not deduplication.

Writes `<out_dir>/<file>_uniprot_domains.fasta`, a `<file>_uniprot_domains.tsv`
fragment map and a `<file>_repeat_clusters.tsv` cluster map, plus
`_no_annotation.txt`, `_no_features.txt`, `_isoforms_skipped.txt` and
`_large_unannotated.tsv` where non-empty. Headers are
`{acc}_{gene}_{tag}_{start}-{end}` with 1-based inclusive coordinates and the
*real* gene symbol from UniProt, so unlike `slice` output this needs no
`cpt fasta reformat-windows` pass (and that subcommand will reject it, since it
requires a comma in the header).

Isoform accessions are skipped: UniProt serves them with a sequence but with no
feature annotation at all, and canonical coordinates do not apply to an isoform
sequence, so borrowing them would shift every boundary silently.

Annotation is cached one JSON per accession under
`<out_dir>/_annotation_cache`, so reruns are fully offline.

> Needs no optional dependencies, unlike `cpt fasta domains` and `cpt fasta dedup`.

#### Reformatting Windowed-Slice Headers
`cpt fasta clean` names each header after the UniProt *entry mnemonic*
(`>ACC_NAME`), which often diverges from the real gene symbol (`CAD23` is the
mnemonic for `CDH23`). After slicing, that mismatch propagates into
`>ACC_NAME,START-END`. This looks the true gene symbol up by accession from an
uncleaned reference FASTA and rewrites headers to `ACC.wSTART-END GN=GENE`:
```bash
cpt fasta reformat-windows \
    --file file_cleaned_sliced.fasta \
    --gene_map file.fasta
```
`--gene_map` is the *original*, uncleaned FASTA (the one `clean_fasta` was run
on) — it still has `GN=` tags. Writes `<file>_windows.fasta` plus a manifest
CSV (`window_id, accession, gene, start, end, window_len`). The period-joined
window tag (rather than underscore-joined) keeps `accession_gene__ligand`
job-name parsing (`cpt boltz inputs`, `cpt analyze`) working unmodified, since
an accession may itself contain a hyphen (isoform suffixes, e.g. `Q99102-10`)
but never a period.

#### Deduplicating Near-Identical Sequences
Slicing a giant tandem-repeat protein (LDL-repeat receptors, cadherin-repeat
proteins, mucin VNTR domains, ...) into fixed windows produces many
near-identical windows, one per repeat copy — folding every copy independently
wastes compute on the same job run twice. This clusters by sequence identity
(`mmseqs easy-cluster`) and keeps one representative per cluster:
```bash
cpt fasta dedup \
    --file file_windows.fasta \
    --min_seq_id 0.95 \
    --coverage 0.8 \
    --cov_mode 1
```
Writes `<file>_dedup.fasta` (representatives) and
`<file>_dedup_clusters.tsv` (representative↔member map, so a dropped
duplicate can be traced back to the window it collapsed into). Proteins
without close-repeat structure keep effectively all of their windows.

> Requires the `mmseqs` binary on `PATH` (not a pip package — install via
> `conda install -c conda-forge -c bioconda mmseqs2`).

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
