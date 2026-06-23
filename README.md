# cofolding-pulldown-tools

A command-line toolkit designed to facilitate high-throughput protein sequence preparation for co-folding and PPI pulldown workflows, such as [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown).

## Key Features

- **Automated FASTA Retrieval**: Fetch sequences directly from UniProt via query strings or accession lists.
- **Header Standardization**: Convert complex UniProt headers into simplified, machine-readable formats.
- **Sequence Slicing**: Automatically slice long sequences into fixed-length overlapping windows to suit modeling constraints.
- **Complex Generation**: Generate bait-prey pair files for high-throughput interaction screens.

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
| `cpt fasta clean` | Simplify FASTA headers and format sequence lines |
| `cpt fasta slice` | Slice large sequences into overlapping windows |
| `cpt fasta complex` | Generate a bait-prey pairing file |

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
    --bait 'protein1;protein2' \
    --double_count
```
Generates a `.txt` file mapping bait proteins to all sequences in the input FASTA.

---
*Note: All outputs are currently written to the working directory.*

## License
MIT License

## References
- Ahmad S, et al. "The UniProt website API: facilitating programmatic access to protein knowledge." Nucleic Acids Research, 2025.
- AlphaPulldown2: https://github.com/KosinskiLab/AlphaPulldown
- uv: https://github.com/astral-sh/uv
