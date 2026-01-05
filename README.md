# cofolding-pulldown-tools

A small command-line toolkit for fetching and workign with protein sequence files, designed to support co-folding and pulldown workflows (such as [AlphaPulldown](https://github.com/KosinskiLab/AlphaPulldown)).

The package currently provides:
- FASTA retrieval from UniProt (by query or accession list)
- FASTA header cleaning (UniProt → simplified headers)
- FASTA slicing into fixed-length overlapping windows

Additional functionality (input generation, analysis and plotting of results etc.) is planned.
## Installation (users)

You can install directly from github with **uv tool**, which will install this package in an isolated environment and add the **cpt** command to your global path. 

```bash
uv tool install git+https://github.com/ayushd718/cofolding-pulldown-tools.git
```

Alternative you can install directly from github in a clean virtual environment of your choice with pip.

```
pip install git+https://github.com/ayushd718/cofolding-pulldown-tools.git
```

## Command line usage

### Fetch FASTA from UniProt

#### by query

```bash
cpt fetch query \
    --query p53 \
    --taxonomy_id 9606 \
    --reviewed true
```
- --query takes a UniProt query string
- --taxonomy_id takes an NCBI taxonomy_ID (e.g. 9606 for humans)
- --reviewed takes true or false

Any subset of these argumens can be provided, but at least one is required. Outputs for the above command would be p53_9606_reviewed.fasta

#### by accession list

```bash
cpt fetch acc --file accessions.txt
```
The accessions.txt file must contain UniProt accessions (e.g. P04637 for human p53) separated by line. Outputs for the above command would be accessions.fasta 

### FASTA processing

#### cleaning headers and ensuring FASTA sequenes
```bash
cpt fasta clean --file file.fasta
```
Converts headers from: 
```
>sp|P04637|P53_HUMAN ...
```
into:
```
>P04637_P53_HUMAN
```
The output, file_cleaned.fasta, will have each FASTA sequence on a single line (not wrapped).
Inputs are assumed to have UniProt-style pipe delimited headers.

#### slicing FASTA sequences
```bash
cpt fasta fetch \
    --file file.fasta
    --max_slice 500
    --window 50
```
The above command will go through file.fasta and output file_sliced.fasta, which will contain evenly slices fasta sequences for sequences larger than 500 with overlapping windows of length 50 between slices. 

The main function of this tool is to allow for easy pre-processing of fasta files for large in-silico cofolding screens. 

Currently all outputs are written to the working directory. 

## Installation (development)

This project uses **uv** for dependency and environment management. You can find installation instructions for it [here.](https://github.com/astral-sh/uv)

Clone the repository and sync dependencies:

```bash
git clone https://github.com/ayushd718/cofolding-pulldown-tools
cd cofolding-pulldown-tools
uv sync
```
Run CLI and look at subcommands with:

```bash
uv run cpt --help
```
## License

MIT License 

## References

- Ahmad S, da Costa Gonzales L J, Bowler-Barnett E H, Rice D L, Kim M, Wijerathne S, Luciani A, Kandasaamy S, Luo J, Watkins X, Turner E, Martin M J, UniProt Consortium The UniProt website API: facilitating programmatic access to protein knowledge Nucleic Acids Research, gkaf394 (2025)

- AlphaPulldown2—a general pipeline for high-throughput structural modeling: https://github.com/KosinskiLab/AlphaPulldown

- uv: https://github.com/astral-sh/uv
