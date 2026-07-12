"""Deduplicate near-identical sequences by clustering with mmseqs2.

Built for `slice_fasta`/`reformat_sliced_windows` output: a giant tandem-repeat
protein (LDL-repeat receptors, cadherin-repeat proteins, mucin VNTR domains,
...) sliced into fixed windows produces many near-identical windows, one per
repeat copy. Folding every copy independently wastes compute on structures
that are (up to the repeat's natural sequence drift) the same job run twice.
Clustering by sequence identity and keeping one representative per cluster
collapses those down to one fold job per distinct repeat, while proteins
without close-repeat structure keep effectively all of their windows.

Needs the `mmseqs` binary on PATH. It's a compiled tool, not a pip package --
install via conda/mamba:
    conda install -c conda-forge -c bioconda mmseqs2
"""
import csv
import os
import shutil
import subprocess
import tempfile


def _require_mmseqs():
    if shutil.which("mmseqs") is None:
        raise SystemExit(
            "Deduplication needs the mmseqs binary on PATH. Install with:\n"
            "    conda install -c conda-forge -c bioconda mmseqs2"
        )


def _count_records(fasta_path):
    with open(fasta_path) as fh:
        return sum(1 for line in fh if line.startswith(">"))


def dedup_sequences(fasta_file: os.PathLike, min_seq_id: float = 0.95,
                     coverage: float = 0.8, cov_mode: int = 1):
    """Cluster `fasta_file` by sequence identity (`mmseqs easy-cluster`) and
    keep one representative sequence per cluster.

    `min_seq_id` is the identity threshold for two sequences to cluster
    together (0.95 = 95% identical); `coverage`/`cov_mode` control how much
    of each sequence must align (default cov_mode=1: coverage of the shorter
    of the two sequences, so a window fully contained in a longer duplicate
    still clusters).

    Writes `{root}_dedup{ext}` (representative sequences) and
    `{root}_dedup_clusters.tsv` (representative<TAB>member, one row per input
    sequence -- use this to trace a dropped duplicate back to its
    representative). Returns (rep_fasta_path, cluster_tsv_path).
    """
    _require_mmseqs()

    abs_path = os.path.abspath(fasta_file)
    root, ext = os.path.splitext(abs_path)
    out_prefix = f"{root}_dedup"

    with tempfile.TemporaryDirectory(prefix="mmseqs_tmp_") as tmp_dir:
        result = subprocess.run(
            ["mmseqs", "easy-cluster", abs_path, out_prefix, tmp_dir,
             "--min-seq-id", str(min_seq_id),
             "-c", str(coverage),
             "--cov-mode", str(cov_mode)],
            capture_output=True, text=True,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"mmseqs easy-cluster failed (exit {result.returncode}):\n"
                f"{result.stderr}"
            )

    rep_fasta = f"{out_prefix}_rep_seq.fasta"
    cluster_tsv = f"{out_prefix}_cluster.tsv"
    all_seqs = f"{out_prefix}_all_seqs.fasta"
    if os.path.exists(all_seqs):
        os.remove(all_seqs)

    n_in = _count_records(abs_path)
    n_out = _count_records(rep_fasta)
    print(f"deduped fasta written to {rep_fasta} ({n_in} -> {n_out} sequences, "
          f"min_seq_id={min_seq_id})")
    print(f"cluster membership written to {cluster_tsv}")
    return rep_fasta, cluster_tsv
