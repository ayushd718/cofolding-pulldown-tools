import os
import shutil

import pytest

from cofolding_pulldown_tools.cluster import dedup_sequences

pytestmark = pytest.mark.skipif(
    shutil.which("mmseqs") is None, reason="mmseqs binary not on PATH"
)

# Three near-identical windows (a repeated domain sliced 3x, differing by a
# couple of positions -- simulating tandem-repeat window duplication) plus one
# clearly distinct sequence that should survive dedup on its own.
REPEATED_UNIT = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ"
FASTA = "\n".join([
    ">rep_1.w1-67 GN=REPGENE",
    REPEATED_UNIT,
    ">rep_2.w70-137 GN=REPGENE",
    REPEATED_UNIT,
    ">rep_3.w140-207 GN=REPGENE",
    REPEATED_UNIT[:-1] + "A",  # one residue different, still >95% identical
    ">unique_1.w1-56 GN=OTHERGENE",
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVP",
]) + "\n"


def test_dedup_sequences_collapses_near_identical_windows(tmp_path):
    fasta_file = tmp_path / "windows.fasta"
    fasta_file.write_text(FASTA)

    rep_fasta, cluster_tsv = dedup_sequences(fasta_file, min_seq_id=0.95)

    assert os.path.exists(rep_fasta)
    assert os.path.exists(cluster_tsv)

    with open(rep_fasta) as f:
        headers = [line for line in f if line.startswith(">")]
    # the 3 near-identical repeat windows collapse to 1 representative;
    # the unrelated sequence survives as its own representative
    assert len(headers) == 2

    with open(cluster_tsv) as f:
        rows = [line.rstrip("\n").split("\t") for line in f]
    assert len(rows) == 4  # one row per input sequence

    os.remove(rep_fasta)
    os.remove(cluster_tsv)
