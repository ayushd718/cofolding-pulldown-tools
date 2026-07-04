import os
import cofolding_pulldown_tools.fetch as fetch
from cofolding_pulldown_tools.fetch import (
    all_fasta_records,
    first_fasta_record,
    fetch_fasta_by_gene_names,
    _record_primary_gene,
    _select_record_for_gene,
)

RECORD_A2M = (
    ">sp|P01023|A2MG_HUMAN Alpha-2-macroglobulin OS=Homo sapiens OX=9606 GN=A2M PE=1 SV=3\n"
    "MGKNKLLHPSLVLLLLVLLPTDASV\n"
)
RECORD_TR = (
    ">tr|X0X0X0|X0X0X0_HUMAN Unreviewed thing OS=Homo sapiens OX=9606 GN=FOO PE=4 SV=1\n"
    "MSEQUENCEUNREVIEWED\n"
)


def test_first_fasta_record_single():
    assert first_fasta_record(RECORD_A2M) == RECORD_A2M


def test_first_fasta_record_multi_returns_first():
    two = RECORD_A2M + ">sp|P99999|OTHER_HUMAN x\nMMMM\n"
    assert first_fasta_record(two).startswith(">sp|P01023|A2MG_HUMAN")
    assert "P99999" not in first_fasta_record(two)


def test_first_fasta_record_empty():
    assert first_fasta_record("") is None
    assert first_fasta_record("not a fasta") is None


def test_all_fasta_records_splits():
    two = RECORD_A2M + ">sp|P99999|OTHER_HUMAN x\nMMMM\n"
    recs = all_fasta_records(two)
    assert len(recs) == 2
    assert recs[0].startswith(">sp|P01023|A2MG_HUMAN")
    assert recs[1].startswith(">sp|P99999|OTHER_HUMAN")
    assert all_fasta_records("") == []


def test_record_primary_gene():
    assert _record_primary_gene(RECORD_A2M) == "A2M"
    assert _record_primary_gene(">sp|X|Y_HUMAN no gene here\nMM\n") == ""


def test_select_record_prefers_primary_gene_match():
    # gene_exact:CCR4 returns the NOCT entry first, then the real CCR4 entry.
    noct = ">sp|Q9UK39|NOCT_HUMAN Nocturnin OS=Homo sapiens OX=9606 GN=NOCT PE=1 SV=1\nMNOCT\n"
    ccr4 = ">sp|P51679|CCR4_HUMAN C-C chemokine receptor 4 OS=Homo sapiens OX=9606 GN=CCR4 PE=1 SV=1\nMCCR4\n"
    chosen = _select_record_for_gene([noct, ccr4], "CCR4")
    assert chosen.startswith(">sp|P51679|CCR4_HUMAN")


def test_select_record_falls_back_to_first_on_rename():
    # A renamed symbol (no primary-GN match) should take the top hit.
    aars1 = ">sp|P49588|SYAC_HUMAN Alanine--tRNA ligase OS=Homo sapiens OX=9606 GN=AARS1 PE=1 SV=2\nMAARS\n"
    chosen = _select_record_for_gene([aars1], "AARS")
    assert chosen.startswith(">sp|P49588|SYAC_HUMAN")


def test_fetch_by_gene_names_reviewed(tmp_path, monkeypatch):
    genes = tmp_path / "genes.txt"
    genes.write_text("A2M\nMISSING\n")

    def fake(gene, organism_id, reviewed, retries=3):
        if gene == "A2M" and reviewed == "true":
            return RECORD_A2M
        return None

    monkeypatch.setattr(fetch, "_fetch_gene_record", fake)
    fetch_fasta_by_gene_names(str(genes), delay=0)

    root = os.path.splitext(str(genes))[0]
    assert os.path.exists(f"{root}.fasta")
    assert ">sp|P01023|A2MG_HUMAN" in open(f"{root}.fasta").read()

    report = open(f"{root}_report.tsv").read()
    assert "A2M\tP01023\treviewed" in report
    assert "MISSING\t\tfailed" in report

    assert "MISSING" in open(f"{root}_failed.txt").read()


def test_fetch_by_gene_names_unreviewed_fallback(tmp_path, monkeypatch):
    genes = tmp_path / "genes.txt"
    genes.write_text("FOO\n")

    def fake(gene, organism_id, reviewed, retries=3):
        return RECORD_TR if reviewed == "false" else None

    monkeypatch.setattr(fetch, "_fetch_gene_record", fake)
    fetch_fasta_by_gene_names(str(genes), delay=0)

    root = os.path.splitext(str(genes))[0]
    assert "FOO\tX0X0X0\tunreviewed" in open(f"{root}_report.tsv").read()


def test_fetch_by_gene_names_reviewed_only_skips_fallback(tmp_path, monkeypatch):
    genes = tmp_path / "genes.txt"
    genes.write_text("A2M\nFOO\n")

    def fake(gene, organism_id, reviewed, retries=3):
        if gene == "A2M" and reviewed == "true":
            return RECORD_A2M
        return RECORD_TR if reviewed == "false" else None

    monkeypatch.setattr(fetch, "_fetch_gene_record", fake)
    fetch_fasta_by_gene_names(str(genes), reviewed_only=True, delay=0)

    root = os.path.splitext(str(genes))[0]
    report = open(f"{root}_report.tsv").read()
    assert "A2M\tP01023\treviewed" in report
    assert "FOO\t\tfailed" in report  # fallback skipped
