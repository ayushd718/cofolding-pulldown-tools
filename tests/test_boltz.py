import os
import csv
import pytest
from cofolding_pulldown_tools.boltz import (
    parse_smiles,
    write_yaml,
    generate_boltz_inputs,
)

SMILES_FILE = """# a comment
>ligA
CCO
>ligB
c1ccccc1
"""

FASTA = """>sp|P11111|AAA_HUMAN Protein A OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=1
MKTAYIAKQR
>sp|P22222|BBB_HUMAN Protein B OS=Homo sapiens OX=9606 GN=BBB PE=1 SV=1
MSEQVERYLONGENOUGHTOSKIP
"""


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_parse_smiles_first_and_named(tmp_path):
    sm = _write(tmp_path, "l.smiles", SMILES_FILE)
    assert parse_smiles(sm) == ("ligA", "CCO")
    assert parse_smiles(sm, "ligB") == ("ligB", "c1ccccc1")


def test_parse_smiles_missing_raises(tmp_path):
    sm = _write(tmp_path, "l.smiles", SMILES_FILE)
    with pytest.raises(ValueError, match="not found"):
        parse_smiles(sm, "nope")


def test_write_yaml_schema(tmp_path):
    out = tmp_path / "job.yaml"
    write_yaml(str(out), "MKTA", "CCO", msa_field="/x/P1.a3m", affinity=True)
    text = out.read_text()
    assert "version: 1" in text
    assert "      sequence: MKTA" in text
    assert "      msa: /x/P1.a3m" in text
    assert "      smiles: 'CCO'" in text
    assert "  - affinity:" in text and "      binder: B" in text


def test_write_yaml_no_msa_no_affinity(tmp_path):
    out = tmp_path / "job.yaml"
    write_yaml(str(out), "MKTA", "CCO", msa_field=None, affinity=False)
    text = out.read_text()
    assert "msa:" not in text
    assert "affinity" not in text


def test_generate_inputs_precomputed_and_skip(tmp_path):
    fa = _write(tmp_path, "seqs.fasta", FASTA)
    sm = _write(tmp_path, "l.smiles", SMILES_FILE)
    out_dir = tmp_path / "boltz_inputs"
    # max_len=15 -> the 24-residue protein B is skipped, A (10) is written
    generate_boltz_inputs(fa, sm, ligand_name="ligA", out_dir=str(out_dir),
                          msa_mode="precomputed", msa_dir=str(tmp_path / "msa"),
                          max_len=15)

    yamls = sorted(f for f in os.listdir(out_dir) if f.endswith(".yaml"))
    assert yamls == ["P11111_AAA__ligA.yaml"]
    written = (out_dir / "P11111_AAA__ligA.yaml").read_text()
    assert os.path.join(str(tmp_path / "msa"), "P11111.a3m") in written

    rows = list(csv.DictReader(open(out_dir / "manifest.csv")))
    status = {r["id"]: r["status"] for r in rows}
    assert status == {"P11111": "written", "P22222": "skipped_long"}


def test_generate_inputs_msa_modes(tmp_path):
    fa = _write(tmp_path, "seqs.fasta", FASTA)
    sm = _write(tmp_path, "l.smiles", SMILES_FILE)

    empty_dir = tmp_path / "empty_out"
    generate_boltz_inputs(fa, sm, out_dir=str(empty_dir), msa_mode="empty", max_len=0)
    assert "      msa: empty" in (empty_dir / "P11111_AAA__ligA.yaml").read_text()

    server_dir = tmp_path / "server_out"
    generate_boltz_inputs(fa, sm, out_dir=str(server_dir), msa_mode="server", max_len=0)
    assert "msa:" not in (server_dir / "P11111_AAA__ligA.yaml").read_text()
