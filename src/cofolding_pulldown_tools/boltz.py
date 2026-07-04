"""Generate Boltz-2 cofolding input YAMLs: each protein x one ligand.

For every protein in a FASTA this writes one Boltz-2 job pairing it with a chosen
ligand (by SMILES), with the affinity-prediction block enabled by default so hits
can be ranked by predicted binding, not just interface confidence. MSAs are
referenced as <msa-dir>/<id>.a3m (the naming produced by the MMseqs2 MSA script),
or omitted for single-sequence / server modes.
"""
import csv
import os
import re

from .fasta import iter_fasta


def parse_smiles(path, want=None):
    """Return (name, smiles) for ligand `want` (or the first) from a .smiles file
    using the `>name` header convention."""
    name, first = None, None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                name = line[1:].strip()
            elif name is not None:
                if first is None:
                    first = (name, line)
                if want is None or name == want:
                    return name, line
                name = None
    if want is not None:
        raise ValueError(f"Ligand '{want}' not found in {path}")
    if first is None:
        raise ValueError(f"No ligands parsed from {path}")
    return first


def _safe(s):
    return re.sub(r"[^A-Za-z0-9_.-]", "-", s)


def write_yaml(path, sequence, smiles, msa_field=None, affinity=True):
    """Write one Boltz-2 YAML (protein chain A + ligand chain B)."""
    lines = ["version: 1", "sequences:", "  - protein:", "      id: A",
             f"      sequence: {sequence}"]
    if msa_field is not None:
        lines.append(f"      msa: {msa_field}")
    lines += ["  - ligand:", "      id: B", f"      smiles: '{smiles}'"]
    if affinity:
        lines += ["properties:", "  - affinity:", "      binder: B"]
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def generate_boltz_inputs(fasta_file, ligands_file, ligand_name=None,
                          out_dir="boltz_inputs", msa_mode="precomputed",
                          msa_dir="msa", max_len=2500, affinity=True):
    """Write one Boltz-2 YAML per protein in `fasta_file` paired with a ligand
    from `ligands_file`, plus a manifest.csv. Proteins longer than max_len
    (0 = no limit) are skipped and recorded in the manifest.

    msa_mode: 'precomputed' (msa=<msa_dir>/<id>.a3m), 'empty' (single-sequence),
    or 'server' (omit msa; run `boltz predict --use_msa_server`)."""
    lig_name, smiles = parse_smiles(ligands_file, ligand_name)
    os.makedirs(out_dir, exist_ok=True)
    msa_dir_abs = os.path.abspath(msa_dir)

    manifest, n_written, n_skipped = [], 0, 0
    for ident, gene, seq in iter_fasta(fasta_file):
        if max_len and len(seq) > max_len:
            manifest.append(dict(id=ident, gene=gene, ligand=lig_name,
                                 seq_len=len(seq), yaml="", status="skipped_long"))
            n_skipped += 1
            continue
        if msa_mode == "precomputed":
            msa_field = os.path.join(msa_dir_abs, f"{ident}.a3m")
        elif msa_mode == "empty":
            msa_field = "empty"
        else:
            msa_field = None
        job = f"{_safe(ident)}_{_safe(gene)}__{_safe(lig_name)}"
        yaml_path = os.path.join(out_dir, f"{job}.yaml")
        write_yaml(yaml_path, seq, smiles, msa_field, affinity)
        manifest.append(dict(id=ident, gene=gene, ligand=lig_name,
                             seq_len=len(seq), yaml=yaml_path, status="written"))
        n_written += 1

    man_path = os.path.join(out_dir, "manifest.csv")
    with open(man_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["id", "gene", "ligand",
                                                "seq_len", "yaml", "status"])
        writer.writeheader()
        writer.writerows(manifest)

    print(f"ligand {lig_name}: {smiles}")
    print(f"wrote {n_written} Boltz YAMLs to {out_dir}/ (msa_mode={msa_mode})"
          + (f", skipped {n_skipped} > {max_len} aa" if n_skipped else ""))
    print(f"run: boltz predict {out_dir}/ --out_dir boltz_out --output_format pdb"
          + (" --use_msa_server" if msa_mode == "server" else ""))
    return man_path
