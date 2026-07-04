"""Split proteins into domain fragments using AlphaFold DB PAE.

For each UniProt accession this downloads the AlphaFold DB predicted-aligned-
error (PAE) matrix and clusters residues into domains by graph community
detection on the PAE (the pae_to_domains method: residues AlphaFold places
confidently relative to each other become one domain). Fragments are written as
a sliced FASTA in the >{ID}_{GENE}_{start}-{end} convention that the rest of the
toolkit (and Boltz input generation) consumes.

Needs the optional dependencies: pip install "cofolding-pulldown-tools[domains]"
"""
import csv
import json
import os
import time

import requests

from .fasta import iter_fasta

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


def _require_scicomp():
    try:
        import numpy as np
        import networkx as nx
    except ImportError as exc:
        raise SystemExit(
            "Domain splitting needs numpy + networkx. Install with:\n"
            "    pip install 'cofolding-pulldown-tools[domains]'"
        ) from exc
    return np, nx


# ---------------------------------------------------------------------------
# AlphaFold DB access
# ---------------------------------------------------------------------------
def _get(url, retries=3):
    for attempt in range(retries):
        try:
            resp = requests.get(url, timeout=120)
        except requests.RequestException:
            time.sleep(2 ** attempt)
            continue
        if resp.status_code == 200:
            return resp.text
        if resp.status_code == 404:
            return None
        time.sleep(2 ** attempt)
    return None


def fetch_pae(acc, cache_dir, retries=3):
    """Return the NxN PAE matrix for an accession, or None if AFDB has no model.
    The file version varies per entry, so the PAE URL is resolved through the
    prediction API and cached to disk (reruns are then offline)."""
    np, _ = _require_scicomp()
    cache = os.path.join(cache_dir, f"AF-{acc}-F1-pae.json")
    if os.path.exists(cache):
        with open(cache) as fh:
            text = fh.read()
    else:
        api = _get(AFDB_API.format(acc=acc), retries)
        if api is None:
            return None
        entries = json.loads(api)
        if not entries:
            return None
        entry = next((e for e in entries if e.get("entryId") == f"AF-{acc}-F1"), entries[0])
        pae_url = entry.get("paeDocUrl")
        if not pae_url:
            return None
        text = _get(pae_url, retries)
        if text is None:
            return None
        with open(cache, "w") as fh:
            fh.write(text)
    obj = json.loads(text)[0]
    if "predicted_aligned_error" in obj:
        return np.asarray(obj["predicted_aligned_error"], dtype=float)
    if "distance" in obj:                                   # legacy sparse schema
        n = int(max(obj["residue2"]))
        pae = np.zeros((n, n))
        pae[np.asarray(obj["residue1"]) - 1, np.asarray(obj["residue2"]) - 1] = obj["distance"]
        return pae
    raise ValueError(f"Unrecognized PAE JSON schema for {acc}")


# ---------------------------------------------------------------------------
# Clustering + segmentation (pure helpers, so they are unit-testable)
# ---------------------------------------------------------------------------
def cluster_pae(pae, pae_power=1.0, pae_cutoff=5.0, resolution=0.3, seed=0):
    """Community-detect residues into domains from the PAE. Returns a per-residue
    integer label list of length N. Lower resolution -> larger domains; 0.1-0.5
    is a stable range (1.0 tends to over-split single domains)."""
    np, nx = _require_scicomp()
    n = pae.shape[0]
    sym = np.minimum(pae, pae.T)                 # undirected: keep the confident direction
    iu = np.triu_indices(n, k=1)
    vals = sym[iu]
    mask = vals < pae_cutoff
    rows, cols = iu[0][mask], iu[1][mask]
    weights = 1.0 / np.clip(vals[mask], 0.2, None) ** pae_power

    g = nx.Graph()
    g.add_nodes_from(range(n))
    g.add_weighted_edges_from(zip(rows.tolist(), cols.tolist(), weights.tolist()))
    communities = nx.community.louvain_communities(
        g, weight="weight", resolution=resolution, seed=seed)

    labels = [-1] * n
    for cid, members in enumerate(communities):
        for r in members:
            labels[r] = cid
    return labels


def contiguous_segments(labels):
    """Split a per-residue label list into contiguous (start, end, label) runs,
    0-based inclusive."""
    segs = []
    start = 0
    for i in range(1, len(labels) + 1):
        if i == len(labels) or labels[i] != labels[start]:
            segs.append((start, i - 1, int(labels[start])))
            start = i
    return segs


def build_fragments(segs, n, min_domain_len=40, ndr_mode="exclude",
                    min_ndr_len=30, pad_len=10):
    """Classify contiguous segments into domain/ndr fragments and apply the NDR
    policy. Domains are contiguous same-cluster runs >= min_domain_len; every
    residue not in a domain is NDR (coalesced into maximal runs). Returns dicts
    with 0-based inclusive start/end and a 'type' of 'domain' or 'ndr'.

    ndr_mode: 'exclude' (drop NDRs), 'keep' (emit NDR runs >= min_ndr_len as
    their own fragments), or 'pad' (extend domains by pad_len into flanking
    linker, drop the rest)."""
    is_dom = [False] * n
    domains = []
    for (s, e, _lab) in segs:
        if e - s + 1 >= min_domain_len:
            for r in range(s, e + 1):
                is_dom[r] = True
            domains.append(dict(start=s, end=e, length=e - s + 1, type="domain"))

    if ndr_mode == "pad":
        for i, d in enumerate(domains):
            left_limit = domains[i - 1]["end"] + 1 if i > 0 else 0
            right_limit = domains[i + 1]["start"] - 1 if i < len(domains) - 1 else n - 1
            d["start"] = max(left_limit, d["start"] - pad_len)
            d["end"] = min(right_limit, d["end"] + pad_len)
            d["length"] = d["end"] - d["start"] + 1
        return domains

    out = list(domains)
    if ndr_mode == "keep":
        i = 0
        while i < n:
            if is_dom[i]:
                i += 1
                continue
            j = i
            while j < n and not is_dom[j]:
                j += 1
            if j - i >= min_ndr_len:
                out.append(dict(start=i, end=j - 1, length=j - i, type="ndr"))
            i = j
    out.sort(key=lambda f: f["start"])
    return out


def windowize(frag, max_len, overlap):
    """Split an over-long fragment into overlapping windows (0-based inclusive)."""
    if not max_len or frag["length"] <= max_len:
        return [frag]
    pieces, step = [], max_len - overlap
    s = frag["start"]
    while s <= frag["end"]:
        e = min(s + max_len - 1, frag["end"])
        pieces.append(dict(start=s, end=e, length=e - s + 1, type=frag["type"]))
        if e == frag["end"]:
            break
        s += step
    return pieces


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def split_domains_by_afdb(fasta_file, out_dir="domains", cache_dir=None,
                          ndr_mode="exclude", min_domain_len=40, min_ndr_len=30,
                          pad_len=10, max_domain_len=0, overlap=50,
                          pae_power=1.0, pae_cutoff=5.0, resolution=0.3,
                          delay=0.1, limit=0):
    """Domain-split every protein in a FASTA using AlphaFold DB PAE and write a
    sliced FASTA plus a mapping TSV. Proteins whose AFDB model does not cover the
    full sequence (long, multi-fragment entries) are logged, not mis-sliced."""
    _require_scicomp()
    cache_dir = cache_dir or os.path.join(out_dir, "_pae_cache")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(fasta_file))[0]
    fa_out = os.path.join(out_dir, f"{base}_domains.fasta")
    tsv_out = os.path.join(out_dir, f"{base}_domains.tsv")

    proteins = list(iter_fasta(fasta_file))
    if limit:
        proteins = proteins[:limit]

    fa, rows, missing, multi = [], [], [], []
    n_frag = 0
    for i, (acc, gene, seq) in enumerate(proteins, 1):
        pae = fetch_pae(acc, cache_dir)
        if pae is None:
            missing.append(acc)
        elif pae.shape[0] != len(seq):
            multi.append(f"{acc}\tpae_len={pae.shape[0]}\tseq_len={len(seq)}")
        else:
            labels = cluster_pae(pae, pae_power, pae_cutoff, resolution)
            frags = build_fragments(contiguous_segments(labels), len(seq),
                                    min_domain_len, ndr_mode, min_ndr_len, pad_len)
            frags = [w for f in frags for w in windowize(f, max_domain_len, overlap)]
            dom_idx = 0
            for f in frags:
                s1, e1 = f["start"] + 1, f["end"] + 1
                if f["type"] == "domain":
                    dom_idx += 1
                    tag = f"d{dom_idx}"
                else:
                    tag = "ndr"
                frag_id = f"{acc}_{gene}_{tag}_{s1}-{e1}"
                fa.append(f">{frag_id}\n{seq[f['start']:f['end'] + 1]}")
                rows.append(dict(fragment_id=frag_id, accession=acc, gene=gene,
                                 type=f["type"], start=s1, end=e1, length=f["length"]))
                n_frag += 1
        if i % 25 == 0 or i == len(proteins):
            print(f"  {i}/{len(proteins)} proteins | {n_frag} fragments "
                  f"| no_afdb={len(missing)} multifragment={len(multi)}", flush=True)
        time.sleep(delay)

    with open(fa_out, "w") as fh:
        fh.write("\n".join(fa) + ("\n" if fa else ""))
    with open(tsv_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["fragment_id", "accession", "gene",
                                                "type", "start", "end", "length"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    if missing:
        with open(os.path.join(out_dir, "_no_afdb.txt"), "w") as fh:
            fh.write("\n".join(missing) + "\n")
    if multi:
        with open(os.path.join(out_dir, "_multifragment.txt"), "w") as fh:
            fh.write("\n".join(multi) + "\n")

    n_dom = sum(1 for r in rows if r["type"] == "domain")
    n_ndr = sum(1 for r in rows if r["type"] == "ndr")
    print(f"domain-split {len(proteins)} proteins -> {n_frag} fragments "
          f"(domains={n_dom}, ndr={n_ndr}); fasta {fa_out}")
    if missing:
        print(f"  no AFDB model: {len(missing)} (see {out_dir}/_no_afdb.txt)")
    if multi:
        print(f"  multi-fragment, needs separate handling: {len(multi)} "
              f"(see {out_dir}/_multifragment.txt)")
    return fa_out
