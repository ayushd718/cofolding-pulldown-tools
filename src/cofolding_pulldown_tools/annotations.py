"""Split proteins into domain fragments using UniProt sequence annotation.

A counterpart to `domains.py`, which derives domains from the AlphaFold DB PAE
matrix. That approach cannot touch the proteins that most need slicing:
AlphaFold DB has no model above roughly 2700 residues, so FBN1 (2871 aa),
HSPG2 (4391 aa) and MUC16 (14507 aa) return nothing at all, while LAMA5, VWF
and HMCN1 return only short isoform models covering a few hundred residues.
UniProt domain annotation needs no structure model and covers these well: all
37 entries over 2700 aa in a matrisome test set carry annotation, 1174 domain
and repeat features in total.

Fragments are packed rather than cut. Consecutive whole domains accumulate
until the next would exceed the cap, so a boundary always falls between one
domain's end and the next one's start. There is deliberately no attempt to find
a safe residue inside a linker, because there usually is no linker: the median
inter-domain gap is 9 residues, only 46.5% exceed 10, and inside repeat arrays
it collapses to 1-5 residues, with tenascin-C's EGF-like domains overlapping
outright.

Tandem-repeat arrays are grouped by family plus qualifier and collapsed only
above a sequence-identity threshold, and collapsing is opt-in. Measured over
the 37 matrisome proteins above 2700 aa, collapsing at 0.9 removes 4634
residues, 3.1% of the output, but raises the fragment count from 193 to 235
because a dropped copy splits its run in two. More fold jobs for slightly less
sequence is the wrong trade for most pipelines, so it is off by default. That threshold matters more than it looks.
Classical domain arrays are homologous folds with divergent sequence --
measured median pairwise identity is 33% across FBN1's 43 EGF-like copies and
23% across HMCN1's 44 Ig-like copies -- so at the 0.9 default nothing in either
collapses, which is correct rather than a failure. Only mucin-type arrays
(MUC16's SEA domains, median 72%) are genuinely redundant. The large reduction
in fragment count comes from packing, not deduplication.

Pure stdlib plus `requests`; no optional dependencies.
"""
import csv
import json
import os
import re
import time

import requests

from .fasta import iter_fasta
from .domains import windowize

UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{acc}.json"

# JSON, not TSV: in the TSV flattening "; " is both the inter-feature separator
# and appears inside /note="EGF-like 4; calcium-binding", so naive splitting
# destroys exactly the qualifier the dedup key depends on.
FEATURE_FIELDS = ",".join([
    "accession", "gene_primary", "sequence", "ft_domain", "ft_repeat",
    "ft_region", "ft_coiled", "ft_signal", "ft_propep", "ft_disulfid",
])

# Qualifiers that mark a structurally different copy and must never share a
# representative with the plain form. In FBN1, EGF-like 1/2/3/6 lack
# "calcium-binding" while 4/5/7-47 carry it, and cbEGF is a distinct subtype.
QUALIFIERS_DISTINCT = (
    "calcium-binding", "approximate", "truncated", "degenerate",
    "incomplete", "atypical", "extra domain a", "extra domain b",
)

# "<N> X [<M> AA ]['approximate '][tandem ]repeats[ of <consensus>]"
REGION_REPEAT = re.compile(
    r"^(?P<n>\d+)\s*X\s*(?:(?P<unit>\d+)\s*AA\s*)?(?:approximate\s+)?"
    r"(?:tandem\s+)?repeats?(?:\s+of\s+(?P<consensus>.*))?$", re.IGNORECASE)

# Region descriptions that denote a real structural unit. "Disordered" is
# deliberately absent: it is predictor-derived and spans COL1A1's rigid triple
# helix in full, so it must not gate what counts as foldable.
STRUCTURAL_REGIONS = ("triple-helical", "nonhelical", "collagen-like")


# ---------------------------------------------------------------------------
# UniProt access
# ---------------------------------------------------------------------------
def _get(url, params=None, retries=3):
    for attempt in range(retries):
        try:
            resp = requests.get(url, params=params, timeout=60)
        except requests.RequestException:
            time.sleep(2 ** attempt)
            continue
        if resp.status_code == 200:
            return resp.text
        if resp.status_code in (400, 404):
            return None
        time.sleep(2 ** attempt)
    return None


def fetch_annotation(acc, cache_dir, retries=3):
    """Return the UniProt entry payload for an accession, or None if absent.

    Cached to disk so reruns are offline, mirroring `domains.fetch_pae`."""
    cache = os.path.join(cache_dir, f"{acc}-features.json")
    if os.path.exists(cache):
        with open(cache) as fh:
            return json.load(fh)
    text = _get(UNIPROT_ENTRY.format(acc=acc), {"fields": FEATURE_FIELDS},
                retries)
    if text is None:
        return None
    payload = json.loads(text)
    if not payload.get("sequence", {}).get("value"):
        return None
    # Create the directory here too: losing a successful network fetch because
    # the cache directory does not exist is a poor failure mode for a helper
    # that may be called outside the driver.
    os.makedirs(cache_dir, exist_ok=True)
    with open(cache, "w") as fh:
        fh.write(text)
    return payload


# ---------------------------------------------------------------------------
# Description parsing (pure helpers, so they are unit-testable)
# ---------------------------------------------------------------------------
def parse_description(desc):
    """Split a Domain/Repeat description into (family, copy_index, qualifier).

    `copy_index` is None for a singleton domain. `family` comes back empty when
    the description carries nothing but an index -- MUC16's repeats are
    literally "1".."12", as are MUC17's, MUC1's, MUC3A's, hornerin's, MUC12's
    and MUC21's, 365 features in all -- and must then be recovered from the
    enclosing Region. Aggrecan's hierarchical "2-34" form keeps its group in
    the family so its two repeat regions do not merge."""
    desc = (desc or "").strip()
    qual = ""
    if ";" in desc:
        head, _, tail = desc.partition(";")
        desc, qual = head.strip(), tail.strip()

    # A bare index, optionally hierarchical: "1", "1-1", "2-34".
    bare = re.match(r"^(?P<grp>\d+)(?:-(?P<idx>\d+))?$", desc)
    if bare:
        if bare.group("idx"):
            return "#" + bare.group("grp"), int(bare.group("idx")), qual
        return "", int(bare.group("grp")), qual

    # Whitespace before the index is REQUIRED in the named forms, so a digit
    # welded to the family name stays part of it: "Fibrillar collagen NC1" is
    # one singleton domain, not copy 1 of a "Fibrillar collagen NC" family.
    hier = re.match(r"^(?P<fam>.+?)\s+(?P<grp>\d+)-(?P<idx>\d+)$", desc)
    if hier:
        return hier.group("fam").strip(), int(hier.group("idx")), qual

    flat = re.match(r"^(?P<fam>.+?)\s+(?P<idx>\d+)$", desc)
    if flat:
        return flat.group("fam").strip(), int(flat.group("idx")), qual

    return desc, None, qual


def qualifier_class(qual):
    """Normalise a qualifier into a dedup class.

    "first part"/"second part" are an annotation artifact describing one
    discontinuous domain, so they collapse to the base class and the two
    features get merged elsewhere. Everything in QUALIFIERS_DISTINCT marks a
    structurally different copy and is kept apart."""
    q = (qual or "").strip().lower()
    if q in ("first part", "second part"):
        return ""
    for keep in QUALIFIERS_DISTINCT:
        if keep in q:
            return keep
    return ""


def parse_repeat_region(desc):
    """Parse a Region repeat description into (n_copies, unit_len, consensus).

    Returns None when the text is not a repeat annotation. This is the only
    source of repeat period and consensus for the low-complexity proteins,
    covering 49 features across 33 of 1027 matrisome entries, e.g. MUC1's
    "42 X 20 AA approximate tandem repeats of P-A-P-G-S-T-...". """
    m = REGION_REPEAT.match((desc or "").strip())
    if not m:
        return None
    unit = int(m.group("unit")) if m.group("unit") else None
    return int(m.group("n")), unit, (m.group("consensus") or "").strip()


# ---------------------------------------------------------------------------
# Block construction
# ---------------------------------------------------------------------------
def _gene_of(payload, fallback):
    """Prefer UniProt's primary gene symbol over the FASTA header token.

    `iter_fasta` returns the whole token as the gene for a cleaned header, so
    `FBN1_HUMAN` would land in the fragment id and its extra underscore would
    break the `accession_gene__ligand` job-name split downstream. The UniProt
    symbol is also the real one: the entry mnemonic diverges from the gene
    symbol often enough that `reformat_sliced_windows` exists to repair it
    (`CAD23` for `CDH23`)."""
    genes = payload.get("genes") or []
    if genes:
        name = (genes[0].get("geneName") or {}).get("value")
        if name:
            return name
    return fallback


def split_identifier(ident):
    """Recover (accession, label) from a FASTA identifier.

    `iter_fasta` yields the bare accession for a raw UniProt header, but a
    `clean_fasta`d header is `ACC_MNEMONIC` (e.g. `P35555_FBN1_HUMAN`) and
    comes back whole, which would be sent to UniProt as-is and 404. An
    accession may contain a hyphen for an isoform (`Q99102-10`) but never an
    underscore, so the first underscore is the split point."""
    if "_" in ident:
        acc, _, rest = ident.partition("_")
        return acc, rest or acc
    return ident, ident


def _exact(loc):
    return (loc.get("start", {}).get("modifier", "EXACT") == "EXACT"
            and loc.get("end", {}).get("modifier", "EXACT") == "EXACT")


def collect_features(payload, min_repeat_len=40):
    """Turn an entry payload into (blocks, regions, bonds, sequence).

    Coordinates are converted to 0-based inclusive on ingest, matching
    `domains.py`, and converted back on output. Repeat features shorter than
    `min_repeat_len` are held back as solenoid sub-units: SLIT2's 20
    leucine-rich repeats are 22 residues each but form a single fold, and PRG4
    has 59 repeats of 8 residues.

    Signal peptides and propeptides come back separately and are never units.
    A signal peptide is cleaved from the mature protein, so including it would
    add a hydrophobic helix that nucleates spurious interfaces in cofolding."""
    seq = payload.get("sequence", {}).get("value", "")
    blocks, short_repeats, regions, bonds, cleaved = [], [], [], [], []
    for f in payload.get("features", []):
        loc = f.get("location", {})
        if not _exact(loc):
            continue
        s = loc.get("start", {}).get("value")
        e = loc.get("end", {}).get("value")
        if s is None or e is None or e < s:
            continue
        s0, e0 = s - 1, e - 1
        ftype = f.get("type")
        desc = f.get("description", "") or ""
        if ftype == "Disulfide bond":
            bonds.append((s0, e0))
            continue
        if ftype in ("Signal", "Propeptide"):
            cleaved.append(dict(start=s0, end=e0, type=ftype.lower()))
            continue
        if ftype == "Region":
            regions.append(dict(start=s0, end=e0, description=desc))
            continue
        if ftype == "Coiled coil":
            regions.append(dict(start=s0, end=e0, description="Coiled coil"))
            continue
        if ftype not in ("Domain", "Repeat"):
            continue
        fam, idx, qual = parse_description(desc)
        block = dict(start=s0, end=e0, length=e0 - s0 + 1, type="domain",
                     family=fam, index=idx, qualifier=qual,
                     qclass=qualifier_class(qual), names=[desc])
        if ftype == "Repeat" and block["length"] < min_repeat_len:
            short_repeats.append(block)
        else:
            blocks.append(block)
    return blocks, short_repeats, regions, bonds, cleaved, seq


def aggregate_short_repeats(short_repeats, gap=None):
    """Coalesce runs of sub-fold repeat features into single blocks.

    A 22-residue leucine-rich repeat is not a foldable fragment; the solenoid
    the repeats build is. Consecutive features of the same family within `gap`
    residues merge into one block."""
    if not short_repeats:
        return []
    out = []
    for b in sorted(short_repeats, key=lambda x: x["start"]):
        limit = gap if gap is not None else max(b["length"] * 2, 20)
        if out:
            prev = out[-1]
            same = prev["family"] == b["family"]
            near = b["start"] - prev["end"] - 1 <= limit
            if same and near:
                prev["end"] = max(prev["end"], b["end"])
                prev["length"] = prev["end"] - prev["start"] + 1
                prev["names"].append(b["names"][0])
                prev["index"] = None
                continue
        merged = dict(b)
        merged["type"] = "repeat_array"
        merged["names"] = list(b["names"])
        out.append(merged)
    return out


def merge_split_domains(blocks):
    """Merge "first part"/"second part" feature pairs into one block.

    Perlecan annotates Laminin EGF-like 1 as a 10-residue first part at
    521-530 and a second part at 731-763, interrupted by Laminin IV type A 1.
    Treated independently the first part becomes a useless 10-residue
    fragment; merged, the span is one unit that must not be split."""
    by_key = {}
    for b in blocks:
        q = (b["qualifier"] or "").strip().lower()
        if q in ("first part", "second part"):
            by_key.setdefault((b["family"], b["index"]), []).append(b)
    merged_ids = set()
    out = []
    for key, parts in by_key.items():
        if len(parts) < 2:
            continue
        parts.sort(key=lambda x: x["start"])
        head = dict(parts[0])
        head["start"] = min(p["start"] for p in parts)
        head["end"] = max(p["end"] for p in parts)
        head["length"] = head["end"] - head["start"] + 1
        head["qualifier"] = ""
        head["qclass"] = ""
        head["names"] = [n for p in parts for n in p["names"]]
        out.append(head)
        merged_ids.update(id(p) for p in parts)
    out.extend(b for b in blocks if id(b) not in merged_ids)
    return sorted(out, key=lambda x: (x["start"], x["end"]))


def recover_families(blocks, regions):
    """Name bare-index repeat blocks from the enclosing Region description.

    Without this the 365 bare-index features collapse into one nameless
    family, so MUC16's 12 SEA repeats and aggrecan's two distinct repeat
    regions would be deduplicated against each other."""
    repeat_regions = [r for r in regions if parse_repeat_region(r["description"])]
    for b in blocks:
        if b["family"]:
            continue
        host = next((r for r in repeat_regions
                     if r["start"] <= b["start"] and b["end"] <= r["end"]), None)
        b["family"] = host["description"] if host else "unnamed repeat"
    return blocks


def disjoin_blocks(blocks):
    """Make blocks disjoint without merging neighbours.

    Nested blocks are absorbed into their container; a partial overlap has the
    later block's start pushed past the earlier block's end. Merging on any
    overlap would be wrong: tenascin-C's 15 EGF-like domains each overlap the
    next by one residue, and merging would fuse the whole array into a single
    block."""
    out = []
    for b in sorted(blocks, key=lambda x: (x["start"], -x["length"])):
        if not out:
            out.append(dict(b))
            continue
        prev = out[-1]
        if b["end"] <= prev["end"]:
            prev["names"].extend(b["names"])
            continue
        nb = dict(b)
        if nb["start"] <= prev["end"]:
            nb["start"] = prev["end"] + 1
            nb["length"] = nb["end"] - nb["start"] + 1
        if nb["length"] > 0:
            out.append(nb)
    return out


# ---------------------------------------------------------------------------
# Repeat similarity and clustering
# ---------------------------------------------------------------------------
def shift_identity(a, b, max_shift=15, min_overlap=15):
    """Best ungapped identity over a sweep of shifts, normalised by overlap.

    Repeat copies differ by one- and two-residue indels -- MUC16's SEA copies
    range 155-158 residues -- so a single fixed alignment underestimates
    identity and a fixed period goes out of phase after two copies. Sweeping
    the shift recovers the real figure.

    k-mer Jaccard was tried first and discarded: on FBN1's EGF-like array,
    which is genuinely 33% identical, 3-mer Jaccard returns a median of 0.03
    and gives zero discrimination at any threshold."""
    if not a or not b:
        return 0.0
    best = 0.0
    for d in range(-max_shift, max_shift + 1):
        x, y = (a[d:], b) if d >= 0 else (a, b[-d:])
        n = min(len(x), len(y))
        if n < min_overlap:
            continue
        same = sum(1 for i in range(n) if x[i] == y[i])
        if same / n > best:
            best = same / n
    return best


def cluster_copies(copies, threshold=0.9, max_shift=15):
    """Greedily cluster (key, sequence) copies by identity, preserving order.

    Returns {key: (representative_key, identity_to_representative)}. The first
    and last copy of an array are always their own representatives, because
    terminal copies diverge from the interior.

    At the 0.9 default most arrays collapse nothing. That is the measured and
    intended behaviour: FBN1's 43 calcium-binding EGF-like copies have zero
    pairs at or above 80% identity, so folding each separately is legitimate
    rather than wasteful."""
    assign = {}
    reps = []
    forced = set()
    if copies:
        forced = {copies[0][0], copies[-1][0]}
    for key, seq in copies:
        if key in forced:
            reps.append((key, seq))
            assign[key] = (key, 1.0)
            continue
        hit = None
        for rkey, rseq in reps:
            score = shift_identity(seq, rseq, max_shift)
            if score >= threshold:
                hit = (rkey, score)
                break
        if hit is None:
            reps.append((key, seq))
            assign[key] = (key, 1.0)
        else:
            assign[key] = hit
    return assign


def collapse_repeat_blocks(blocks, seq, threshold=0.9, max_shift=15):
    """Identify redundant copies within each (family, qualifier class) group.

    Returns (dropped_starts, cluster_rows). Grouping on family alone would be
    wrong: FBN1's plain and calcium-binding EGF-like copies are distinct
    subtypes and must not share a representative."""
    groups = {}
    for b in blocks:
        if b["index"] is None:
            continue
        groups.setdefault((b["family"], b["qclass"]), []).append(b)

    dropped = set()
    rows = []
    for (family, qclass), members in groups.items():
        if len(members) < 2:
            continue
        members.sort(key=lambda x: x["start"])
        copies = [(b["start"], seq[b["start"]:b["end"] + 1]) for b in members]
        assign = cluster_copies(copies, threshold, max_shift)
        for b in members:
            rep_start, score = assign[b["start"]]
            kept = rep_start == b["start"]
            if not kept:
                dropped.add(b["start"])
            rows.append(dict(family=family, qualifier=qclass,
                             copy_index=b["index"],
                             start=b["start"] + 1, end=b["end"] + 1,
                             representative_start=rep_start + 1,
                             identity=round(score, 3),
                             kept="yes" if kept else "no"))
    return dropped, rows


# ---------------------------------------------------------------------------
# Packing and unannotated regions
# ---------------------------------------------------------------------------
def contiguous_groups(blocks, dropped_starts):
    """Split the ordered block list wherever a dropped repeat copy intervenes.

    Without this, deduplication saves nothing. A packed run spans its first
    block's start to its last block's end, so removing copy 7 of a 16-copy
    array leaves its residues sitting inside a fragment that still spans copies
    1 to 10. Measured before this fix: collapsing 59 copies across the 37 long
    matrisome proteins made the emitted output 63 residues LONGER, not shorter.

    Breaking the list at each dropped copy means a dropped copy genuinely
    shortens the output, while every emitted fragment stays one contiguous
    slice of the parent protein, so the START-END in its header remains
    truthful. Fragments are never built by splicing non-adjacent sequence."""
    groups, cur = [], []
    for b in blocks:
        if b["start"] in dropped_starts:
            if cur:
                groups.append(cur)
                cur = []
            continue
        cur.append(b)
    if cur:
        groups.append(cur)
    return groups


def pack_blocks(blocks, max_len):
    """Group ordered disjoint blocks into runs whose span fits `max_len`.

    A run spans its first block's start to its last block's end, so every
    boundary falls between one block's end and the next one's start. A single
    block longer than `max_len` becomes a run of its own and is tiled by the
    caller -- the only case where a fold gets split."""
    if max_len <= 0:
        raise ValueError("max_len must be greater than 0")
    runs, cur = [], []
    for b in blocks:
        if not cur:
            cur = [b]
            continue
        if b["end"] - cur[0]["start"] + 1 <= max_len:
            cur.append(b)
        else:
            runs.append(cur)
            cur = [b]
    if cur:
        runs.append(cur)
    return runs


def respect_disulfides(runs, bonds, max_len):
    """Shift a run boundary by one block when a disulfide bond straddles it.

    FBN1 carries 155 disulfide bonds and HMCN1 76; separating the two
    cysteines of a pair leaves both fragments unable to form the bond.
    Returns (runs, n_unresolved)."""
    if not bonds or len(runs) < 2:
        return runs, 0
    runs = [list(r) for r in runs]
    unresolved = 0
    i = 0
    while i < len(runs) - 1:
        left, right = runs[i], runs[i + 1]
        edge_l, edge_r = left[-1]["end"], right[0]["start"]
        straddling = [(a, b) for a, b in bonds if a <= edge_l < edge_r <= b]
        if not straddling:
            i += 1
            continue
        moved = False
        if len(right) > 1:
            span = right[0]["end"] - left[0]["start"] + 1
            if span <= max_len:
                left.append(right.pop(0))
                moved = True
        if not moved and len(left) > 1:
            span = right[-1]["end"] - left[-1]["start"] + 1
            if span <= max_len:
                right.insert(0, left.pop())
                moved = True
        if not moved:
            unresolved += len(straddling)
            i += 1
    return [r for r in runs if r], unresolved


def find_period(seq, min_period=20, max_period=400, min_identity=0.5):
    """Best tandem period in a sequence by offset-identity scan.

    Returns (period, identity), with period None when nothing passes
    `min_identity`. This correctly finds MUC16's annotated SEA array at a
    period of 156 and 61% identity, and correctly finds nothing in MUC16's
    12066-residue N-terminal region, where the best offset identity is 15.8%,
    barely above composition background. Exact k-mer search was rejected: no
    25-mer recurs even once there, and in the COL1A1 triple helix only the
    Gly-X-Y register is detectable, which is not a fragment boundary.

    Cost is about 4.6M character comparisons on a 12kaa sequence."""
    n = len(seq)
    best_d, best_f = None, 0.0
    for d in range(min_period, min(max_period, n // 2) + 1):
        m = n - d
        if m < min_period:
            break
        same = sum(1 for i in range(m) if seq[i] == seq[i + d])
        frac = same / m
        if frac > best_f:
            best_d, best_f = d, frac
    if best_f >= min_identity:
        return best_d, best_f
    return None, best_f


def unannotated_runs(spans, n, min_len):
    """Maximal 0-based inclusive spans not covered by any input span.

    Call this with the EMITTED fragment spans, not the raw annotation blocks.
    A packed run covers its interior gaps, so measuring against blocks would
    re-emit an enclosed linker that is already inside a domain fragment --
    FBN1's proline-rich 390-448 gap sits inside the 81-1445 fragment."""
    covered = [False] * n
    for b in spans:
        for r in range(max(0, b["start"]), min(n, b["end"] + 1)):
            covered[r] = True
    runs, i = [], 0
    while i < n:
        if covered[i]:
            i += 1
            continue
        j = i
        while j < n and not covered[j]:
            j += 1
        if j - i >= min_len:
            runs.append(dict(start=i, end=j - 1, length=j - i, type="ndr"))
        i = j
    return runs


def subtract_covered(start, end, blocks, min_len):
    """Return the parts of [start, end] that no block already covers.

    0-based inclusive on both input and output."""
    pieces = []
    cursor = start
    for b in sorted(blocks, key=lambda x: x["start"]):
        if b["end"] < start or b["start"] > end:
            continue
        if b["start"] > cursor:
            lo, hi = cursor, min(b["start"] - 1, end)
            if hi - lo + 1 >= min_len:
                pieces.append((lo, hi))
        cursor = max(cursor, b["end"] + 1)
        if cursor > end:
            break
    if cursor <= end and end - cursor + 1 >= min_len:
        pieces.append((cursor, end))
    return pieces


def structural_region_spans(regions, covered, min_len=40):
    """Cover residues no domain covers, using structurally meaningful Regions.

    The span is added by SUBTRACTION, never by merging. COL6A3 carries a
    2013-residue "Nonhelical region" spanning 19 of its domains; added whole it
    would absorb every one of them and leave a single uncuttable block.
    Subtracted against what the domains already cover, it contributes only the
    genuinely uncovered leftovers and the 19 domains survive.

    "Disordered" is excluded on purpose: it is predictor-derived and covers
    COL1A1's rigid triple helix in full, so trusting it would discard a folded
    assembly as unstructured. The real unit there is "Triple-helical region",
    and matching is a substring test because qualified forms exist such as
    "Triple-helical region (COL2) with 1 imperfection"."""
    out = []
    for r in regions:
        d = (r["description"] or "").lower()
        if not any(k in d for k in STRUCTURAL_REGIONS):
            continue
        for lo, hi in subtract_covered(r["start"], r["end"], covered, min_len):
            out.append(dict(start=lo, end=hi, length=hi - lo + 1,
                            type="region", family=r["description"],
                            index=None, qualifier="", qclass="",
                            names=[r["description"]]))
    return out


def handle_ndr(runs, seq, mode, max_len, overlap, min_identity=0.5):
    """Apply the non-domain-residue policy, returning emittable fragments.

    `exclude`, `keep` and `pad` mirror the AlphaFold path. `tile` additionally
    runs a period scan first and, on a hit, emits one representative unit plus
    its flanks instead of the whole stretch."""
    if mode == "exclude":
        return []
    out = []
    for run in runs:
        if mode in ("keep", "pad"):
            out.extend(windowize(dict(run), max_len, overlap))
            continue
        sub = seq[run["start"]:run["end"] + 1]
        period, score = find_period(sub, min_identity=min_identity)
        if period and run["length"] > 2 * period:
            unit = dict(start=run["start"], end=run["start"] + period - 1,
                        length=period, type="ndr_repeat")
            out.append(unit)
            tail_start = run["end"] - period + 1
            if tail_start > unit["end"]:
                out.append(dict(start=tail_start, end=run["end"],
                                length=run["end"] - tail_start + 1,
                                type="ndr_repeat"))
        else:
            out.extend(windowize(dict(run), max_len, overlap))
    return out


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def split_domains_by_annotation(fasta_file, out_dir="uniprot_domains",
                                max_slice=1400, cache_dir=None, ndr_mode="tile",
                                min_domain_len=40, min_ndr_len=30, pad_len=10,
                                overlap=50, repeat_identity=0.9, big_ndr=1000,
                                collapse_repeats=False, delay=0.1, limit=0):
    """Annotation-split every protein in a FASTA and write a sliced FASTA plus
    two mapping tables. Proteins with no usable annotation are logged, not
    mis-sliced."""
    if max_slice <= 0:
        raise ValueError("max_slice must be greater than 0")
    if not 0.0 < repeat_identity <= 1.0:
        raise ValueError("repeat_identity must be in (0, 1]")

    cache_dir = cache_dir or os.path.join(out_dir, "_annotation_cache")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)
    base = os.path.splitext(os.path.basename(fasta_file))[0]
    fa_out = os.path.join(out_dir, f"{base}_uniprot_domains.fasta")
    tsv_out = os.path.join(out_dir, f"{base}_uniprot_domains.tsv")
    clu_out = os.path.join(out_dir, f"{base}_repeat_clusters.tsv")

    proteins = list(iter_fasta(fasta_file))
    if limit:
        proteins = proteins[:limit]

    fa, rows, cluster_rows = [], [], []
    missing, bare, isoforms, big_runs = [], [], [], []
    n_frag = n_unresolved = 0
    for i, (ident, gene, seq) in enumerate(proteins, 1):
        acc, label = split_identifier(ident)
        if gene == ident:
            gene = label
        if "-" in acc:
            # Isoform accessions carry a sequence but no features at all, so
            # there is nothing to slice on. Canonical coordinates do not apply
            # to an isoform sequence, so borrowing them would shift every
            # boundary silently.
            isoforms.append(f"{acc}\tlen={len(seq)}")
            payload = None
            missing_handled = True
        else:
            payload = fetch_annotation(acc, cache_dir)
            missing_handled = False
        if payload is None:
            if not missing_handled:
                missing.append(acc)
        else:
            blocks, shorts, regions, bonds, cleaved, ref = collect_features(
                payload, min_domain_len)
            seq = ref or seq
            gene = _gene_of(payload, gene)
            blocks = merge_split_domains(blocks)
            blocks.extend(aggregate_short_repeats(shorts))
            blocks = recover_families(blocks, regions)
            blocks = disjoin_blocks(blocks)
            blocks.extend(structural_region_spans(regions, blocks,
                                                  min_domain_len))
            blocks = disjoin_blocks(blocks)
            if not blocks:
                bare.append(f"{acc}\tlen={len(seq)}")

            dropped = set()
            if collapse_repeats and blocks:
                dropped, crows = collapse_repeat_blocks(
                    blocks, seq, repeat_identity)
                for r in crows:
                    r.update(accession=acc, gene=gene)
                cluster_rows.extend(crows)

            runs = []
            for group in contiguous_groups(blocks, dropped):
                grun = pack_blocks(group, max_slice)
                grun, bad = respect_disulfides(grun, bonds, max_slice)
                n_unresolved += bad
                runs.extend(grun)

            frags = []
            for run in runs:
                span = dict(start=run[0]["start"], end=run[-1]["end"],
                            length=run[-1]["end"] - run[0]["start"] + 1,
                            type="domain")
                pieces = windowize(span, max_slice, overlap)
                for p in pieces:
                    p["names"] = [n for b in run for n in b["names"]]
                    p["family"] = run[0]["family"]
                    p["qualifier"] = run[0]["qclass"]
                frags.extend(pieces)

            # A collapsed repeat copy leaves a hole that is not "unannotated"
            # -- it was deliberately dropped. Without excluding it here the NDR
            # policy tiles it straight back in, and deduplication ends up
            # making the output longer rather than shorter.
            suppressed = list(cleaved) + [b for b in blocks
                                          if b["start"] in dropped]
            raw_ndr = unannotated_runs(frags, len(seq), min_ndr_len)
            trimmed = []
            for run in raw_ndr:
                for lo, hi in subtract_covered(run["start"], run["end"],
                                               suppressed, min_ndr_len):
                    trimmed.append(dict(start=lo, end=hi,
                                        length=hi - lo + 1, type="ndr"))
            ndr = handle_ndr(trimmed, seq, ndr_mode, max_slice, overlap)
            # A large aperiodic stretch tiled blindly is the expensive case and
            # should never be silent. MUC16 alone contributes 12066 such
            # residues, 83% of the protein, which become 9 fragments of
            # glycan-covered mucin sequence under the default policy.
            for run in trimmed:
                if run["length"] >= big_ndr and ndr_mode != "exclude":
                    pct = 100.0 * run["length"] / max(len(seq), 1)
                    print(f"  WARN {acc} {gene}: {run['length']} residues "
                          f"({pct:.1f}%) unannotated and aperiodic at "
                          f"{run['start'] + 1}-{run['end'] + 1}; tiled under "
                          f"--ndr {ndr_mode}", flush=True)
                    big_runs.append(f"{acc}\t{gene}\t{run['start'] + 1}\t"
                                    f"{run['end'] + 1}\t{run['length']}")
            for p in ndr:
                p.setdefault("names", [])
                p.setdefault("family", "")
                p.setdefault("qualifier", "")
            frags.extend(ndr)
            frags.sort(key=lambda f: f["start"])

            counters = {}
            for f in frags:
                kind = f["type"]
                counters[kind] = counters.get(kind, 0) + 1
                tag = ("d" if kind == "domain" else kind) + str(counters[kind])
                s1, e1 = f["start"] + 1, f["end"] + 1
                frag_id = f"{acc}_{gene}_{tag}_{s1}-{e1}"
                fa.append(f">{frag_id}\n{seq[f['start']:f['end'] + 1]}")
                rows.append(dict(fragment_id=frag_id, accession=acc, gene=gene,
                                 type=kind, family=f.get("family", ""),
                                 qualifier=f.get("qualifier", ""),
                                 start=s1, end=e1, length=f["length"],
                                 n_domains=len(f.get("names", [])),
                                 domains="|".join(f.get("names", []))))
                n_frag += 1
        if i % 25 == 0 or i == len(proteins):
            print(f"  {i}/{len(proteins)} proteins | {n_frag} fragments "
                  f"| no_annotation={len(missing)} isoform={len(isoforms)}",
                  flush=True)
        time.sleep(delay)

    with open(fa_out, "w") as fh:
        fh.write("\n".join(fa) + ("\n" if fa else ""))
    with open(tsv_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "fragment_id", "accession", "gene", "type", "family", "qualifier",
            "start", "end", "length", "n_domains", "domains"])
        writer.writeheader()
        writer.writerows(rows)
    with open(clu_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "accession", "gene", "family", "qualifier", "copy_index", "start",
            "end", "representative_start", "identity", "kept"])
        writer.writeheader()
        writer.writerows(cluster_rows)
    if missing:
        with open(os.path.join(out_dir, "_no_annotation.txt"), "w") as fh:
            fh.write("\n".join(missing) + "\n")
    if bare:
        with open(os.path.join(out_dir, "_no_features.txt"), "w") as fh:
            fh.write("\n".join(bare) + "\n")
    if isoforms:
        with open(os.path.join(out_dir, "_isoforms_skipped.txt"), "w") as fh:
            fh.write("\n".join(isoforms) + "\n")
    if big_runs:
        with open(os.path.join(out_dir, "_large_unannotated.tsv"), "w") as fh:
            fh.write("accession\tgene\tstart\tend\tlength\n")
            fh.write("\n".join(big_runs) + "\n")

    n_dropped = sum(1 for r in cluster_rows if r["kept"] == "no")
    print(f"annotation-split {len(proteins)} proteins -> {n_frag} fragments; "
          f"fasta {fa_out}")
    print(f"  fragment map {tsv_out}")
    print(f"  repeat clusters {clu_out} ({n_dropped} copies collapsed)")
    if missing:
        print(f"  no UniProt annotation: {len(missing)} "
              f"(see {out_dir}/_no_annotation.txt)")
    if bare:
        print(f"  no domain/repeat features: {len(bare)} "
              f"(see {out_dir}/_no_features.txt)")
    if isoforms:
        print(f"  isoform accessions skipped (no annotation exists): "
              f"{len(isoforms)} (see {out_dir}/_isoforms_skipped.txt)")
    if big_runs:
        print(f"  large unannotated stretches tiled: {len(big_runs)} "
              f"(see {out_dir}/_large_unannotated.tsv); consider --ndr exclude")
    if n_unresolved:
        print(f"  disulfide pairs straddling a boundary: {n_unresolved}")
    return fa_out
