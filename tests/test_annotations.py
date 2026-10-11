import pytest

from cofolding_pulldown_tools.annotations import (
    aggregate_short_repeats,
    cluster_copies,
    collapse_repeat_blocks,
    collect_features,
    disjoin_blocks,
    find_period,
    handle_ndr,
    merge_split_domains,
    pack_blocks,
    parse_description,
    parse_repeat_region,
    qualifier_class,
    recover_families,
    respect_disulfides,
    shift_identity,
    split_identifier,
    unannotated_runs,
)


def _block(start, end, family="EGF-like", index=None, qualifier="",
           btype="domain"):
    return dict(start=start, end=end, length=end - start + 1, type=btype,
                family=family, index=index, qualifier=qualifier,
                qclass=qualifier_class(qualifier), names=[f"{family} {index}"])


# ---------------------------------------------------------------------------
# Identifier and description parsing
# ---------------------------------------------------------------------------
def test_split_identifier_handles_cleaned_headers():
    """A clean_fasta'd header is ACC_MNEMONIC; sending it to UniProt whole 404s.
    An isoform accession carries a hyphen but never an underscore."""
    assert split_identifier("P35555_FBN1_HUMAN") == ("P35555", "FBN1_HUMAN")
    assert split_identifier("P35555") == ("P35555", "P35555")
    assert split_identifier("Q99102-10_MUC4_HUMAN") == ("Q99102-10", "MUC4_HUMAN")


def test_parse_description_family_index_qualifier():
    """Tandem copies are named '<FAMILY> <N>', and the qualifier after the
    semicolon distinguishes structurally different copies."""
    assert parse_description("EGF-like 1") == ("EGF-like", 1, "")
    assert parse_description("EGF-like 14; calcium-binding") == (
        "EGF-like", 14, "calcium-binding")
    assert parse_description("Ig-like C2-type 22") == ("Ig-like C2-type", 22, "")
    assert parse_description("TB 1") == ("TB", 1, "")


def test_parse_description_singleton_has_no_index():
    """A domain with no tandem copies carries no index and must be left alone."""
    assert parse_description("VWFA") == ("VWFA", None, "")
    assert parse_description("Fibrillar collagen NC1") == (
        "Fibrillar collagen NC1", None, "")
    assert parse_description("Ig-like C2-type") == ("Ig-like C2-type", None, "")


def test_parse_description_bare_index_yields_empty_family():
    """MUC16's repeat descriptions are literally '1'..'12' (365 such features
    corpus-wide). The family must come back empty so recover_families fills it,
    rather than every mucin repeat collapsing into one nameless group."""
    assert parse_description("1") == ("", 1, "")
    assert parse_description("1; approximate") == ("", 1, "approximate")


def test_parse_description_hierarchical_index_keeps_group():
    """Aggrecan uses a two-level '<region>-<copy>' scheme. Naive stripping gives
    family '1', which would merge its two distinct repeat regions."""
    assert parse_description("1-1") == ("#1", 1, "")
    assert parse_description("2-34") == ("#2", 34, "")


def test_qualifier_class_collapses_split_parts_only():
    """first/second part is an annotation artifact for one discontinuous domain,
    so it must not create a separate dedup group. calcium-binding must, because
    cbEGF and plain EGF are different subtypes."""
    assert qualifier_class("first part") == ""
    assert qualifier_class("second part") == ""
    assert qualifier_class("calcium-binding") == "calcium-binding"
    assert qualifier_class("truncated") == "truncated"
    assert qualifier_class("") == ""


def test_parse_repeat_region_grammar():
    """Region repeat text is the only source of period and consensus for the
    low-complexity proteins."""
    assert parse_repeat_region(
        "42 X 20 AA approximate tandem repeats of P-A-P-G-S-T-A-P") == (
        42, 20, "P-A-P-G-S-T-A-P")
    assert parse_repeat_region("12 X approximate tandem repeats") == (
        12, None, "")
    assert parse_repeat_region("59 X 8 AA repeats of K-X-P-X-P-T-T-X") == (
        59, 8, "K-X-P-X-P-T-T-X")
    assert parse_repeat_region("Disordered") is None
    assert parse_repeat_region("Triple-helical region") is None


# ---------------------------------------------------------------------------
# Block construction
# ---------------------------------------------------------------------------
def test_collect_features_skips_fuzzy_coordinates():
    """A location modifier of UNKNOWN or OUTSIDE means the coordinate is not
    trustworthy and the feature must not drive a cut."""
    payload = dict(sequence=dict(value="A" * 300), features=[
        dict(type="Domain", description="EGF-like 1",
             location=dict(start=dict(value=1, modifier="EXACT"),
                           end=dict(value=60, modifier="EXACT"))),
        dict(type="Domain", description="EGF-like 2",
             location=dict(start=dict(value=70, modifier="UNKNOWN"),
                           end=dict(value=130, modifier="EXACT"))),
    ])
    blocks, shorts, regions, bonds, cleaved, seq = collect_features(payload)
    assert [b["names"][0] for b in blocks] == ["EGF-like 1"]
    assert len(seq) == 300


def test_collect_features_separates_short_repeats():
    """A 22-residue leucine-rich repeat is a solenoid sub-unit, not a foldable
    fragment; SLIT2 has 20 of them forming a single fold."""
    payload = dict(sequence=dict(value="A" * 300), features=[
        dict(type="Repeat", description="LRR 1",
             location=dict(start=dict(value=1, modifier="EXACT"),
                           end=dict(value=22, modifier="EXACT"))),
        dict(type="Domain", description="EGF-like 1",
             location=dict(start=dict(value=100, modifier="EXACT"),
                           end=dict(value=160, modifier="EXACT"))),
    ])
    blocks, shorts, _, _, _, _ = collect_features(payload, min_repeat_len=40)
    assert [b["names"][0] for b in blocks] == ["EGF-like 1"]
    assert [b["names"][0] for b in shorts] == ["LRR 1"]


def test_aggregate_short_repeats_coalesces_a_solenoid():
    """Twenty abutting 22-residue repeats become one block, not twenty
    fragments."""
    shorts = [_block(i * 22, i * 22 + 21, family="LRR", index=i + 1)
              for i in range(20)]
    out = aggregate_short_repeats(shorts)
    assert len(out) == 1
    assert (out[0]["start"], out[0]["end"]) == (0, 439)
    assert out[0]["type"] == "repeat_array"
    assert len(out[0]["names"]) == 20


def test_merge_split_domains_rescues_a_ten_residue_fragment():
    """Perlecan annotates Laminin EGF-like 1 as a 10-residue first part at
    521-530 plus a second part at 731-763, interrupted by another domain.
    Emitted independently the first part is a useless 10-residue fragment."""
    blocks = [
        _block(520, 529, family="Laminin EGF-like", index=1,
               qualifier="first part"),
        _block(537, 729, family="Laminin IV type A", index=1),
        _block(730, 762, family="Laminin EGF-like", index=1,
               qualifier="second part"),
    ]
    out = merge_split_domains(blocks)
    merged = [b for b in out if b["family"] == "Laminin EGF-like"]
    assert len(merged) == 1
    assert (merged[0]["start"], merged[0]["end"]) == (520, 762)
    assert merged[0]["length"] == 243
    assert not any(b["length"] == 10 for b in out)


def test_disjoin_blocks_does_not_fuse_an_overlapping_array():
    """Tenascin-C's EGF-like domains overlap the next by one residue. Merging on
    any overlap would fuse the whole array into a single block, so a partial
    overlap must only push the later start."""
    blocks = [_block(185, 216, index=1), _block(216, 247, index=2),
              _block(247, 278, index=3)]
    out = disjoin_blocks(blocks)
    assert [(b["start"], b["end"]) for b in out] == [
        (185, 216), (217, 247), (248, 278)]


def test_disjoin_blocks_absorbs_a_nested_block():
    """FBN1 nests a Hybrid domain inside TB 1; the container wins and records
    the absorbed name."""
    blocks = [_block(183, 235, family="TB", index=1),
              _block(194, 220, family="Hybrid domain", index=1)]
    out = disjoin_blocks(blocks)
    assert len(out) == 1
    assert (out[0]["start"], out[0]["end"]) == (183, 235)
    assert "Hybrid domain 1" in out[0]["names"]


def test_recover_families_names_bare_repeats_from_region():
    """Without this, MUC16's 12 bare-index repeats and aggrecan's two distinct
    repeat regions would deduplicate against each other."""
    blocks = [_block(100, 250, family="", index=1),
              _block(260, 410, family="", index=2)]
    regions = [dict(start=90, end=500,
                    description="12 X approximate tandem repeats"),
               dict(start=0, end=50, description="Disordered")]
    out = recover_families(blocks, regions)
    assert {b["family"] for b in out} == {"12 X approximate tandem repeats"}


# ---------------------------------------------------------------------------
# Similarity and clustering
# ---------------------------------------------------------------------------
def test_shift_identity_recovers_identity_through_an_indel():
    """Repeat copies differ by one- and two-residue indels (MUC16's SEA copies
    span 155-158 aa), so a fixed alignment underestimates identity."""
    a = "MKTAYIAKQRQISFVKSHFSRQ"
    assert shift_identity(a, a, min_overlap=10) == 1.0
    # same sequence with one residue inserted at the front
    assert shift_identity("X" + a, a, min_overlap=10) == 1.0
    assert shift_identity(a, "QQQQQQQQQQQQQQQQQQQQQQ", min_overlap=10) < 0.3


def test_shift_identity_handles_empty_input():
    assert shift_identity("", "ACDEFGHIKLMNPQRST") == 0.0


def test_cluster_copies_collapses_only_near_duplicates():
    """Mucin-type arrays are genuinely redundant; classical domain arrays are
    not. Two near-identical copies share a representative, a divergent one does
    not."""
    a = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
    b = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKM"      # 1 substitution
    c = "QQWWEERRTTYYUUIIOOPPAASSDDFFGG"      # unrelated
    mid = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
    copies = [(0, a), (100, b), (200, c), (300, mid)]
    assign = cluster_copies(copies, threshold=0.9)
    assert assign[100][0] == 0                 # b folds into a
    assert assign[200][0] == 200               # c stands alone
    assert assign[300][0] == 300               # last copy always kept


def test_cluster_copies_always_keeps_first_and_last():
    """Terminal copies of an array diverge from the interior, so they are never
    represented by another copy."""
    s = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
    copies = [(0, s), (10, s), (20, s), (30, s)]
    assign = cluster_copies(copies, threshold=0.9)
    assert assign[0][0] == 0
    assert assign[30][0] == 30
    assert assign[10][0] == 0                  # interior collapses


def test_collapse_repeat_blocks_keeps_qualifier_subtypes_apart():
    """FBN1's plain and calcium-binding EGF-like copies are different subtypes.
    Grouping on family alone would let one represent the other."""
    unit = "ACDEFGHIKLMNPQRSTVWYACDEFGHIKL"
    seq = unit * 4
    blocks = [
        _block(0, 29, index=1),
        _block(30, 59, index=2),
        _block(60, 89, index=3, qualifier="calcium-binding"),
        _block(90, 119, index=4, qualifier="calcium-binding"),
    ]
    dropped, rows = collapse_repeat_blocks(blocks, seq, threshold=0.9)
    groups = {(r["family"], r["qualifier"]) for r in rows}
    assert groups == {("EGF-like", ""), ("EGF-like", "calcium-binding")}
    # every group has exactly two members, so both are first-and-last and
    # nothing can be dropped even though all four sequences are identical
    assert dropped == set()


# ---------------------------------------------------------------------------
# Packing
# ---------------------------------------------------------------------------
def test_pack_blocks_fills_to_the_cap_without_splitting():
    """Three 100-residue domains with a cap of 250 give two runs, and a run's
    span runs from its first block's start to its last block's end."""
    blocks = [_block(0, 99, index=1), _block(100, 199, index=2),
              _block(200, 299, index=3)]
    runs = pack_blocks(blocks, 250)
    assert [len(r) for r in runs] == [2, 1]
    assert (runs[0][0]["start"], runs[0][-1]["end"]) == (0, 199)
    assert (runs[1][0]["start"], runs[1][-1]["end"]) == (200, 299)


def test_pack_blocks_isolates_an_oversized_block():
    """A single domain larger than the cap becomes its own run, the only case
    where a fold gets split, and the caller tiles it."""
    blocks = [_block(0, 49, index=1), _block(50, 400, index=2)]
    runs = pack_blocks(blocks, 100)
    assert [len(r) for r in runs] == [1, 1]
    assert runs[1][0]["length"] == 351


def test_pack_blocks_rejects_a_nonpositive_cap():
    with pytest.raises(ValueError, match="max_len must be greater than 0"):
        pack_blocks([_block(0, 10)], 0)


def test_respect_disulfides_shifts_a_boundary():
    """FBN1 carries 155 disulfide bonds; separating the two cysteines of a pair
    leaves neither fragment able to form it."""
    blocks = [_block(0, 99, index=1), _block(100, 199, index=2),
              _block(200, 299, index=3)]
    runs = [[blocks[0], blocks[1]], [blocks[2]]]
    # a bond spanning the 199/200 boundary
    runs2, unresolved = respect_disulfides(runs, [(150, 250)], max_len=400)
    assert unresolved == 0
    assert (runs2[0][-1]["end"], runs2[1][0]["start"]) != (199, 200)


def test_respect_disulfides_reports_when_it_cannot_fix():
    """With no slack in the cap the boundary cannot move, and that is reported
    rather than silently ignored."""
    blocks = [_block(0, 99, index=1), _block(100, 199, index=2)]
    runs = [[blocks[0]], [blocks[1]]]
    _, unresolved = respect_disulfides(runs, [(50, 150)], max_len=100)
    assert unresolved == 1


# ---------------------------------------------------------------------------
# Unannotated regions
# ---------------------------------------------------------------------------
def test_unannotated_runs_ignores_gaps_inside_a_fragment():
    """Measured against emitted fragments, not raw blocks: FBN1's proline-rich
    390-448 gap sits inside the 81-1445 fragment and must not be re-emitted."""
    spans = [dict(start=80, end=1444)]
    runs = unannotated_runs(spans, 2871, min_len=30)
    assert [(r["start"], r["end"]) for r in runs] == [(0, 79), (1445, 2870)]


def test_unannotated_runs_respects_the_minimum_length():
    spans = [dict(start=0, end=99), dict(start=110, end=299)]
    assert unannotated_runs(spans, 300, min_len=30) == []
    assert len(unannotated_runs(spans, 300, min_len=5)) == 1


def test_find_period_detects_a_tandem_repeat():
    """An offset-identity scan finds a real period where exact k-mer search
    fails on degenerate copies."""
    unit = "ACDEFGHIKLMNPQRSTVWY"
    period, score = find_period(unit * 8, min_period=10, max_period=60)
    assert period == len(unit)
    assert score > 0.95


def test_find_period_returns_nothing_on_aperiodic_sequence():
    """MUC16's 12066-residue N-terminal region has no recoverable periodicity;
    the best offset identity there is 15.8%, barely above background. Claiming
    a period would mis-slice it."""
    import random
    random.seed(0)
    seq = "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(1200))
    period, score = find_period(seq, min_period=20, max_period=200)
    assert period is None
    assert score < 0.5


def test_handle_ndr_exclude_drops_everything():
    runs = [dict(start=0, end=199, length=200, type="ndr")]
    assert handle_ndr(runs, "A" * 200, "exclude", 100, 20) == []


def test_handle_ndr_tile_collapses_a_periodic_stretch():
    """A periodic unannotated stretch yields a representative unit plus its
    terminal copy, instead of a blind tiling of every copy."""
    unit = "ACDEFGHIKLMNPQRSTVWY"
    seq = unit * 10
    runs = [dict(start=0, end=len(seq) - 1, length=len(seq), type="ndr")]
    out = handle_ndr(runs, seq, "tile", max_len=1000, overlap=10)
    assert all(f["type"] == "ndr_repeat" for f in out)
    assert out[0]["length"] == len(unit)
    assert len(out) < 10


def test_handle_ndr_tile_falls_back_to_windows():
    """With no period to find, the policy degrades to the existing overlapping
    window behaviour rather than emitting one giant fragment."""
    import random
    random.seed(1)
    seq = "".join(random.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(600))
    runs = [dict(start=0, end=599, length=600, type="ndr")]
    out = handle_ndr(runs, seq, "tile", max_len=200, overlap=50)
    assert len(out) > 1
    assert all(f["length"] <= 200 for f in out)
