import pytest
from cofolding_pulldown_tools.fasta import iter_fasta
from cofolding_pulldown_tools.domains import (
    contiguous_segments,
    build_fragments,
    windowize,
)

# labels: 0-49 domain, 50-59 short linker, 60-99 domain
SEGS = [(0, 49, 0), (50, 59, 1), (60, 99, 2)]
N = 100


def test_iter_fasta_uniprot_and_sliced(tmp_path):
    fa = tmp_path / "x.fasta"
    fa.write_text(
        ">sp|P12345|AAA_HUMAN Protein OS=Homo sapiens OX=9606 GN=AAA PE=1 SV=1\n"
        "MKTA\nYIAK\n"                       # wrapped sequence -> should concatenate
        ">P12345_AAA_1-50\n"                 # already-sliced header (no pipes / GN)
        "MSEQ\n"
    )
    recs = list(iter_fasta(str(fa)))
    assert recs[0] == ("P12345", "AAA", "MKTAYIAK")
    assert recs[1] == ("P12345_AAA_1-50", "P12345_AAA_1-50", "MSEQ")


def test_contiguous_segments():
    labels = [0, 0, 0, 1, 1, 2]
    assert contiguous_segments(labels) == [(0, 2, 0), (3, 4, 1), (5, 5, 2)]


def test_build_fragments_exclude():
    frags = build_fragments(SEGS, N, min_domain_len=40, ndr_mode="exclude")
    assert [(f["start"], f["end"], f["type"]) for f in frags] == [
        (0, 49, "domain"), (60, 99, "domain")]


def test_build_fragments_keep_respects_min_ndr_len():
    # linker is 10 long: kept when min_ndr_len<=10, dropped when >10
    kept = build_fragments(SEGS, N, min_domain_len=40, ndr_mode="keep", min_ndr_len=5)
    assert (50, 59, "ndr") in [(f["start"], f["end"], f["type"]) for f in kept]

    dropped = build_fragments(SEGS, N, min_domain_len=40, ndr_mode="keep", min_ndr_len=30)
    assert all(f["type"] == "domain" for f in dropped)


def test_build_fragments_keep_coalesces_microclusters():
    # two short adjacent clusters between domains should merge into one NDR run
    segs = [(0, 49, 0), (50, 59, 1), (60, 74, 2), (75, 114, 3)]
    frags = build_fragments(segs, 115, min_domain_len=40, ndr_mode="keep", min_ndr_len=20)
    ndrs = [(f["start"], f["end"]) for f in frags if f["type"] == "ndr"]
    assert ndrs == [(50, 74)]                # 50-59 and 60-74 coalesced


def test_build_fragments_pad_no_overlap():
    frags = build_fragments(SEGS, N, min_domain_len=40, ndr_mode="pad", pad_len=5)
    spans = [(f["start"], f["end"]) for f in frags]
    assert spans == [(0, 54), (55, 99)]      # linker split between neighbors, gap-free


def test_windowize():
    frag = dict(start=0, end=99, length=100, type="domain")
    pieces = windowize(frag, max_len=60, overlap=10)
    assert [(p["start"], p["end"]) for p in pieces] == [(0, 59), (50, 99)]
    # short fragment is returned unchanged
    assert windowize(frag, max_len=0, overlap=10) == [frag]


def test_cluster_pae_two_blocks():
    np = pytest.importorskip("numpy")
    pytest.importorskip("networkx")
    from cofolding_pulldown_tools.domains import cluster_pae

    n = 80
    pae = np.full((n, n), 30.0)      # high PAE everywhere (unrelated)
    pae[:40, :40] = 1.0             # block 1: confident
    pae[40:, 40:] = 1.0             # block 2: confident
    np.fill_diagonal(pae, 0.0)

    labels = cluster_pae(pae, resolution=0.3)
    assert len(set(labels)) == 2
    domains = [s for s in contiguous_segments(labels) if s[1] - s[0] + 1 >= 40]
    assert len(domains) == 2
