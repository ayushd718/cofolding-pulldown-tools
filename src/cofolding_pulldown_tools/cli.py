import argparse
from .fasta import clean_fasta, slice_fasta, generate_bait_prey, reformat_sliced_windows
from .fetch import fetch_fasta_by_query, fetch_fasta_by_accession, fetch_fasta_by_gene_names
from .analysis import generate_csv
from .domains import split_domains_by_afdb
from .boltz import generate_boltz_inputs
from .cluster import dedup_sequences

def main():
    parser = argparse.ArgumentParser(prog="cpt", description="FASTA utilities for Uniprot workflows")

    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch_parser = subparsers.add_parser("fetch", help="Use 'query' or 'acc' subcommands to fetch fasta files from uniprot.")
    fetch_sub = fetch_parser.add_subparsers(dest="fetch_mode", required=True)
    fasta_parser = subparsers.add_parser("fasta", help="Use 'clean', 'slice', 'complex', 'domains', 'reformat-windows' or 'dedup' subcommands to process fasta files")
    fasta_sub = fasta_parser.add_subparsers(dest="fasta_mode", required=True)

    boltz_parser = subparsers.add_parser("boltz", help="Generate Boltz-2 cofolding inputs")
    boltz_sub = boltz_parser.add_subparsers(dest="boltz_mode", required=True)

    analyze_parser = subparsers.add_parser("analyze", help="Analyze cofolding results and generate a summary CSV")
    analyze_parser.add_argument("--dir", required=True, help="Directory containing cofolding output folders")

    fetch_query = fetch_sub.add_parser("query")
    fetch_query.add_argument("--query")
    fetch_query.add_argument("--taxonomy_id", type=int)
    fetch_query.add_argument("--reviewed")

    fetch_acc = fetch_sub.add_parser("acc")
    fetch_acc.add_argument("--file", required=True)

    fetch_genes = fetch_sub.add_parser("genes")
    fetch_genes.add_argument("--file", required=True,
                             help="Text file with one gene symbol per line")
    fetch_genes.add_argument("--organism_id", type=int, default=9606)
    fetch_genes.add_argument("--reviewed_only", action="store_true",
                             help="Skip unreviewed (TrEMBL) fallback")
    fetch_genes.add_argument("--delay", type=float, default=0.15,
                             help="Politeness delay between requests (seconds)")

    fasta_clean = fasta_sub.add_parser("clean")
    fasta_clean.add_argument("--file", required=True)

    fasta_slice = fasta_sub.add_parser("slice")
    fasta_slice.add_argument("--file", required=True)
    fasta_slice.add_argument("--max_slice", type=int, required=True)
    fasta_slice.add_argument("--window", type=int)

    fasta_complex = fasta_sub.add_parser("complex")
    fasta_complex.add_argument("--file", required=True)
    fasta_complex.add_argument("--bait", required=True)
    fasta_complex.add_argument("--double_count", action='store_true')

    fasta_domains = fasta_sub.add_parser("domains",
        help="Split proteins into domain fragments using AlphaFold DB PAE")
    fasta_domains.add_argument("--file", required=True, help="Input FASTA")
    fasta_domains.add_argument("--out_dir", default="domains")
    fasta_domains.add_argument("--ndr", choices=["exclude", "keep", "pad"], default="exclude",
        help="Non-domain residues: exclude (drop), keep (own fragments), pad (absorb into domains)")
    fasta_domains.add_argument("--min_domain_len", type=int, default=40)
    fasta_domains.add_argument("--min_ndr_len", type=int, default=30)
    fasta_domains.add_argument("--pad_len", type=int, default=10)
    fasta_domains.add_argument("--max_domain_len", type=int, default=0,
        help="If >0, split domains longer than this into overlapping windows")
    fasta_domains.add_argument("--overlap", type=int, default=50)
    fasta_domains.add_argument("--resolution", type=float, default=0.3,
        help="Louvain resolution; 0.1-0.5 stable, higher over-splits domains")
    fasta_domains.add_argument("--limit", type=int, default=0, help="Process only first N (debug)")

    fasta_reformat_windows = fasta_sub.add_parser("reformat-windows",
        help="Rewrite `slice` output headers (ACC.wSTART-END GN=GENE) with real gene symbols")
    fasta_reformat_windows.add_argument("--file", required=True, help="FASTA produced by `cpt fasta slice`")
    fasta_reformat_windows.add_argument("--gene_map", required=True,
        help="Uncleaned UniProt FASTA (has GN= tags) to look up real gene symbols by accession")

    fasta_dedup = fasta_sub.add_parser("dedup",
        help="Cluster near-identical sequences (e.g. tandem-repeat windows) and keep one representative per cluster")
    fasta_dedup.add_argument("--file", required=True, help="Input FASTA")
    fasta_dedup.add_argument("--min_seq_id", type=float, default=0.95,
        help="Identity threshold to cluster two sequences together (default 0.95)")
    fasta_dedup.add_argument("--coverage", type=float, default=0.8,
        help="Minimum alignment coverage (mmseqs -c)")
    fasta_dedup.add_argument("--cov_mode", type=int, default=1,
        help="mmseqs --cov-mode (default 1: coverage of the shorter sequence)")

    boltz_inputs = boltz_sub.add_parser("inputs",
        help="Write one Boltz-2 YAML per protein paired with a ligand")
    boltz_inputs.add_argument("--file", required=True, help="Input FASTA")
    boltz_inputs.add_argument("--ligands", required=True, help=".smiles file (>name / SMILES)")
    boltz_inputs.add_argument("--ligand_name", default=None, help="Ligand header to use (default: first)")
    boltz_inputs.add_argument("--out_dir", default="boltz_inputs")
    boltz_inputs.add_argument("--msa_mode", choices=["precomputed", "empty", "server"],
        default="precomputed")
    boltz_inputs.add_argument("--msa_dir", default="msa")
    boltz_inputs.add_argument("--max_len", type=int, default=2500,
        help="Skip proteins longer than this (0 = no limit)")
    boltz_inputs.add_argument("--no_affinity", action="store_true")

    args = parser.parse_args()
    
    if args.command == "fetch":
        if args.fetch_mode == "query":
            fetch_fasta_by_query(
                query=args.query,
                taxonomy_id=args.taxonomy_id,
                reviewed=args.reviewed,
            )

        elif args.fetch_mode == "acc":
            fetch_fasta_by_accession(args.file)

        elif args.fetch_mode == "genes":
            fetch_fasta_by_gene_names(
                args.file,
                organism_id=args.organism_id,
                reviewed_only=args.reviewed_only,
                delay=args.delay,
            )

    elif args.command == "fasta":
        if args.fasta_mode == "clean":
            clean_fasta(args.file)

        elif args.fasta_mode == "slice":
            slice_fasta(
                args.file,
                max_length=args.max_slice,
                window=args.window
            )
        elif args.fasta_mode == "complex":
            generate_bait_prey(
                prey_fasta=args.file,
                bait=args.bait,
                double_count=args.double_count
            )
        elif args.fasta_mode == "domains":
            split_domains_by_afdb(
                args.file,
                out_dir=args.out_dir,
                ndr_mode=args.ndr,
                min_domain_len=args.min_domain_len,
                min_ndr_len=args.min_ndr_len,
                pad_len=args.pad_len,
                max_domain_len=args.max_domain_len,
                overlap=args.overlap,
                resolution=args.resolution,
                limit=args.limit,
            )
        elif args.fasta_mode == "reformat-windows":
            reformat_sliced_windows(
                args.file,
                gene_map_fasta=args.gene_map,
            )
        elif args.fasta_mode == "dedup":
            dedup_sequences(
                args.file,
                min_seq_id=args.min_seq_id,
                coverage=args.coverage,
                cov_mode=args.cov_mode,
            )

    elif args.command == "boltz":
        if args.boltz_mode == "inputs":
            generate_boltz_inputs(
                args.file,
                args.ligands,
                ligand_name=args.ligand_name,
                out_dir=args.out_dir,
                msa_mode=args.msa_mode,
                msa_dir=args.msa_dir,
                max_len=args.max_len,
                affinity=not args.no_affinity,
            )

    elif args.command == "analyze":
        generate_csv(args.dir)


if __name__ == "__main__":
    main()

