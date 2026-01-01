import argparse
from .fasta import clean_fasta, slice_fasta
from .fetch import fetch_fasta_by_query, fetch_fasta_by_accession

def main():
    parser = argparse.ArgumentParser(prog="cpt", description="FASTA utilities for Uniprot workflows")

    subparsers = parser.add_subparsers(dest="command", required=True)

    fetch_parser = subparsers.add_parser("fetch", help="Use 'query' or 'acc' subcommands to fetch fasta files from uniprot.")
    fetch_sub = fetch_parser.add_subparsers(dest="fetch_mode", required=True)
    fasta_parser = subparsers.add_parser("fasta", help="Use 'clean' or 'slice' subcommands to process fasta files and associated")
    fasta_sub = fasta_parser.add_subparsers(dest="fasta_mode", required=True)

    fetch_query = fetch_sub.add_parser("query")
    fetch_query.add_argument("--query")
    fetch_query.add_argument("--taxonomy_id")
    fetch_query.add_argument("--reviewed")

    fetch_acc = fetch_sub.add_parser("acc")
    fetch_acc.add_argument("--file", required=True)

    fasta_clean = fasta_sub.add_parser("clean")
    fasta_clean.add_argument("--file", required=True)

    fasta_slice = fasta_sub.add_parser("slice")
    fasta_slice.add_argument("--file", required=True)
    fasta_slice.add_argument("--max_slice", required=True)
    fasta_slice.add_argument("--window")
    
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

    elif args.command == "fasta":
        if args.fasta_mode == "clean":
            clean_fasta(args.file)

        elif args.fasta_mode == "slice":
            slice_fasta(
                args.file,
                max_length=args.max_slice,
                window=args.window,
            )

if __name__ == "__main__":
    main()

