import csv
import os
import math


def iter_fasta(fasta_file: os.PathLike):
    """Yield (identifier, gene, sequence) for each record in a FASTA file.

    Handles UniProt headers (>sp|ACC|NAME ... GN=GENE -> id=ACC, gene=GENE) and
    already-simplified/sliced headers (>ACC_NAME_1-308 -> id and gene both the
    first whitespace token), concatenating wrapped sequence lines.
    """
    ident, gene, seq = None, None, []
    with open(os.path.abspath(fasta_file)) as fh:
        for line in fh:
            if line.startswith(">"):
                if ident is not None:
                    yield ident, gene, "".join(seq)
                ident = line.split("|")[1] if "|" in line else line[1:].split()[0]
                gene = next((t[3:] for t in line.split() if t.startswith("GN=")), ident)
                seq = []
            else:
                seq.append(line.strip())
        if ident is not None:
            yield ident, gene, "".join(seq)


def clean_fasta(fasta_file: os.PathLike):

    abs_path = os.path.abspath(fasta_file)
    root, ext = os.path.splitext(abs_path)

    new_file = []
    seq = ""
    with open(abs_path, "r") as file:
        for line in file:
            if line.startswith(">"):
                if seq:
                    new_file.append(seq)
                seq = ""
                subline = line[line.find("|")+1:]
                acc_id = subline[:subline.find("|")]
                entry_name = subline[subline.find("|")+1:subline.find(" ")]
                new_header = f">{acc_id}_{entry_name}"
                new_file.append(new_header)
            else:
                seq += line.rstrip()
        if seq:
            new_file.append(seq)
    
    output_path = f"{root}_cleaned{ext}"
    with open(output_path, 'w') as file:
        file.write("\n".join(new_file))
    
    print(f"cleaned fasta file written to {output_path}")
    return output_path

def slice_fasta(fasta_file: os.PathLike, max_length: int, window: int | None):

    abs_path = os.path.abspath(fasta_file)
    root, ext = os.path.splitext(abs_path)

    if not window:
        window = 0

    max_length = int(max_length)
    
    window = int(window)

    if max_length <= 0:
        raise ValueError("max_length must be greater than 0")
    if window < 0:
        raise ValueError("window must be greater than or equal to 0")
    if window >= max_length:
        raise ValueError("window must be smaller than max_length")
    
    new_file = []
    current_header = ""
    sliced_header = ""
    sliced_fasta = ""
    seen_seq = False
    with open(abs_path, "r") as file:
        for line in file:
            line = line.rstrip()
            if line.startswith(">"):
                seen_seq = False
                current_header = line
                continue
            if len(line) <= max_length:
                if seen_seq:
                    raise Exception("The submitted fasta file does not have single line sequences")
                seen_seq = True
                new_file.append(current_header)
                new_file.append(line)
                continue 
            
            if seen_seq:
               raise Exception("The submitted fasta file does not have single line sequences")
            seen_seq = True 

            min_chunks = math.ceil(len(line)/max_length)
            step = int(math.ceil((len(line)-window)/min_chunks))
            chunk_len = step + window
            while chunk_len > max_length:
                min_chunks += 1
                step = int(math.ceil((len(line)-window)/min_chunks))
                chunk_len = step + window
            for k in range(min_chunks):
                start = k * step
                end = start + chunk_len
                if k == min_chunks-1:
                    sliced_fasta = line[start:len(line)]
                    sliced_header = f"{current_header},{start+1}-{len(line)}"
                    new_file.append(sliced_header)
                    new_file.append(sliced_fasta)
                else:
                    sliced_fasta = line[start:end]
                    sliced_header = f"{current_header},{start+1}-{end}"
                    new_file.append(sliced_header)
                    new_file.append(sliced_fasta) 

    output_path = f"{root}_sliced{ext}"
    with open(output_path, 'w') as file:
        file.write("\n".join(new_file))
    
    print(f"sliced fasta file written to {output_path}")
    return output_path

def reformat_sliced_windows(sliced_fasta: os.PathLike, gene_map_fasta: os.PathLike):
    """Rewrite `slice_fasta` output headers into a form downstream tooling
    (this toolkit's `boltz.generate_boltz_inputs`, and job-name parsing that
    splits `accession_gene__ligand`) can consume directly.

    `slice_fasta` (fed a `clean_fasta`d input) names each window after the
    UniProt *entry mnemonic* baked into the cleaned header
    (`>ACC_NAME,START-END`, e.g. `>Q9H251_CAD23_HUMAN,1-1802`) -- not the real
    gene symbol (`CAD23` is the mnemonic for `CDH23`; mnemonics and gene
    symbols frequently diverge). This looks the true gene symbol up by
    accession from `gene_map_fasta` (an uncleaned UniProt FASTA that still
    has `GN=` tags -- typically the same FASTA `clean_fasta`'s input came
    from) and rewrites each header to `{ACC}.w{START}-{END} GN={GENE}`. The
    period-joined window tag (rather than underscore-joined) keeps
    `accession_gene__ligand`-style job-name splitting working unmodified,
    since a UniProt accession may itself contain a hyphen (isoform suffixes,
    e.g. `Q99102-10`) but never a period.

    Writes `{root}_windows{ext}` plus a manifest CSV
    (`{root}_windows_manifest.csv`: window_id, accession, gene, start, end,
    window_len) and returns (fasta_path, manifest_path).
    """
    gene_map = {ident: gene for ident, gene, _ in iter_fasta(gene_map_fasta)}

    abs_path = os.path.abspath(sliced_fasta)
    root, ext = os.path.splitext(abs_path)

    records = []
    header, seq = None, ""
    with open(abs_path) as file:
        for line in file:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, seq))
                header, seq = line[1:], ""
            else:
                seq += line
        if header is not None:
            records.append((header, seq))

    rows = []
    new_file = []
    for header, seq in records:
        # header shape: "ACC_NAME,START-END" (ACC may itself contain a "-"
        # for isoforms, e.g. "Q99102-10_MUC4_HUMAN,4419-6419")
        acc_name, coord = header.rsplit(",", 1)
        start, end = coord.split("-")
        acc = acc_name.split("_", 1)[0]
        gene = gene_map.get(acc, acc)
        window_id = f"{acc}.w{start}-{end}"
        new_file.append(f">{window_id} GN={gene}")
        new_file.append(seq)
        rows.append(dict(window_id=window_id, accession=acc, gene=gene,
                         start=start, end=end, window_len=len(seq)))

    out_fasta = f"{root}_windows{ext}"
    with open(out_fasta, "w") as file:
        file.write("\n".join(new_file))

    manifest_path = f"{root}_windows_manifest.csv"
    with open(manifest_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["window_id", "accession", "gene",
                                           "start", "end", "window_len"])
        w.writeheader()
        w.writerows(rows)

    print(f"reformatted windows written to {out_fasta} ({len(rows)} windows)")
    print(f"manifest written to {manifest_path}")
    return out_fasta, manifest_path


def generate_bait_prey(prey_fasta: os.PathLike, bait: str, double_count: bool = False):
    """
    Generates a bait-prey pairing file in 'prey;bait' format.
    Validates that all bait proteins are present in the provided FASTA file.
    """ 
    base_name, ext = os.path.splitext(os.path.basename(prey_fasta))
    root, _ = os.path.splitext(os.path.abspath(prey_fasta))

    if ',' in bait: 
        bait_prots = {b.strip() for b in bait.split(',')}
    else:
        bait_prots = {bait.strip()}

    # Extract all unique prey headers from the FASTA file
    prey_prots = set()
    with open(prey_fasta, 'r') as rfile:
        for line in rfile:
            if line.startswith('>'):
                prey_prots.add(line.strip().lstrip('>'))

    # Validate that all baits exist in the FASTA file
    missing_baits = bait_prots - prey_prots
    if missing_baits:
        raise ValueError(f"The following bait proteins were not found in the FASTA file: {', '.join(missing_baits)}")

    output_path = f'{root}_complex.txt'
    with open(output_path, 'w') as wfile:
        for prey in prey_prots:
            for b in bait_prots:
                if not double_count and b == prey:
                    continue
                wfile.write(f'{prey};{b}\n')

    print(f'bait;prey file written to {output_path}')
    return output_path 
