import requests
import os
import csv
import time

def fetch_fasta_by_query(query: str | None, taxonomy_id: int | None, reviewed: str | None):
    
    terms = []
    
    filename = []

    if query:
        terms.append("+".join(query.split()))
        filename.append(query.strip())

    if taxonomy_id is not None:
        terms.append(f"organism_id:{taxonomy_id}")
        filename.append(f"{taxonomy_id}")

    if reviewed:
        user_rev = reviewed.lower()
        if user_rev not in {"true", "false"}:
            raise ValueError("reviewed must be 'true' or 'false'")
        terms.append(f"reviewed:{user_rev}")
        if user_rev == "true":
            filename.append("reviewed")
        if user_rev == "false":
            filename.append("unreviewed")

    if not terms:
        raise ValueError("At least one of query, taxonomy_id, or reviewed must be provided")

    query_string = "+".join(terms)
    url = f"https://rest.uniprot.org/uniprotkb/stream?query={query_string}&format=fasta"

    response = requests.get(url=url)

    if response.status_code != 200:
        raise Exception(f"Could not retrieve fasta from {url}, received error {response.status_code}")

    raw_fastas = response.text
    
    success = fasta_url_checker(url, raw_fastas)

    filename_string = "_".join(filename)

    if success:
        with open(f"{filename_string}.fasta", "w") as file:
            file.write(raw_fastas)
    else:
        raise Exception(f"Error fetching fasta at {url}")
    
    return print(f"fasta file fetched and written")
    
def fetch_fasta_by_accession(accession_file: str| os.PathLike):

    abs_path = os.path.abspath(accession_file)
    root, ext = os.path.splitext(abs_path)

    acc_list = []

    acc_fasta_list = []

    failed = []

    with open(abs_path, "r") as file:
        for line in file:
            line = line.strip()
            if line == "":
                continue
            acc_list.append(line)
    
    for i in acc_list:
        url = f"https://rest.uniprot.org/uniprotkb/{i}.fasta"
        response = requests.get(url=url, timeout=60)
        raw_fasta = response.text
        success = fasta_url_checker(url, raw_fasta)
        if success:
            acc_fasta_list.append(raw_fasta)
            print(f"fasta successfuly fetched from {url}")
        else:
            print(f"could not fetch fasta from {url}, written to {root}_failed{ext}")
            failed.append(f"{i} {url} {response.status_code}")
        #sleep_time = 1
        #time.sleep(sleep_time)

    if failed:
        with open(f"{root}_failed{ext}", "w") as file:
            file.write("\n".join(failed))
    
    if not acc_fasta_list:
        raise Exception("No queried accessions were valid or returned fastas")

    with open(f"{root}.fasta", 'w') as file:
        file.write("".join(acc_fasta_list))
    

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def all_fasta_records(text):
    """Split a multi-record FASTA string into a list of individual records
    (each 'header\\n...sequence...\\n'). Empty list if none."""
    text = text.strip()
    if not text.startswith(">"):
        return []
    return [(">" + chunk).rstrip() + "\n" for chunk in ("\n" + text).split("\n>")[1:]]


def first_fasta_record(text):
    """Return the first FASTA record (header + sequence) from a possibly
    multi-record string, or None if the text does not start with a record."""
    records = all_fasta_records(text)
    return records[0] if records else None


def _record_primary_gene(record):
    """Return the primary gene name (GN= field) from a UniProt FASTA record,
    or '' if absent."""
    header = record.splitlines()[0]
    for token in header.split():
        if token.startswith("GN="):
            return token[3:]
    return ""


def _select_record_for_gene(records, gene):
    """Pick the record whose UniProt primary gene name matches `gene`
    (case-insensitive); otherwise fall back to the first (highest-ranked)
    record. This avoids mis-resolving a symbol to an entry where it is only a
    synonym (e.g. gene_exact:CCR4 matching the NOCT entry)."""
    for record in records:
        if _record_primary_gene(record).upper() == gene.upper():
            return record
    return records[0] if records else None


def _fetch_gene_record(gene, organism_id, reviewed, retries=3):
    """Query UniProtKB for a gene symbol and return the FASTA record whose
    primary gene name matches (preferred), else the top hit; None if no hit.
    `reviewed` is 'true'/'false'."""
    query = f"gene_exact:{gene} AND organism_id:{organism_id} AND reviewed:{reviewed}"
    params = {"query": query, "format": "fasta", "size": 5}
    for attempt in range(retries):
        try:
            response = requests.get(UNIPROT_SEARCH, params=params, timeout=60)
        except requests.RequestException:
            time.sleep(2 ** attempt)
            continue
        if response.status_code == 200:
            record = _select_record_for_gene(all_fasta_records(response.text), gene)
            if record and fasta_url_checker(UNIPROT_SEARCH, record):
                return record
            return None
        if response.status_code in (429, 500, 502, 503):
            time.sleep(2 ** attempt)
            continue
        return None
    return None


def fetch_fasta_by_gene_names(gene_file: str | os.PathLike, organism_id: int = 9606,
                              reviewed_only: bool = False, delay: float = 0.15):
    """Fetch one canonical sequence per gene symbol listed in `gene_file`.

    For each symbol the reviewed (Swiss-Prot) entry is preferred; unless
    `reviewed_only`, symbols with no reviewed entry fall back to the best
    unreviewed (TrEMBL) entry and are flagged. Writes, next to the input:
      <root>.fasta         combined sequences (raw UniProt headers preserved)
      <root>_report.tsv    gene -> accession, status, header
      <root>_failed<ext>   gene symbols with no hit (re-runnable through this fn)
    """
    abs_path = os.path.abspath(gene_file)
    root, ext = os.path.splitext(abs_path)

    with open(abs_path, "r") as file:
        genes = [line.strip() for line in file if line.strip()]

    records = []
    report_rows = []
    failed = []
    counts = {"reviewed": 0, "unreviewed": 0, "failed": 0}

    for i, gene in enumerate(genes, 1):
        record = _fetch_gene_record(gene, organism_id, "true")
        status = "reviewed"
        if record is None and not reviewed_only:
            record = _fetch_gene_record(gene, organism_id, "false")
            status = "unreviewed"
        if record is None:
            status = "failed"
            failed.append(gene)
            header, accession = "", ""
        else:
            header = record.splitlines()[0]
            accession = header.split("|")[1] if "|" in header else ""
            records.append(record)
        counts[status] += 1
        report_rows.append({"gene": gene, "accession": accession,
                            "status": status, "header": header})
        if i % 25 == 0 or i == len(genes):
            print(f"  {i}/{len(genes)} reviewed={counts['reviewed']} "
                  f"unreviewed={counts['unreviewed']} failed={counts['failed']}",
                  flush=True)
        time.sleep(delay)

    if not records:
        raise Exception("No gene symbols resolved to a UniProt sequence")

    with open(f"{root}.fasta", "w") as file:
        file.write("".join(records))

    with open(f"{root}_report.tsv", "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["gene", "accession", "status", "header"],
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows)

    if failed:
        with open(f"{root}_failed{ext}", "w") as file:
            file.write("\n".join(failed) + "\n")

    print(f"fasta written to {root}.fasta "
          f"(reviewed={counts['reviewed']} unreviewed={counts['unreviewed']} "
          f"failed={counts['failed']})")


def fasta_url_checker(url, raw_fasta):
    if len(raw_fasta.strip()) > 0:
        first_char = raw_fasta.strip()[0]
    else:
        return False

    if first_char == '>':
        return True
    else:
        return False