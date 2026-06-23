import requests 
import os

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
    

def fasta_url_checker(url, raw_fasta):
    if len(raw_fasta.strip()) > 0:
        first_char = raw_fasta.strip()[0]
    else:
        return False

    if first_char == '>':
        return True
    else:
        return False