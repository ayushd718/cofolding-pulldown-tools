import os
import math

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
