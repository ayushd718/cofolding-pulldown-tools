import os
import csv
import json
import subprocess
from pathlib import Path

def process_directory(subdir_path: Path):
    # Find PDB/CIF and PAE files
    pdb_file = None
    pae_file = None
    for f in subdir_path.iterdir():
        if f.suffix in [".pdb", ".cif"]:
            pdb_file = f
        elif f.suffix in [".json", ".npz"]:
            # Need a better heuristic? 
            # Looking at ipsae.py usage, JSON is for PDB/CIF, NPZ for Boltz
            pae_file = f
    
    if not pdb_file or not pae_file:
        print(f"Skipping {subdir_path}: Missing PDB/CIF or PAE file.")
        return

    # Call ipsae.py
    # Use default cutoffs 10, 15
    subprocess.run(["python", "public/ipsae.py", str(pae_file), str(pdb_file), "10", "15"], check=True)

def generate_csv(output_dir: os.PathLike):
    output_dir = Path(output_dir)
    data = []
    
    for subdir in output_dir.iterdir():
        if not subdir.is_dir():
            continue
            
        process_directory(subdir)
        
        # ---- Parse ranking_debug.json ----
        ranking_file = subdir / "ranking_debug.json"
        max_iptm_ptm = ""
        max_iptm = ""
        min_iptm = ""
        if ranking_file.exists():
            try:
                with open(ranking_file, "r") as f:
                    ranking_data = json.load(f)
                    if "iptm+ptm" in ranking_data and isinstance(ranking_data["iptm+ptm"], dict):
                        vals = list(ranking_data["iptm+ptm"].values())
                        if vals: max_iptm_ptm = max(vals)
                    if "iptm" in ranking_data and isinstance(ranking_data["iptm"], dict):
                        iptm_vals = list(ranking_data["iptm"].values())
                        if iptm_vals:
                            max_iptm = max(iptm_vals)
                            min_iptm = min(iptm_vals)
            except json.JSONDecodeError:
                pass

        # ---- parse the generated .txt file ----
        ipsae_vals = []
        pdockq_vals = []
        
        # Look for _10_15.txt (based on cutoff 10 15)
        for f in subdir.glob("*_10_15.txt"):
            with open(f, "r") as file:
                lines = [ln.strip() for ln in file if ln.strip()]
            
            if len(lines) < 2:
                continue
                
            header = lines[0].split()
            try:
                type_idx = header.index("Type")
                ipsae_idx = header.index("ipSAE")
                pdockq_idx = header.index("pDockQ")
            except ValueError:
                continue

            for line in lines[1:]:
                cols = line.split()
                if len(cols) <= max(type_idx, ipsae_idx, pdockq_idx):
                    continue
                if cols[type_idx] != "max":
                    continue
                
                try:
                    ipsae_vals.append(float(cols[ipsae_idx]))
                    pdockq_vals.append(float(cols[pdockq_idx]))
                except ValueError:
                    continue

        if not ipsae_vals:
            continue

        protein1, protein2 = subdir.name.split("_and_", 1) if "_and_" in subdir.name else (subdir.name, "Unknown")
        
        data.append([
            protein1,
            protein2,
            max_iptm_ptm,
            max_iptm,
            min_iptm,
            max(ipsae_vals),
            min(ipsae_vals),
            max(pdockq_vals),
            min(pdockq_vals),
        ])

    # Write CSV
    with open(output_dir / "results_summary.csv", "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Protein 1", "Protein 2", 
            "Max iPTM+PTM", "Max iPTM", "Min iPTM", 
            "Max ipSAE", "Min ipSAE", "Max pDockQ", "Min pDockQ"
        ])
        writer.writerows(data)