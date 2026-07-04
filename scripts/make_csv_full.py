import os
import json
import csv

parent_directory = "./structures_output/"
output_csv = "STK11_scores_summary.csv"

data = []

for subdir in os.listdir(parent_directory):
    subdir_path = os.path.join(parent_directory, subdir)
    if not os.path.isdir(subdir_path):
        continue

    ranking_file = os.path.join(subdir_path, "ranking_debug.json")
    max_iptm_ptm = ""
    max_iptm = ""
    min_iptm = ""

    # ---- parse ranking_debug.json for iptm+ptm and iptm ----
    if os.path.exists(ranking_file):
        with open(ranking_file, "r") as f:
            try:
                ranking_data = json.load(f)

                # max iptm+ptm
                if "iptm+ptm" in ranking_data and isinstance(ranking_data["iptm+ptm"], dict):
                    vals = list(ranking_data["iptm+ptm"].values())
                    if vals:
                        max_iptm_ptm = max(vals)

                # max and min iptm
                if "iptm" in ranking_data and isinstance(ranking_data["iptm"], dict):
                    iptm_vals = list(ranking_data["iptm"].values())
                    if iptm_vals:
                        max_iptm = max(iptm_vals)
                        min_iptm = min(iptm_vals)

            except json.JSONDecodeError:
                print(f"Error reading JSON in {subdir}")

    # ---- parse ipSAE output text files for ipSAE and pDockQ ----
    ipsae_vals = []
    pdockq_vals = []

    # assume ipSAE outputs like: unrelaxed_model_1_multimer_v3_pred_0_15_15.txt
    for fname in os.listdir(subdir_path):
        if not fname.endswith("_15_15.txt"):
            continue
        ipsae_file = os.path.join(subdir_path, fname)

        with open(ipsae_file, "r") as f:
            lines = [ln.strip() for ln in f if ln.strip()]
        if len(lines) < 2:
            continue

        header = lines[0].split()
        # robust index lookup
        try:
            type_idx = header.index("Type")
            ipsae_idx = header.index("ipSAE")
            pdockq_idx = header.index("pDockQ")
        except ValueError:
            # header not in expected form
            print(f"Unexpected header in {ipsae_file}: {header}")
            continue

        # use only rows where Type == "max"
        for line in lines[1:]:
            cols = line.split()
            if len(cols) <= max(type_idx, ipsae_idx, pdockq_idx):
                continue
            if cols[type_idx] != "max":
                continue

            try:
                ipsae_val = float(cols[ipsae_idx])
                pdockq_val = float(cols[pdockq_idx])
            except ValueError:
                continue

            ipsae_vals.append(ipsae_val)
            pdockq_vals.append(pdockq_val)

    max_ipsae = max(ipsae_vals) if ipsae_vals else ""
    min_ipsae = min(ipsae_vals) if ipsae_vals else ""
    max_pdockq = max(pdockq_vals) if pdockq_vals else ""
    min_pdockq = min(pdockq_vals) if pdockq_vals else ""

    # ---- protein names from directory name ----
    if "_and_" in subdir:
        protein1, protein2 = subdir.split("_and_", 1)
    else:
        protein1, protein2 = subdir, "Unknown"

    # store one row per directory
    data.append([
        protein1,
        protein2,
        max_iptm_ptm,
        max_iptm,
        min_iptm,
        max_ipsae,
        min_ipsae,
        max_pdockq,
        min_pdockq,
    ])

# write CSV
with open(output_csv, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow([
        "Protein 1",
        "Protein 2",
        "Max iPTM+PTM",
        "Max iPTM",
        "Min iPTM",
        "Max ipSAE",
        "Min ipSAE",
        "Max pDockQ",
        "Min pDockQ",
    ])
    writer.writerows(data)

print(f"CSV file saved: {output_csv}")

