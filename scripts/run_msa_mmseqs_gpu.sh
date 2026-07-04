#!/usr/bin/env bash
# =============================================================================
# GPU-accelerated MMseqs2 MSA generation for the surfaceome cofolding screen.
#
# Produces one ColabFold-style a3m per input sequence, named <accession>.a3m so
# it drops straight into the Boltz-2 YAMLs written by make_boltz_inputs.py
# (--msa-mode precomputed --msa-dir msa).
#
# It drives MMseqs2's GPU search path (`--gpu 1`, released in MMseqs2 15+) via
# ColabFold's search wrapper, which is the maintained, correct way to reproduce
# the AlphaFold/Boltz UniRef30 + environmental MSA. A raw `mmseqs` GPU pipeline
# is documented at the bottom for reference.
#
# REQUIREMENTS
#   * NVIDIA GPU + recent driver (CUDA); check with `nvidia-smi`.
#   * mmseqs2 >= 15 built with GPU support (`mmseqs version`; the GPU build
#     prints CUDA in `mmseqs --help`). Install: conda install -c conda-forge -c
#     bioconda "mmseqs2>=15" ; or use the localcolabfold bundle.
#   * colabfold_search on PATH (from localcolabfold / pip install "colabfold").
#   * ColabFold databases downloaded once (~2 TB) into $DB_DIR:
#       colabfold_search --help  # see setup_databases.sh from ColabFold
#     then GPU-pad them (one-time, big speedup):
#       mmseqs makepaddedseqdb $DB_DIR/uniref30_2302_db          $DB_DIR/uniref30_2302_db_gpu
#       mmseqs makepaddedseqdb $DB_DIR/colabfold_envdb_202108_db $DB_DIR/colabfold_envdb_202108_db_gpu
#
# USAGE
#   ./run_msa_mmseqs_gpu.sh [FASTA] [DB_DIR] [OUT_DIR]
#   FASTA    input sequences        (default: surfaceome_seqs.fasta)
#   DB_DIR   ColabFold database dir (default: $HOME/colabfold_db)
#   OUT_DIR  a3m output directory   (default: msa)
#   env vars: THREADS (default nproc), GPUS (default 1), DB_LOAD_MODE (default 2)
# =============================================================================
set -euo pipefail

FASTA="${1:-sequences.fasta}"
DB_DIR="${2:-$HOME/colabfold_db}"
OUT_DIR="${3:-msa}"
THREADS="${THREADS:-$(nproc 2>/dev/null || echo 8)}"
GPUS="${GPUS:-1}"
DB_LOAD_MODE="${DB_LOAD_MODE:-2}"   # 2 = keep DB in RAM across the run

WORK="${OUT_DIR}/_work"
QFASTA="${WORK}/queries.fasta"      # headers normalized to bare accessions

# ---- preflight checks -------------------------------------------------------
command -v mmseqs >/dev/null            || { echo "ERROR: mmseqs not on PATH"; exit 1; }
command -v colabfold_search >/dev/null  || { echo "ERROR: colabfold_search not on PATH"; exit 1; }
[[ -f "$FASTA" ]]                       || { echo "ERROR: FASTA not found: $FASTA"; exit 1; }
[[ -d "$DB_DIR" ]]                      || { echo "ERROR: DB_DIR not found: $DB_DIR"; exit 1; }
if ! command -v nvidia-smi >/dev/null || ! nvidia-smi >/dev/null 2>&1; then
    echo "ERROR: no usable NVIDIA GPU detected (nvidia-smi failed)."; exit 1
fi
echo "mmseqs: $(mmseqs version 2>/dev/null || echo '?')  | GPUs: ${GPUS}  | threads: ${THREADS}"

mkdir -p "$WORK" "$OUT_DIR"

# ---- normalize FASTA headers to bare accessions (>sp|ACC|NAME -> >ACC) -------
# This makes the query order == accession order, so output a3m map cleanly.
awk '/^>/{n=split($1,a,"|"); print ">"(n>=2?a[2]:substr($1,2)); next} {print}' \
    "$FASTA" > "$QFASTA"
mapfile -t ACCS < <(grep "^>" "$QFASTA" | sed 's/^>//')
echo "Queries: ${#ACCS[@]} sequences"

# ---- GPU MSA search ---------------------------------------------------------
# colabfold_search wraps: createdb -> (GPU) search vs UniRef30 -> expandaln ->
# align -> result2msa, then the same vs the environmental DB, merged to a3m.
echo "Running colabfold_search on GPU..."
colabfold_search \
    --gpu 1 \
    --db-load-mode "$DB_LOAD_MODE" \
    --threads "$THREADS" \
    --use-env 1 \
    --db2 colabfold_envdb_202108_db \
    "$QFASTA" "$DB_DIR" "$WORK/raw"

# ---- rename output a3m -> <accession>.a3m -----------------------------------
# colabfold_search emits either <accession>.a3m (header naming) or 0-based
# <index>.a3m in query order. Handle both deterministically.
echo "Collecting a3m into ${OUT_DIR}/ ..."
n=0
for i in "${!ACCS[@]}"; do
    acc="${ACCS[$i]}"
    if   [[ -f "$WORK/raw/${acc}.a3m" ]]; then src="$WORK/raw/${acc}.a3m"
    elif [[ -f "$WORK/raw/${i}.a3m"   ]]; then src="$WORK/raw/${i}.a3m"
    else echo "  WARN: no a3m for ${acc} (index ${i})"; continue; fi
    cp -f "$src" "$OUT_DIR/${acc}.a3m"
    n=$((n+1))
done
echo "Done: wrote ${n}/${#ACCS[@]} a3m files to ${OUT_DIR}/"
[[ "$n" -eq "${#ACCS[@]}" ]] || echo "WARNING: ${#ACCS[@]}-n=$(( ${#ACCS[@]} - n )) sequences missing an MSA."

# =============================================================================
# REFERENCE: raw MMseqs2 GPU pipeline (if you prefer not to use colabfold_search)
# -----------------------------------------------------------------------------
# UNIREF=$DB_DIR/uniref30_2302_db_gpu ; ENVDB=$DB_DIR/colabfold_envdb_202108_db_gpu
# mmseqs createdb "$QFASTA" "$WORK/qdb"
# mmseqs search "$WORK/qdb" "$UNIREF" "$WORK/res_u" "$WORK/tmp" \
#     --gpu 1 --num-iterations 3 -a --e-profile 1e-1 --max-seqs 10000 -s 8
# mmseqs expandaln "$WORK/qdb" "$UNIREF" "$WORK/res_u" "$UNIREF" "$WORK/res_u_exp" --gpu 1
# mmseqs result2msa "$WORK/qdb" "$UNIREF" "$WORK/res_u_exp" "$WORK/uniref.a3m" \
#     --msa-format-mode 6
# # ...repeat search/expandaln/result2msa vs $ENVDB -> env.a3m, then:
# mmseqs mergedbs "$WORK/qdb" "$WORK/final.a3m" "$WORK/uniref.a3m" "$WORK/env.a3m"
# mmseqs unpackdb "$WORK/final.a3m" "$OUT_DIR" --unpack-name-mode 1 --unpack-suffix .a3m
# # (--unpack-name-mode 1 names files by the query FASTA header == accession)
# =============================================================================
