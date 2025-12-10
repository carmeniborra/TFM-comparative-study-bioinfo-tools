#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --mem=300G
#SBATCH --time=30:00:00
#SBATCH --output=kaiju_db_%j.log
#SBATCH --partition=middle_idx
#SBATCH --chdir /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/kaiju

set -euo pipefail
module load singularity

DB=/data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/kaiju
IMG=/data/ucct/bi/pipelines/singularity-images/kaiju-1.10.0--h43eeafb_0.img
mkdir -p "$DB"

# Kaiju usa el CWD; trabajamos dentro de $DB y montamos la MISMA ruta
singularity exec --cleanenv \
  -B "$DB":"$DB" \
  "$IMG" \
  bash -lc "cd '$DB' && kaiju-makedb -s refseq -t ${SLURM_CPUS_PER_TASK:-16} -v"