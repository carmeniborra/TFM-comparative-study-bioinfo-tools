#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 16
#SBATCH --mem 200G
#SBATCH --time 24:00:00
#SBATCH --partition middle_idx
#SBATCH --output 20250924_ganon_build_species1.log
#SBATCH --chdir /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/ganon

set -euo pipefail

module purge
module load singularity

# Carpeta donde guardarás las bases
export DB=/data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/ganon/db
mkdir -p "$DB"/{cache,tmp}

# Variables
export XDG_CACHE_HOME="$DB/cache"
export TMPDIR="/data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/ganon/tmp"
mkdir -p "$TMPDIR"
export IMG=/data/ucct/bi/pipelines/singularity-images/singularity-ganon-2.0.0--py39ha35b9be_0.img

# Asegura variables dentro del contenedor pese a --cleanenv
export SINGULARITYENV_XDG_CACHE_HOME="$XDG_CACHE_HOME"
export SINGULARITYENV_TMPDIR="$TMPDIR"

# Ejecutar ganon build en el contenedor
singularity exec --cleanenv   -B "$DB":"$DB"   -B "$TMPDIR":"$TMPDIR"   -H "$PWD"   "$IMG"   ganon build --source refseq               --organism-group archaea bacteria fungi viral               --complete-genomes               --genome-updater "-A 'species:1'"               --threads "${SLURM_CPUS_PER_TASK:-16}"               --db-prefix "$DB/abfv_rs_t1s"               --tmp-dir "$TMPDIR"
