#!/bin/bash
#SBATCH -J dl_microbialdb
#SBATCH -p middle_idx
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -t 16:00:00
#SBATCH -o dl_microbialdb.out
#SBATCH -e dl_microbialdb.err
#SBATCH --chdir /scratch/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/krakenuniq

set -euo pipefail

DEST="/data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/krakenuniq/microbialdb"
URL_KDB="https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/database.kdb"
URL_TAR="https://genome-idx.s3.amazonaws.com/kraken/uniq/krakendb-2023-08-08-MICROBIAL/kuniq_microbialdb_minus_kdb.20230808.tgz"

mkdir -p "$DEST"

echo "[INFO] Descargando .tar.gz y .kdb a: $DEST"

WGET="wget -c --no-verbose --tries=0 --read-timeout=30 --timeout=30 --waitretry=30 --retry-connrefused"

$WGET -P "$DEST" "$URL_KDB"
$WGET -P "$DEST" "$URL_TAR"

TARFILE="$(basename "$URL_TAR")"
if [[ "$TARFILE" == *.tar.gz ]]; then
  echo "[INFO] Extrayendo $TARFILE ..."
  tar -xzf "$TARFILE"
fi

echo "[OK] Descarga y extracción completadas en: $DEST"
