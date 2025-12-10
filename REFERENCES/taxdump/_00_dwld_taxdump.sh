#!/bin/bash

set -euo pipefail
#con e si cualquier comando devuelve error, el script se detiene, con u, usar una variable no definida provoca error,
# y con o, en una tuberia, el exitcode sera el del primer comando que falle
 
echo " Iniciando descarga de taxdump ..."

# Parámetros
DATE=2025-09-01
BASE="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump_archive"
SNAP="new_taxdump_${DATE}.zip"

# Carpeta destino
DEST="./${DATE}_taxdump"
mkdir -p "$DEST"
cd "$DEST"

# Requisitos
command -v unzip >/dev/null || { echo "Falta 'unzip' (sudo apt-get install unzip)"; exit 1; }
command -v wget  >/dev/null || { echo "Falta 'wget' (sudo apt-get install wget)"; exit 1; }

# 1) Descargar el ZIP del snapshot
echo " Descargando ${SNAP} ..."
wget -q "${BASE}/${SNAP}" -O "${SNAP}"

# 2) Extraer SOLO lo necesario para Taxpasta (son .dmp, no .gz)
echo " Extrayendo names.dmp, nodes.dmp, merged.dmp, delnodes.dmp ..."
unzip -jo "${SNAP}" names.dmp nodes.dmp merged.dmp delnodes.dmp

echo "Archivos extraídos:"
ls -lh names.dmp nodes.dmp merged.dmp delnodes.dmp || true

echo " ✅ Listo. Usa esta carpeta con: --taxpasta_taxonomy_dir $(pwd)"
