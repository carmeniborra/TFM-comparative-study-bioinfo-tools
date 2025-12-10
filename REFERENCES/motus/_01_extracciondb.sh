#!/bin/bash
set -euo pipefail

# 1. Verificar tamano del archivo descargado (~3 GB)
echo "Tamano de la base de datos:"
ls -lh /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/tmp_motus/db_mOTU.tar.gz

# 2. Extraer manualmente (fuera del contenedor)
echo "Extrayendo la base de datos..."
tar -xzf /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/tmp_motus/db_mOTU.tar.gz \
    -C /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/motus

echo "TODO LISTO :) !! "
