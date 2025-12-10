#!/bin/bash
set -euo pipefail
# 1. Obtener el software de mOTUs

# Ruta de la imagen
IMG="/data/ucct/bi/pipelines/singularity-images/motus_3.1.0.sif"

# Solo descargar si no existe
if [ ! -f "$IMG" ]; then
	echo "Descargando imagen de mOTUs..."
	cd /data/ucct/bi/pipelines/singularity-images
	wget "https://depot.galaxyproject.org/singularity/motus%3A3.1.0--pyhdfd78af_0" \
	-O motus_3.1.0.sif
	echo "Imagen guardada correctamente"
else
	echo "Imagen ya existe en $IMG"
fi

# 2. Descargar la base de datos
	# 2.1 crear carpetas necesarias
mkdir -p /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/tmp_motus
touch /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/tmp_motus/db_mOTU.tar.gz

	# 2.2 descargar el tar.gz dentro del contenedor
module load singularity

echo "Descargando la base de datos de mOTUS..."
echo "Ejecuta _01_extracciondb.sh para solucionar el fallo cuando salga"
singularity exec \
        --bind /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/tmp_motus/db_mOTU.tar.gz:/usr/local/lib/python3.9/site-packages/motus/db_mOTU.tar.gz \
        /data/ucct/bi/pipelines/singularity-images/motus_3.1.0.sif \
        motus downloadDB \
        --db /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/db_mOTU
