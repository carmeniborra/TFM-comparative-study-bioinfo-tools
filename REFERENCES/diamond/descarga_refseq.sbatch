#!/bin/bash
#SBATCH --job-name=RefSeq_Diamond_DB
#SBATCH --nodes=1
#SBATCH --partition middle_idx
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=diamond_makedb_%j.out
#SBATCH --error=diamond_makedb_%j.err

# --- CONFIGURACION ---
set -euo pipefail
module load singularity

# Directorio base de NCBI FTP
NCBI_FTP="ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/"
# Nombre del archivo de base de datos final
OUTPUT_DB="refseq_prokaryotic_proteins.faa.gz"
# Directorio temporal de descarga (creado por wget -m)
DOWNLOAD_DIR="ftp.ncbi.nlm.nih.gov/refseq/release/complete/"
#db de diamond
DB_DIR=/data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/diamond/
#imagen
IMG_DIAMOND="/data/ucct/bi/pipelines/singularity-images/diamond-2.1.12--hdb4b4cc_1.img"


echo "Iniciando la descarga de la base de datos de proteínas de RefSeq..."

# 1. Descargar todos los fragmentos (chunks) FASTA de proteína
# Se utiliza la ruta -m para descargar de forma recursiva al directorio DOWNLOAD_DIR
wget -A "complete.wp_protein.*.protein.faa.gz" -m "$NCBI_FTP"
if [ $? -ne 0 ]; then
    echo "ERROR: Falló la descarga de los archivos de proteínas."
    exit 1
fi
echo "Descarga de fragmentos FASTA completada."

# 2. Concatenar los fragmentos en el archivo final
echo "Concatenando fragmentos en $OUTPUT_DB..."
cat "$DOWNLOAD_DIR"/complete.wp_protein.*.protein.faa.gz > "$OUTPUT_DB"

if [ $? -ne 0 ]; then
    echo "ERROR: Fallo la concatenación de archivos."
    exit 1
fi
echo "Concatenación completada. Tamaño de $OUTPUT_DB: $(du -sh $OUTPUT_DB | awk '{print $1}')"

# 3. Eliminar el directorio temporal de descarga para ahorrar espacio (¡CRÍTICO!)
echo "Limpiando archivos fragmentados..."
rm -r ftp.ncbi.nlm.nih.gov
echo "Limpieza completada."

# 4. Descargar los archivos de taxonomía necesarios para DIAMOND
echo "Descargando archivos de mapeo y taxonomía de NCBI..."
# Mapeo Accesión a Taxonomía
wget "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz"
# Archivos de la estructura Taxonómica (nodes y names)
wget "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip"
unzip taxdmp.zip nodes.dmp names.dmp 
rm taxdmp.zip # Limpiar el zip

echo "¡Todos los archivos necesarios para la base de datos de DIAMOND han sido descargados y preparados!"
echo "Archivos listos: $OUTPUT_DB, prot.accession2taxid.gz, nodes.dmp, names.dmp"

#creación de la base de datos de DIAMOND, que puede ir aquí:
echo "Creando la base de datos de DIAMOND..."

singularity exec --cleanenv \
  -B "$DB_DIR":"$DB_DIR" \
  "$IMG_DIAMOND" \
  bash -lc "cd '$DB_DIR' && \
    diamond makedb \
      --in $OUTPUT_DB \
      --db diamond_refseq_prok \
      --taxonmap prot.accession2taxid.gz \
      --taxonnodes nodes.dmp \
      --taxonnames names.dmp \
      --threads ${SLURM_CPUS_PER_TASK:-16}"

if [ $? -ne 0 ]; then
    echo "ERROR: Falló la creación de la base de datos de DIAMOND (v 2.1.12)."
    exit 1
fi

echo "Base de datos de DIAMOND completada. El archivo es: $DB_DIR/diamond_refseq_prok.dmnd"
