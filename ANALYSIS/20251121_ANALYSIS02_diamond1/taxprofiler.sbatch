#!/bin/bash
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 4
#SBATCH --mem 6G
#SBATCH --time 120:00:00
#SBATCH --partition long_idx
#SBATCH --output 20251204_taxprofiler.log
#SBATCH --chdir /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/ANALYSIS/20251121_ANALYSIS16_ERP015409_diamond

# module load Nextflow singularity
module purge
module load Nextflow/25.04.6
module load singularity/3.7.1

mkdir -p ./results

export NXF_OPTS="-Xms1G -Xmx6G"

nextflow run /data/ucct/bi/pipelines/nf-core-taxprofiler/nf-core-taxprofiler-1.2.4 \
    -profile singularity \
    -c ../../DOC/taxprofiler.config \
    --input samplesheet.csv \
    --outdir ./results \
    --databases ../../DOC/databasesheet.csv \
    --multiqc_title ERP015409 \
    --preprocessing_qc_tool fastqc \
    --perform_shortread_qc true --shortread_qc_tool fastp --shortread_qc_mergepairs false --shortread_qc_includeunmerged false --shortread_qc_minlength 50 \
    --perform_shortread_complexityfilter true --shortread_complexityfilter_tool bbduk --shortread_complexityfilter_entropy 0.5 --shortread_complexityfilter_bbduk_windowsize 50 --save_complexityfiltered_reads true \
    --perform_shortread_hostremoval true --hostremoval_reference /data/ucct/bi/references/eukaria/homo_sapiens/hg38/NCBI/genome/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz \
    --shortread_hostremoval_index ../../REFERENCES/bowtie2/ \
    --run_kraken2 false \
    --run_bracken false \
    --run_centrifuge false \
    --run_kaiju false --kaiju_taxon_rank genus \
    --run_metaphlan false \
    --run_ganon false --ganon_report_rank genus \
    --run_motus false \
    --run_diamond true \
    --run_krakenuniq false \
    --run_malt false \
    --run_kmcp false \
    --run_profile_standardisation \
    --standardisation_taxpasta_format tsv \
    --taxpasta_taxonomy_dir ../../REFERENCES/taxdump/2025-09-01_taxdump/\
    --taxpasta_add_name --taxpasta_add_rank \
    --run_krona true --krona_taxonomy_directory /data/ucct/bi/research/20250625_TFM_CIBORRA_IC-SM_T/REFERENCES/krona/taxonomy_db --run_krona_combined \
    -resume
