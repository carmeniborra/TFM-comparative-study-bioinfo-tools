srun \
        --job-name metrics \
        --cpus-per-task 1 \
        --mem 2G \
        --partition short_idx \
        --time 00:30:00 \
        bash -lc 'module load Python; python ../../DOC/metricas.py --gs aligned_tables/gold_standard.aligned.tsv --tools_dir aligned_tables --out metricas_por_muestra.tsv'
