srun \
        --job-name test_no_parametricos \
        --cpus-per-task 1 \
        --mem 2G \
        --partition short_idx \
        --time 00:30:00 \
        bash -lc 'module load Python; python ../../DOC/test_no_parametricos.py metricas_por_muestra.tsv --outdir tests_metricas'
