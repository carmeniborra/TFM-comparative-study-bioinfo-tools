srun \
        --job-name normalidad_test \
        --cpus-per-task 1 \
        --mem 2G \
        --partition short_idx \
        --time 00:30:00 \
        bash -lc 'module load Python; python ../../DOC/normalidad_test.py .'
