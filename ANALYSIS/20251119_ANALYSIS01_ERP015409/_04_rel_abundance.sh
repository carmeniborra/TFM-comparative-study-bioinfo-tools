srun \
        --job-name RELABUNDANCE \
        --cpus-per-task 2 \
        --mem 2G \
        --partition short_idx \
        --time 01:00:00 \
        bash -lc 'module load Python; python ../../DOC/rel_abundance.py results/taxpasta/filtered/ --pattern "*.tsv" --outdir results/taxpasta/rel_abundance --level genus --study-wide --sample-trim "_ERR"'
