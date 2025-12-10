srun \
        --job-name MOTUSTOGENUS \
        --cpus-per-task 2 \
        --mem 2G \
        --partition short_idx \
        --time 01:00:00 \
        bash -lc 'module load Python; python ../../DOC/motus_to_genus.py results/taxpasta/motus_db8.tsv --out ./results/taxpasta/filtered/motus_db8.tsv.genus.tsv --taxdump-dir ../../REFERENCES/taxdump/2025-09-01_taxdump'
