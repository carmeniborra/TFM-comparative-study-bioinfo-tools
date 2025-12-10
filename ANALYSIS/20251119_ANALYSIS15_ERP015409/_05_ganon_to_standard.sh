srun \
        --job-name GANONTOSTD \
        --cpus-per-task 2 \
        --mem 2G \
        --partition short_idx \
        --time 01:00:00 \
        bash -lc 'module load Python; python ../../DOC/ganon_to_standard.py results/ganon/ganon_db11_combined_reports.txt  --taxdump-dir ../../REFERENCES/taxdump/2025-09-01_taxdump  --drop-root   -o results/taxpasta/rel_abundance/ganon_db11.tsv.genus.tsv.study.tsv'
