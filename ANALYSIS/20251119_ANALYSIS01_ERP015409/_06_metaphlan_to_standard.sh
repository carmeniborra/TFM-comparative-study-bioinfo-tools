srun \
        --job-name METAPHLANTOSTD \
        --cpus-per-task 2 \
        --mem 2G \
        --partition short_idx \
        --time 01:00:00 \
        bash -lc 'module load Python; python ../../DOC/metaphlan_to_standard.py results/metaphlan/metaphlan_db4_combined_reports.txt --taxdump-dir ../../REFERENCES/taxdump/2025-09-01_taxdump -o results/taxpasta/rel_abundance/metaphlan_db4.tsv.genus.tsv.study.tsv'
