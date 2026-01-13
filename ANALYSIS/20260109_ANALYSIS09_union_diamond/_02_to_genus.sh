srun \
	--job-name TOGENUS \
	--cpus-per-task 2 \
  	--mem 2G \
  	--partition short_idx \
  	--time 01:00:00 \
  	bash -lc 'module load Python; python ../../DOC/to_genus.py results/taxpasta/ --pattern "*.tsv" --outdir results/taxpasta/filtered --taxdump-dir ../../REFERENCES/taxdump/2025-09-01_taxdump'
