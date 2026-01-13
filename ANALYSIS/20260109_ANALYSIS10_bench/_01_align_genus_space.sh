srun \
	--job-name aligntables \
	--cpus-per-task 1 \
  	--mem 2G \
  	--partition short_idx \
  	--time 01:00:00 \
  	bash -lc 'module load Python; \
		  python ../../DOC/align_genus_space.py --inputs taxpasta/rel_abundance/*study.tsv gold_standard.tsv --outdir aligned_tables && \
		  rm -rf aligned_tables/kraken2_db2*'
