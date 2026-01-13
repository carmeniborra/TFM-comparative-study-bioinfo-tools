srun \
	--job-name supdatatogs \
	--cpus-per-task 1 \
  	--mem 2G \
 	--partition short_idx \
  	--time 01:00:00 \
  	bash -lc 'module load Python; python ../../DOC/supdata_to_gs.py ../../RAW/ERP015409/sup_data_3_ERP015409.xlsx gold_standard.tsv'
