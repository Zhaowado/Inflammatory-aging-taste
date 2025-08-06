# if csv is empty,check whether the ${name}.sh.o is finished
for sample in `ls ../01.rds2loom/*.loom`
do
    name=`basename ${sample} | cut -d '.' -f 1-2`
    echo -e "/dellfsqd2/ST_OCEAN/USER/hankai/software/miniconda/envs/st-pipe/bin/pyscenic ctx \\
	../02.network_inference/${name}.adj.tsv \\
	../00.database/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather \\
	../00.database/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather \\
	--annotations_fname ../00.database/motifs-v9-nr.mgi-m0.001-o0.0.tbl \\
	--expression_mtx_fname ${sample} \\
	--mode "dask_multiprocessing" \\
	--output ${name}.reg.csv \\
	--num_workers 20 \\
	--mask_dropouts">${name}.sh
#	if [[ "$name" != "M12_CV" && "$name" != "M12_FL" ]]; then
#		qsub -cwd -l vf=60g,num_proc=20 -P P22Z25400N0209 -binding linear:20 -q st.q ${name}.sh
#	else
#		echo "Skipping ${name} as it is M12_CV or M12_FL"
#	fi
#qsub -cwd -l vf=60g,num_proc=20 -P P22Z25400N0209 -binding linear:20 -q st.q ${name}.sh
qsub -cwd -l vf=100g,num_proc=20 -P st_supermem -binding linear:20 -q st_supermem.q ${name}.sh
done

#skipped files are both 100m+
#qsub -cwd -l vf=130g,num_proc=20 -P st_supermem -binding linear:20 -q st_supermem.q M12_CV.sh
#qsub -cwd -l vf=130g,num_proc=20 -P st_supermem -binding linear:20 -q st_supermem.q M12_FL.sh
