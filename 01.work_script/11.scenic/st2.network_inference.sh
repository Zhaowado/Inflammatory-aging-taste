for sample in `ls ../01.rds2loom/*.loom`
do
    name=`basename ${sample} | cut -d '.' -f 1-2`
    echo -e "/dellfsqd2/ST_OCEAN/USER/hankai/software/miniconda/envs/st-pipe/bin/pyscenic grn \\
	--num_workers 20 \\
	--output ${name}.adj.tsv \\
	--method grnboost2 \\
	${sample} \\
	../00.database/allTFs_mm.txt">${name}.sh
qsub -cwd -l vf=60g,num_proc=20 -P P22Z25400N0209 -binding linear:20 ${name}.sh
done
