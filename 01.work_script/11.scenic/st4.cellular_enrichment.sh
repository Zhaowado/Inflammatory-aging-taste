for sample in `ls ../01.rds2loom/*.loom`
do
    name=`basename ${sample} | cut -d '.' -f 1-2`
    echo -e "/dellfsqd2/ST_OCEAN/USER/hankai/software/miniconda/envs/st-pipe/bin/pyscenic aucell \\
	${sample} \\
	../03.regulon_prediction/${name}.reg.csv \\
	--output ${name}.scenic.loom \\
	--num_workers 20">${name}.sh
#	qsub -cwd -l vf=60g,num_proc=20 -P P22Z25400N0209 -binding linear:20 -q st.q ${name}.sh
	qsub -cwd -l vf=100g,num_proc=20 -P st_supermem -binding linear:20 -q st_supermem.q ${name}.sh
done
