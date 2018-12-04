## run on HPC ##
nohup snakemake --cluster "qsub -q fat" -j 100 -rp --latency-wait 3600 >> nohup.log 2>&1 &
