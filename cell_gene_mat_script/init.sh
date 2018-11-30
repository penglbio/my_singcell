source activate Scell
conda env export > doc/environment.yml

if [ ! -d step7 ]; then
	mkdir -p step1 step2 step3 step4 step5 step6 step7 fastqc raw_data barcode 
	ls raw_data/*.fq|rush 'seqkit split -p 100 -O {/}/{%@.*(R.).fq} -f {} '
fi
