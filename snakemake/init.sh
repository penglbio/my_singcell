source activate Scell
conda env export > doc/environment.yml

if [ ! -d step7 ]; then
	mkdir -p step1 step2 step3 step4 step5 step6 step7 fastqc raw_data barcode figure stat 
fi

ln -s $1 R1.fastq
ln -s $2 R2.fastq
