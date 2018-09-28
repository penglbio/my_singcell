cat ../barcode/barcode1.fa |grep -v '>'|cut -f 3-8|rush -k 'seqkit grep -s -r -p '{}'$ raw_data/SRR6750057_2.fastq > step1/{}.fastq' &
ls step1/*|rush 'python check_QC_UMI.py step2/{%.}.fq {}' &
ls step2/*|rush -k 'python check_QC_barcode.py {%.} ../barcode/barcode_2_3.fa {} step3/{%}' --dry-run
STAR --runThreadN 30 --genomeDir ../mm10_and_hg19/hg38_mm10/ --readFilesIn  raw_data/SRR6750057_1.fastq --outFileNamePrefix  step5/STAR_to_genome &
