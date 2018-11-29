
prefetch -v SRR6750056.sra &
prefetch -v SRR6750057.sra &
fastq-dump -B -I --split-spot --split-files --skip-technical --gzip SRR6750056.sra &
fastq-dump -B -I --split-spot --split-files --skip-technical --gzip SRR6750057.sra &
