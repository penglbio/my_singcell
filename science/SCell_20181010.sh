#mkdir 0.raw_data 1.fastqc 2.map_last6 3.clean_UMI 4.clean_barcode_all_split script 5.extract_read1 6.map 7.sortmap stat figure

#seqkit subseq -r 1:94 0.raw_data/A1_20180919_CAGATCAT_S0_L006_R2_001.fastq.gz >0.raw_data/A1_R2.fastq &
#seqkit subseq -r 11:18 0.raw_data/A1_20180919_CAGATCAT_S0_L006_R2_001.fastq >0.raw_data/b3.fq &
#seqkit subseq -r 19:48 0.raw_data/A1_20180919_CAGATCAT_S0_L006_R2_001.fastq >0.raw_data/b3_to_b2.fq &
#seqkit subseq -r 49:56 0.raw_data/A1_20180919_CAGATCAT_S0_L006_R2_001.fastq >0.raw_data/b2.fq &
#seqkit subseq -r 57:86 0.raw_data/A1_20180919_CAGATCAT_S0_L006_R2_001.fastq >0.raw_data/b2_to_b1.fq &
#seqkit subseq -r 87:94 0.raw_data/A1_20180919_CAGATCAT_S0_L006_R2_001.fastq >0.raw_data/b1.fq &

#grep -v '>' 0.barcode/barcode.fa |cut -f 3-8|rush -k 'seqkit grep -s -r -p '{}'$ 0.raw_data/A1_R2.fastq > 2.map_last6/{}.fastq'

#python script/check_QC_barcode_con.py 0.barcode/barcode.fa 0.raw_data/A1_R2.fastq > 0.raw_data/all_barcode.fa &

#ls 2.map_last6/*|rush 'python script/check_QC_UMI.py 3.clean_UMI/{%.}.fq {}' &

#ls 3.clean_UMI/*|rush -k 'python script/check_QC_barcode.py  0.barcode/barcode.fa {} 4.clean_barcode/{%}' &

#cat 3.clean_UMI/* >all.fq
#seqkit split -p 100 all.fq
#ls all.fq.split/*|rush -k 'python script/check_QC_barcode_con.py  0.barcode/barcode.fa {} 4.clean_barcode_all_split/{%}'

#ls 4.clean_barcode_all_split/*.fq|rush -k "python script/exact_replace_id.py {} 0.raw_data/A1_R1.fastq 5.extract_read/{%}"

cat 5.extract_read/all.*.fq >5.extract_read/A1_R1.extract.fastq 

trimmomatic SE -threads 10 -phred33 5.extract_read/A1_R1.extract.fastq 5.extract_read/A1_R1.extract_clean.fastq  ILLUMINACLIP:adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

fastqc 5.extract_read/A1_R1.extract_trimmed.fastq -t 30

cutadapt -a "NNNNNNNNCGAATGCTCTGGCCT;max_error_rate=0.1;min_overlap=23" --minimum-length 36 5.extract_read/A1_R1.extract_clean.fastq -o 5.extract_read/A1_R1.extract_trim_3adapter.fastq -j 30 1>5.extract_read/cutadapt_log_1 &
##remove rRNA sortmerna
cd ~/projects/rRNA/rRNA_databases/
indexdb_rna --ref rfam-5.8s-database-id98.fasta,rfam-5.8s-database-id98.idx:rfam-5s-database-id98.fasta,rfam-5s-database-id98.idx:silva-arc-16s-id95.fasta,silva-arc-16s-id95.idx:silva-arc-23s-id98.fasta,silva-arc-23s-id98.idx:silva-bac-16s-id90.fasta,silva-bac-16s-id90.idx:silva-bac-23s-id98.fasta,silva-bac-23s-id98.idx:silva-euk-18s-id95.fasta,silva-euk-18s-id95.idx:silva-euk-28s-id98.fasta,silva-euk-28s-id98.idx  -m 10240 -v

sortmerna --ref rfam-5.8s-database-id98.fasta,rfam-5.8s-database-id98.idx:rfam-5s-database-id98.fasta,rfam-5s-database-id98.idx:silva-arc-16s-id95.fasta,silva-arc-16s-id95.idx:silva-arc-23s-id98.fasta,silva-arc-23s-id98.idx:silva-bac-16s-id90.fasta,silva-bac-16s-id90.idx:silva-bac-23s-id98.fasta,silva-bac-23s-id98.idx:silva-euk-18s-id95.fasta,silva-euk-18s-id95.idx:silva-euk-28s-id98.fasta,silva-euk-28s-id98.idx -a 60 --aligned ../../../04_JSK_SC_20181024/6.map/unmapped_reads_rrna_map.sam --fastx --reads ../../../04_JSK_SC_20181024/6.map/mix_map150_cleanUnmapped.cutadapt4 -m 1032094


##star 
STAR --runThreadN 30 --genomeDir ~/projects/genome/hg38_mm10/hg38_mm10_index_hang150/ --readFilesIn 5.extract_read/A1_R1.extract_clean.fastq --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --outReadsUnmapped Fastx --outFileNamePrefix 6.map/mix_map150_clean > 6.map/star_log.txt
(Scell) [lpeng@node2 04_JSK_SC_20181024]$ STAR --runThreadN 30 --genomeDir ~/projects/genome/hg38_mm10/hg38_mm10_index_hang150/ --readFilesIn 5.extract_read/A1_R1.extract_trim_3adapter1.fastq --outFilterScoreMinOverLread 0.33 --outFilterMatchNminOverLread 0.33 --outReadsUnmapped Fastx  --outFileNamePrefix 6.map/mix_map150_clean 1>6.map/star_log.txt &
#bowtie
bowtie2 -x /cluster/home/lpeng/projects/genome/hg38_mm10/hg38_mm10_index_bowtie/hg38_mm10 -U 5.extract_read/A1_R1.extract.fastq -S 6.map_bowtie_2/mix_mapped.sam -p 60 --phred33 --reorder &

samtools sort 6.map/mix_map_hang150Aligned.out.sam -@ 60 -o 7.sortmap/mix_map_hang150Aligned.out.sort.bam & 

mv 7.sortmap/mix_map_hang10Aligned.out.sort.bam 7.sortmap/mix_mapAligned.out.sort.bam

TagReadWithGeneExon I=7.sortmap/mix_mapAligned.out.sort.bam O=7.sortmap/Taged.bam SUMMARY=7.sortmap/sum.File ANNOTATIONS_FILE=/cluster/home/lpeng/projects/genome/hg38_mm10/hg38_mm10_gene_name_replaced_codgene.94.gtf 

samtools view -@ 30 -F 256 -h 7.sortmap/Taged.bam -b |bedtools bamtobed -i - >7.sortmap/Taged_primary.bed 

samtools view -@ 30 -F 256 -h 7.sortmap/Taged.bam -b >7.sortmap/Taged_primary.bam

samtools view -@ 30 7.sortmap/Taged_primary.bam |cut -f12 >7.sortmap/TAG_list.txt

paste 7.sortmap/Taged_primary.bed 7.sortmap/TAG_list.txt >7.sortmap/Taged_primary_TAG.bed

bedtools intersect -a 7.sortmap/Taged_primary_TAG.bed -b ~/projects/genome/hg38_mm10/hg38_mm10_gene_list_6.bed -s -wa -wb -f 0.9 >7.sortmap/Taged_primary_TAG_gene_list_6.bed

#python script/barcodeID_to_seq.py 0.barcode/barcode.fa 7.sortmap/Taged_annote.bed >7.sortmap/Taged_annote_UBC_trans.bed

cut -f3 -d"_" 7.sortmap/Taged_annote.bed|sed -E 's/\t.*\t//g' >7.sortmap/UBC+UMI &

starcode-v1_3-x86_64_linux -i 7.sortmap/UBC+UMI -o 7.sortmap/UBC+UMI_cluster -q -t 30 -s -d0 --seq-id &
starcode-v1_3-x86_64_linux -i 7.sortmap/UBC+UMI -o 7.sortmap/UBC+UMI_cluster1 -q -t 30 -s -d0 --seq-id --print-cluters &
starcode-v1_3-x86_64_linux -i 7.sortmap/UBC+UMI -o 7.sortmap/UBC+UMI_cluster2 -q -t 30 -s -d1 --seq-id &
starcode-v1_3-x86_64_linux -i 7.sortmap/UBC+UMI -o 7.sortmap/UBC+UMI_cluster3 -q -t 30 -s -d1 --seq-id --print-clusters &

python script/UMI_UBC_identify.py 7.sortmap/Taged_annote.bed 7.sortmap/UBC+UMI_cluster1 7.sortmap/Taged_annote_UBC_trans_most+freq.bed 7.sortmap/UM_UBC_clustergene_list 7.sortmap/UMI_UBC_keep_gene_cluster &
##stat
#python script/check_QC_barcode_con_sepa.py 0.barcode/barcode.fa 0.raw_data/A1_R2.fastq 86 94 >stat/b1_barcode_check_ratio.txt
#python script/check_QC_barcode_con_sepa.py 0.barcode/barcode.fa 0.raw_data/A1_R2.fastq 48 56 >stat/b2_barcode_check_ratio.txt
#python script/check_QC_barcode_con_sepa.py 0.barcode/barcode.fa 0.raw_data/A1_R2.fastq 10 18 >stat/b3_barcode_check_ratio.txt
#python script/check_QC_barcode_con_sepa.py 0.barcode/insert.fa 0.raw_data/A1_R2.fastq 56 86 >stat/b1_b2_barcode_check_ratio.txt &
#python script/check_QC_barcode_con_sepa.py 0.barcode/insert.fa 0.raw_data/A1_R2.fastq 18 48 >stat/b2_b3_barcode_check_ratio.txt &

#seqkit stat 2.map_last6/* >stat/last6_stat.tsv
#seqkit stat 3.clean_UMI/* >stat/clean_UMI_stat.tsv &
seqkit stat 4.clean_barcode_all_split/*.fq >stat/clean_barcode_stat.tsv &
less -S stat/clean_barcode_stat.tsv|sed -E 's/ +/\t/g'|cut -f 4|sed 's/,//g'|sed '1d'|awk '{sum += $1};END {print sum}'
starcode-v1_3-x86_64_linux -i 7.sortmap/UBC+UMI -o 7.sortmap/UBC+UMI_cluster -q -t 30 -s -d0 --seq-id &
starcode-v1_3-x86_64_linux -i 7.sortmap/UBC+UMI -o 7.sortmap/UBC+UMI_cluster1 -q -t 30 -s -d0 --seq-id --print-cluters &
starcode-v1_3-x86_64_linux -i 7.sortmap/UBC+UMI -o 7.sortmap/UBC+UMI_cluster2 -q -t 30 -s -d1 --seq-id &
starcode-v1_3-x86_64_linux -i 7.sortmap/UBC+UMI -o 7.sortmap/UBC+UMI_cluster3 -q -t 30 -s -d1 --seq-id --print-clusters &

seqkit seq -n 5.extract_read1/A1_R1.extract.fastq | cut -f 3 -d'_' > stat/bar_raw
sort stat/bar_raw|uniq -c |sort -k1nr>stat/bar_uniq
less stat/bar_uniq |head -4000|awk '$1>10000' >stat/barcode_select

