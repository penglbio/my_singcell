mkdir step1 step2 step3 step4 step5 step6 step7 step8 step9
#step1 last6
grep -v '>' ../barcode/barcode.fa |cut -f 3-8|rush -k 'seqkit grep -s -r -p '{}'$ raw_data/SRR6750057_2.fastq > step1/{}.fastq'
seqkit stat step1/* > stat/step1_stat.tsv

#step2 phred(UMI)<10 <=1bp
ls step1/*|rush 'python script/1.check_QC_UMI.py step2/{%.}.fq {}' 
seqkit stat step2/* >stat/step2_stat.tsv

#step3 process barcode
ls step2/*|rush -k 'python script/2.check_QC_barcode.py ../barcode/barcode.fa {} step3/{%}' 
seqkit stat step3/* >stat/step3_stat.tsv

#step4 merge split R2(barcode) umi+b3+b2+b1 and select cells
cat step3/* >step4/all_clean_bar_umi.fq
umi_tools whitelist --stdin step4/all_clean_bar_umi.fq --plot-prefix step4/cell_sel --subset-reads 3000000000 --timeit umi_cell_time --bc-pattern=NNNNNNNNNNCCCCCCCCCCCCCCCCCCCCCCCC --set-cell-number=312 --log2stderr > step4/whitelist_312.txt  
#step5 inputfiles whitelist all_clean_bar_umi.fq(R2) raw_R1(R1) seleced R1(change the header)
python script/3.extract_R1_base_selcell_and_clean_R2.py step4/whitelist_312.txt step4/all_clean_bar_umi.fq step5/SRR6750056_1_extract_ubc.fastq

cutadapt -a "NNNNNNNNCGAATGCTCTGGCCT;max_error_rate=0.1;min_overlap=23"  -a "AAAAAAAAAAAAAAAAA;max_error_rate=0.1;min_overlap=6"  step5/SRR6750056_1_extract_ubc.fastq --minimum-length 36  -o 5.step5/SRR6750056_1_extract_ubc_cutadapter.fastq 

trimmomatic SE -threads 10 -phred33 step5/SRR6750056_1_extract_ubc_cutadapter.fastq step5/SRR6750056_1_extract_ubc_clean.fastq  ILLUMINACLIP:adapter.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36


##step6 STAR map
STAR --runThreadN 30 --genomeDir ~/projects/genome/hg38_mm10/hg38_mm10_index_hang65/ --readFilesIn step5/SRR6750056_1_extract_ubc.fastq  --outFileNamePrefix step6/map_ --outSAMtype BAM SortedByCoordinate &

samtools view -F 256 step6/map_Aligned.sortedByCoord.out.bam -h -b>step6/map_primary.bam

TagReadWithGeneFunction INPUT=step6/map_primary.bam OUTPUT=step6/Tag.bam ANNOTATIONS_FILE= ~/projects/genome/hg38_mm10/hg38_mm10_gene_name_replaced.94.gtf

samtools view step6/Tag.bam|grep -v "XF:Z:INTERGENIC" |awk -v OFS="\t" '{if($2==16)print $1,"-",$12,$17,$18,$19;else print $1,"+",$12,$17,$18,$19}'|sed 's/..:Z://g'|sed 's/UTR/CODING/g' >step6/Tag_Reads_gene.tsv

python script/4.Read_to_GENE.py step6/Tag_Reads_gene.tsv > step6/read_to_gene.tsv

###step7 cell gene matrix
cut -f2,4 -d"_" step6/read_to_gene.tsv|sed -E 's/\t.*//g' |sed 's/_//g' > step7/UMI+UBC_seq

starcode-v1_3-x86_64_linux -i step7/UMI+UBC_seq -o step7/UMI_UBC_seq_cluster -q -t 30 -s -d1 --seq-id

python script/5.UMI_UBC_identify.py step6/read_to_gene.tsv step7/UMI_UBC_seq_cluster step7/UMI_UBC_gene

sed -E "s/,/\t/g" step7/UBC+UMI_gene|awk '{if(($3-$5)>0){print $0}}' |awk '{print $1"\t"$2}'> step7/UMI_UBC_gene1

python script/6.barcode_seq_to_id.py barcode/barcode.fa step7/UMI_UBC_gene1 step7/UMI_UBC_gene2

awk -v OFS="\t" '{print $1,$3,$2}' step7/UMI_UBC_gene2  | sort -k1 -k2 | bedtools groupby -g 1,2 -c 3 -o count > step7/cell_gene_count.tsv

Rscript script/df_to_tbl.R step7/cell_gene_count.tsv step7/cell_gene_mat.tsv

awk -v OFS="\t" '{print $1,$3,$2}' step7/UMI_UBC_gene2  | sort -k1 -k2 | bedtools groupby -g 1,2 -c 3 -o count > step7/cell_gene_count.tsv
