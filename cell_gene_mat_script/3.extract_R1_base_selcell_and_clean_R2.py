import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
bc={}
sel_R1=open(sys.argv[4],"w")
##read whitelist
sel_bar={}
for bar in open(sys.argv[1]).readlines():
	barseq=bar.strip().split("\t\t")[0]
	sel_bar[barseq]=1	

##read R2_all_clean_bar_umi
for record in SeqIO.parse(open(sys.argv[2]), "fastq"):
	if str(record.seq[10:34]) in sel_bar.keys():
		id=re.sub(".2$",".1",record.id)
		bc[str(id)]=record.seq[0:10]+"_"+record.seq[10:34]
print("R2 save complete !")
##selected R1_all add barcode and umi to id
for R1 in SeqIO.parse(open(sys.argv[3]), "fastq"):
	if R1.id in bc.keys():
		R1.id=str(R1.id)+"_"+str(bc[R1.id])
		SeqIO.write(R1, sel_R1, "fastq")
