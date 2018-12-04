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
count=0
seq=""
bs=""
with open(sys.argv[2], 'r') as f:
	for line in f:
		count=count+1
		if count%4==1:
			fid=line.split()[0]
		elif count%4==2:
			bs=line[10:34]
			seq=line[0:10]+"_"+line[10:34]
		if bs in sel_bar.keys():	
			bc[str(fid)]=seq
print("R2 save complete !")
##selected R1_all add barcode and umi to id
n=4
with open(sys.argv[3], 'r') as fh:
	lines = []
	for line in fh:
		lines.append(line.strip())
		if len(lines) == n:
			r1_id=lines[0].split()[0]
			if r1_id in bc.keys():
				lines[0]=r1_id+"_"+bc[r1_id]
				for i in lines:
					sel_R1.write(i+"\n")
			lines = []
