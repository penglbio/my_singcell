import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
seq=[]
for b in open(sys.argv[1]).read().split('>')[1:]:
	seq.append(b.split("\n")[1])
for record in SeqIO.parse(open(sys.argv[2]), "fastq"):
	[b3,b2,b1]=[str(record.seq[10:18]),str(record.seq[48:56]),str(record.seq[86:94])]
	if b1 in seq and b2 in seq and b3 in seq:
		print(b1+b2+b3)
		
