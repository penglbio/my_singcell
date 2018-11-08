import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
seq=[]
count=num_b1=0
for b in open(sys.argv[1]).read().split('>')[1:]:
	seq.append(b.split("\n")[1])
for record in SeqIO.parse(open(sys.argv[2]), "fastq"):
	b1=str(record.seq)
	count+=1
	if b1 in seq:
		num_b1+=1
	
print(count,num_b1)	
