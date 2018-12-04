import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
start=int(sys.argv[3])
end=int(int(sys.argv[4]))
seq=[]
total=num_b=num_Nr=num_errb=0
for b in open(sys.argv[1]).read().split('>')[1:]:
	seq.append(b.split("\n")[1])
for record in SeqIO.parse(open(sys.argv[2]), "fastq"):
	total+=1
	bc_seq=str(record.seq[start:end])
	if bc_seq in seq:
		num_b+=1
	elif re.match("N",bc_seq):
		num_Nr+=1
	else:
		num_errb+=1
print("%d=%d+%d+%d" %(total,num_b,num_Nr,num_errb))	
