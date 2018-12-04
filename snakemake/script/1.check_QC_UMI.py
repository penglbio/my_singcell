import sys
from Bio import SeqIO
from collections import Counter
UMI = open(sys.argv[1], "w")
for record in SeqIO.parse(open(sys.argv[2]), "fastq"):
	count=0
	for sr in record.letter_annotations["phred_quality"][:10]:
		if(sr<10):
			count +=1
			if(count>1):
				break
	if(count<1):
		SeqIO.write(record, UMI, "fastq")	
