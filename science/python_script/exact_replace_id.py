import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
bc={}
sel_R1=open(sys.argv[3],"w")
for record in SeqIO.parse(open(sys.argv[1]), "fastq"):
	[id,seq]=record.id.split("_",1)
	id=id.replace(".2",".1")
	bc[id]=seq

for R1 in SeqIO.parse(open(sys.argv[2]), "fastq"):
	if R1.id in bc.keys():
		R1.id=R1.id+"_"+bc[R1.id]
		SeqIO.write(R1, sel_R1, "fastq")
