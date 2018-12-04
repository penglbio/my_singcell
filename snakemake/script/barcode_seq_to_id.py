import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
BC={}
for b in open(sys.argv[1]).read().split('>')[1:]:
	id=b.split("\n")[0]
	seq=b.split("\n")[1]
	BC[seq]=id
BCR=open(sys.argv[3],"w")
for line in open(sys.argv[2]).readlines():
	b1id=BC[(line[0:8])]
	b2seq=line[8:16]
	b3seq=line[16:24]
	if int(b1id) >48:
		line=b3seq+b2seq+"_"+str((int(b1id)-48))+"\t"+line[24:-1]
	else:
		line=b3seq+b2seq+"_"+str(b1id)+"\t"+line[24:-1]
	BCR.write(line+"\n")
