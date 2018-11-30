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
	b3id=BC[(line[10:18])]
	b2id=BC[line[18:26]]
	b1id=BC[line[26:34]]
	if int(b1id) >48:
		line=str((int(b1id)-48))+"+"+b2id+"+"+b3id+"\t"+line[0:10]+"\t"+line[34:-1]
	else:
		line=b1id+"+"+b2id+"+"+b3id+"\t"+line[0:10]+"\t"+line[34:-1]
	BCR.write(line+"\n")
