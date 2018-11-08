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
for record in SeqIO.parse(open(sys.argv[2]), "fastq"):
		[b3q,b2q,b1q]=[str(record.seq[10:18]),str(record.seq[48:56]),str(record.seq[86:94])]
		if b1q in BC.keys() and b2q in BC.keys() and b3q in BC.keys():
			b1id=BC[b1q]
			b2id=BC[b2q]
			b3id=BC[b3q]
			if int(b1id) >48:
				record.id=record.id+"_"+str(b1q)+str(b2q)+str(b3q)+"_"+str(int(b1id)-48)+"+"+str(b2id)+"+"+str(b3id)+"_"+str(record.seq[:10])
			else:
				record.id=record.id+"_"+str(b1q)+str(b2q)+str(b3q)+"_"+str(b1id)+"+"+str(b2id)+"+"+str(b3id)+"_"+str(record.seq[:10])
			SeqIO.write(record, BCR, "fastq")
