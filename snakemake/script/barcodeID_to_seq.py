import sys
from Bio import SeqIO
from Bio.Seq import Seq
import re
BC={}
for b in open(sys.argv[1]).read().split('>')[1:]:
	id=b.split("\n")[0]
	seq=b.split("\n")[1]
	BC[id]=seq
#BCR=open(sys.argv[3],"w")
for line in open(sys.argv[2]).readlines():
	[id1,id2,id3]=line.strip().split("_")[1].split("+")
	print(line.strip()+"\t"+str(BC[id1])+str(BC[id2])+str(BC[id3]))
