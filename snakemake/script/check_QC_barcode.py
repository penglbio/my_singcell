def barcode_mismatch_check(genome, posn, seq):
	errors = 0
	for i in range(8):
		if genome[posn+i] != seq[i]:
			errors += 1
			if errors >= 2:
				break
	if(errors<=1):
		return(seq) 

if __name__ == '__main__':
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
		b1q=None
		b2q=None
		b3q=None
		for bseq in BC.keys():
			if(b1q is None):
				b1q=barcode_mismatch_check(record.seq,86,bseq)
				b1id=BC[bseq]
			if(b2q is None):	
				b2q=barcode_mismatch_check(record.seq,48,bseq)
				b2id=BC[bseq]
			if(b3q is None):
				b3q=barcode_mismatch_check(record.seq,10,bseq)
				b3id=BC[bseq]
		if b1q is not None and b2q is not None and b3q is not None:
			if int(b1id) >48:
				record.id=record.id+"_"+str(b1q)+str(b2q)+str(b3q)+"_"+str(int(b1id)-48)+"+"+str(b2id)+"+"+str(b3id)+"_"+str(record.seq[:10])
			else:
				record.id=record.id+"_"+str(b1q)+str(b2q)+str(b3q)+"_"+str(b1id)+"+"+str(b2id)+"+"+str(b3id)+"_"+str(record.seq[:10])
			SeqIO.write(record, BCR, "fastq")
