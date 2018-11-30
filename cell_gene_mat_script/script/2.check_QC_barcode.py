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
	BCq={}
	for b in open(sys.argv[1]).read().split('>')[1:]:
		id=b.split("\n")[0]
		seq=b.split("\n")[1]
		BC[seq]=id
		BCq[id]=seq
	BCR=open(sys.argv[3],"w")
	for record in SeqIO.parse(open(sys.argv[2]), "fastq"):
		b1q=record[86:94]
		b2q=record[48:56]
		b3q=record[10:18]
		umi=record[0:10]
		record_clean=umi+b1q+b2q+b3q	
		b1qc=b2qc=b3qc=None
		for bseq in BC.keys():
			if(b1qc is None):
				b1qc=barcode_mismatch_check(b1q.seq,0,bseq)
				b1id=BC[bseq]
			if(b2qc is None):	
				b2qc=barcode_mismatch_check(b2q.seq,0,bseq)
				b2id=BC[bseq]
			if(b3qc is None):
				b3qc=barcode_mismatch_check(b3q.seq,0,bseq)
				b3id=BC[bseq]
		if b1qc is not None and b2qc is not None and b3qc is not None:
			if int(b1id) >48:
				recor_clean.id=record.id+str(b3id)+str(b2id)+str(int(b1id)-48))
				record_clean.seq=umi.seq+b3qc+b2qc+BCq[str(int(b1id)-48)]
			else:
				recor_clean.id=record.id+str(b3id)+str(b2id)+str(b1id)
				record_clean.seq=umi.seq+b3qc+b2qc+b1qc
			SeqIO.write(record_clean, BCR, "fastq")
