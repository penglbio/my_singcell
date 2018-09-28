def barcode_mismatch_check(genome, posn, sequence):
	errors = 0
	for i in range(8):
		if genome[posn+i] != sequence[i]:
			errors += 1
			if errors >= 2:
				break
	if(errors<=1):
		return("ok") 

if __name__ == '__main__':
	import sys
	from Bio import SeqIO
	from Bio.Seq import Seq
	import re
	b1=sys.argv[1] 
	BCR=open(sys.argv[4],"w")
	for record in SeqIO.parse(open(sys.argv[3]), "fastq"):
		b3_qc=[]
		b2_qc=[]
		b1_qc=[]
		
		b1_qc.append(barcode_mismatch_check(record.seq,86,b1))
		for b in SeqIO.parse(open(sys.argv[2]), "fasta"):
			if re.match("barcode2",b.id):
				b2_qc.append(barcode_mismatch_check(record.seq,48,b.seq))
			else:
				b3_qc.append(barcode_mismatch_check(record.seq,10,b.seq))
		if "ok" in b1_qc and "ok" in b2_qc and "ok" in b3_qc:
			SeqIO.write(record, BCR, "fastq")
