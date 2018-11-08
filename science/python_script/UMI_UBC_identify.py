import sys
import re
index=0
maf={}
reads={}
for line in open(sys.argv[1]).readlines():
	index+=1
	reads[index]=line
	gene=line.split()[9]
	maf[index]=gene
READS_KEEP_INX=open(sys.argv[3], "w")
UBC_UMI_GENE=open(sys.argv[4], "w")
UBC_UMI_GENE_INX=open(sys.argv[5], "w")

for umi_c in open(sys.argv[2]).readlines():
	UMI={}
	gene_idx={}
	gene_list=[]
	id=umi_c.split("\t")[0]
	cluster=umi_c.split("\t")[2]
	for inx in cluster.strip().split(","):
		gene=maf[int(inx)]
		if(gene in gene_list):
			UMI[gene]+=1
			gene_idx[gene].append(int(inx))	
		else:
			UMI[gene]=1
			gene_list.append(gene)
			gene_idx[gene]=[int(inx)]
	UMI_sort=sorted(UMI.items(),key=lambda item:item[1],reverse=True)
	Most_assign_gene=UMI_sort[0][0]
	Most_assign_gene_str=str(gene_idx[Most_assign_gene]).replace("[","").replace("]","")	
	#UBC_UMI_GENE_INX.write(umi_c.strip()+"\t"+Most_assign_gene+"\t"+str(gene_idx[Most_assign_gene])+"\n")
	UBC_UMI_GENE_INX.write(umi_c.strip()+"\t"+Most_assign_gene+"\t"+Most_assign_gene_str+"\n")
	GENE_list=re.sub(r"\[|\(|\]|\)","",str(UMI_sort)).replace("dict_items","")
	UBC_UMI_GENE.write(id+"\t"+GENE_list+"\n")
	for idx in gene_idx[Most_assign_gene]:
		READS_KEEP_INX.write(reads[idx])
