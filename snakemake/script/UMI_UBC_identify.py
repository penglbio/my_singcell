import sys
import re
index=0
maf={}
reads={}
for line in open(sys.argv[1]).readlines():
	index+=1
	reads[index]=line
	gene=line.split()[10]
	maf[index]=gene
UBC_UMI_GENE=open(sys.argv[3], "w")

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
	GENE_list=re.sub(r"\[|\(|\]|\)","",str(UMI_sort)).replace("dict_items","")
	UBC_UMI_GENE.write(id+"\t"+GENE_list+"\n")
