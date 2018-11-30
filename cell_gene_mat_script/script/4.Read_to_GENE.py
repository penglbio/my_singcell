import sys
import re

##read the 
for line in open(sys.argv[1]).readlines():
#	print(line)
	[read,strand,xf,gf,gn,gs]=line.strip().split("\t")
	gene={}
	gf_list=gf.split(",")
	gn_list=gn.split(",")
	gs_list=gs.split(",")
	###gene={A:[exon,intron] B:[exon]}
	for idx in range(0,len(gs_list)):
	#	print(gs_list[idx])
		if gs_list[idx]==strand:
			#print(gs_list[idx])
			g=gn_list[idx]
			#print(g)
			if g in gene.keys():
				gene[g].add(gf_list[idx])
			else:
				gene[g]={gf_list[idx]}
	#print(gene)	
	#xf:coding A:coding B:intron remove B from gene list	
	gene1={}
	for gk,gv in gene.items():
		if xf in gv:
			gene1[gk]=gv
	#one Funtion to one gene xf:conding A:coding
	if len(gene1)==1:
		print("%s\t%s" %(read,list(gene1.keys())[0]))
	#one function to >2 genes	
	else:
		GF=sorted(gene1.items(),key=lambda d:d[1])
		#all genes have same function group 
		if GF[0][1]==GF[1][1]:
			print("%s\tambiguous"%read)
		else:
			#xf:coding GF=[{exon},{exon,intron(intergenic)}] or [{intron},{intron,intergenic}] or [{exon},{exon,intron,intergic}] or [{exon,intron},{exon,intron,intergenic} or []]
			if len(GF[0][1])< len(GF[1][1]):
				print("%s\t%s"%(read,GF[0][0]))
			elif len(GF[0][1])== len(GF[1][1]) and "intergenic" in GF[0][1]:
				print("%s\t%s"%(read,GF[0][0]))
			else:
				print("%s\t%s"%(read,GF[1][0]))
			

					
				
