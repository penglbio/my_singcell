import sys
list_ubc=open(sys.argv[1]).read().split("\n")[:-1]
#print(list_ubc)
k_ubc_r=open(sys.argv[3],"w")
for line in open(sys.argv[2]).readlines():
	id_ubc=line.strip().split("_")[2]
#	print(id_ubc)
	if id_ubc in list_ubc:
		k_ubc_r.write(line)
