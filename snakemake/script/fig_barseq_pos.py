import sys

b_pos=open(sys.argv[1])
#print("*"*4**8*9)
for line in b_pos.readlines():
	[start,end]=line.strip().split("\t")[4:6]
	print("-"*int(start),"*"*(int(end)-int(start)))
