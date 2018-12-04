import sys
seq=open(sys.argv[1],"w")
j=["A","T","C","G"]
seq.write(">seq\n")
for a in j:
	for b in j:
		for c in j:
			for d in j:
				for e in j:
					for f in j:
						for g in j:
							for h in j:
								seq.write(a+b+c+d+e+f+g+h+"*")
