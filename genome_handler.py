import os
dir_name="ref_dir"
output=open("stats.txt","a")
entries = os.listdir(dir_name+"/")
entries.sort()
for entry in entries:
	i=0
	A_count=0
	C_count=0
	G_count=0
	T_count=0
	N_count=0
	space_count=0
	lent=0
	with open(dir_name+"/"+entry) as f:
		for line in f:
			i=i+1
			if line[0]!='>':
				#lent=lent+len(line)
				A_count=A_count+line.count("A")
				C_count=C_count+line.count("C")
				G_count=G_count+line.count("G")
				T_count=T_count+line.count("T")
				N_count=N_count+line.count("N")
				#space_count=space_count+line.count("\n")
				#s=set(line)
				#print(s)

	lent=A_count+C_count+T_count+G_count+N_count
	out_string=str(A_count)+" "+str(C_count)+" "+str(G_count)+" "+str(T_count)+" "+str(lent)+"\n"
	output.write(entry+"\n")
	output.write(out_string)
	print(entry,"done")
