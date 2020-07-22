import os
dir_name="ref_dir"
out_dir_name="ref_dir_"
entries = os.listdir(dir_name+"/")
entries.sort()
for entry in entries:
	with open(dir_name+"/"+entry, 'r') as infile, \
		 open(out_dir_name+"/"+entry, 'w') as outfile:
		data = infile.read()
		data = data.replace("N", "")
		data = data.replace("a", "A")
		data = data.replace("c", "C")
		data = data.replace("t", "T")
		data = data.replace("g", "G")
		outfile.write(data)
	print(entry + " done.")
