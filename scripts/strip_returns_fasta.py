"""
Garret Miller
2/5/2019

Strips returns from the middle of sequences in a fasta file.
If the fasta file has extra info after ">#######", you can retain it by typing "retain" as a second argument

Run program like this:
	python strip_returns_fasta.py <name of fasta>
It will retain the original file and output another file called:
	<name of fasta>_stripped
"""
import sys

input = open(sys.argv[1])
output = open(sys.argv[1]+"_stripped","w")

try:
	if sys.argv[2].lower() == "retain":
		keep_info = "yes"
		print("Retaining everything in the title line.")
	else:
		keep_info = "no"
except:
	keep_info = "no"
	print("If your title line has more info after a tab or space, it will be removed in the output.")
	print("You can keep this by entering <retain> as a second argument.")

first_line_switch = "on"

if keep_info == "yes":
	for i in input:
		line = i.strip()
		if first_line_switch == "on":
			output.write(line+"\n")
			first_line_switch = "off"
		elif ">" in line:
			output.write("\n"+line+"\n")
		else:
			output.write(line)
elif keep_info == "no":
	for i in input:
		line = i.strip()
		line_info = i.strip().split(" ")
		if first_line_switch == "on":
			output.write(line_info[0]+"\n")
			first_line_switch = "off"
		elif ">" in line:
			output.write("\n"+line_info[0]+"\n")
		else:
			output.write(line)

input.close()
output.close()
