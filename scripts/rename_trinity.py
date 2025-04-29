# Garret Miller
# Python script for renaming Trinity outputs with a specified string

import sys

input = open(sys.argv[1])
code = str(sys.argv[2])
output = open(code+"_"+sys.argv[1],"w")

for i in input:
	if ">" in i:
		output.write(">"+code+i[8:])
	else:
		output.write(i)

input.close()
output.close()
