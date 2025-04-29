# Garret Miller
# Python script for collecting counts and TPM data from salmon outputs.

import sys, os

expMatrixPrefix = str(sys.argv[1])

transcript_TPM_Counts = {}
header_switch = "on"

# Second argument and subsequent arguments are salmon output directories
quantFiles = sys.argv[2:]

# Collects expression data from each salmon output
def get_expression(input):
	header_switch = "on"
	for i in input:
		if header_switch == "on":
			header_switch = "off"
			continue
		line = i.strip().split()
		if line[0] not in transcript_TPM_Counts.keys():
			transcript_TPM_Counts[line[0]] = [[line[3]], [line[4]]]
		else:
			transcript_TPM_Counts[line[0]][0].append(line[3])
			transcript_TPM_Counts[line[0]][1].append(line[4])
	header_switch = "on"
for i in quantFiles:
	with open(i) as individualFile:
		get_expression(individualFile)

# Writes contents of transcript_TPM_Counts to separate matrices for Counts and TPM values
for fileType in ["_expressionTPM", "_expressionCounts"]:
	output = open(expMatrixPrefix+fileType+".out", "w")
	output.write("Transcript")
	for i in quantFiles:
		sampleName = i.split(".")[1]
		output.write("\t"+sampleName)
	output.write("\n")
	for i in transcript_TPM_Counts.keys():
		output.write(i)
		if fileType == "_expressionTPM":
			for j in transcript_TPM_Counts[i][0]:
				output.write("\t"+str(j))
			output.write("\n")
		elif fileType == "_expressionCounts":
			for j in transcript_TPM_Counts[i][1]:
				output.write("\t"+str(j))
			output.write("\n")
	output.close()