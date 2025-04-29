# Garret Miller
# Python script meant to take a clustered transcriptome and pull out all associated sequences from a translated peptide file
# This works by making a set of all sequence names present in the clustered transcriptome (First index in title line) and searching for these names in the translated peptide file
# The string stored in the set above needs to be present within the title lines in the translated peptide file (this splits by a "." to remove a ".p1" added by TransDecoder)

import sys

# Input your transcriptome clustered by CD-Hit as the first argument, .pep file from TransDecoder as the second, and give the name of an output file as the third
clusteredTranscriptome = open(sys.argv[1])
translatedPEP = open(sys.argv[2])
outFile = open(sys.argv[3],"w")

# Used for writing sequences when not iterating over a title line
switch = "off"
# Set for storing sequence names present in clustered transcriptome
seqs = set()

# Iterates over the clustered transcriptome "clusteredTranscriptome" to pick out sequence names
for i in clusteredTranscriptome:
	if ">" in i:
		line = i.strip().split()
		seqName = line[0][1:]
		seqs.add(seqName)

# Iterates over the translated peptide file "translatedPEP" to write all sequences present in the set "seqs" to a new file
for i in translatedPEP:
	if ">" in i:
		line = i.strip().split()
		seqName = line[0][1:]
		seqName = seqName.split(".")[0]
		if seqName in seqs:
			outFile.write(i)
			switch = "on"
	elif switch == "on":
		outFile.write(i)
		switch = "off"

clusteredTranscriptome.close()
translatedPEP.close()
outFile.close()
