params.assemblyReads = "Error: Input the name of your fastq file (pair) with --assemblyReads <fileName_{1,2}.fastq>. Ensure that paired data includes the '{1,2}' as shown in the example."
params.mappingReads = "Error: Input the directory containing individual read fastq files for mapping with --mappingReads <read directory>"
params.CODE = "Error: Input a string to be added as a prefix in the assembly and output files with --CODE <code>"
params.scriptDir = "Error: Input the directory containing all necessary scripts for this pipeline with --scriptDir <directory>"

log.info """\
	Values of Required Parameters:
	======================================================================
	assemblyReads: ${params.assemblyReads}
	mappingReads : ${params.mappingReads}
	CODE         : ${params.CODE}
	scriptDir    : ${params.scriptDir}
	======================================================================
	"""
	.stripIndent(true)

// --- General Parameters ---------------------------------------------------
params.outDir = 'outputs'
params.threads = 48
params.memory = 240

// --- Trinity Parameters ---------------------------------------------------
params.trinityOutDir = 'trinity_out'
params.trinityThreads = '48'
params.trinityMemory = params.memory + 'G'
params.bflyHeapSpaceMax = '20G'

// --- TransDecoder Parameters ----------------------------------------------
params.transdecoderOutDir = 'transdecoder_out'
params.transdecoderTemp = 'transdecoder_temp'

// --- CD-HIT Parameters ----------------------------------------------------
params.CDHoutDir = 'cdhit'
params.CDHdescLen = '9999'
params.CDHpercentID = 0.99
params.CDHwordSize = '11'
params.CDHmemory = params.memory * 1000
params.CDHthreads = params.threads
params.CDHalnType = '0'
params.CDHrepAlnLen = '0.15'
params.CDHpercentIDstring = (params.CDHpercentID * 100).toInteger()
params.CDHoutput = "${params.CDHoutDir}/${params.CODE}_Trinity_clustered${params.CDHpercentIDstring}.fasta"

// --- Match Fastas Parameters ----------------------------------------------
params.matchFastasOutput = "${params.CDHoutDir}/${params.CODE}_Trinity.fasta.transdecoder_clustered${params.CDHpercentIDstring}-byCDS.fasta"

// --- Salmon Parameters ----------------------------------------------------
params.salmonOutDir = 'salmon'
params.salmonThreads = params.threads
params.salmonLibType = "A"
params.salmonIndexOutput = "${params.salmonOutDir}/${params.CODE}_index"
params.salmonQuantsDir = "${params.salmonOutDir}/quants/"
params.salmonOutTPM = "${params.CODE}_expressionTPM.out"
params.salmonOutCounts =  "${params.CODE}_expressionCounts.out"




// Assembles a transcriptome from specified fastq file pairs
process TRINITY {
	publishDir "${params.outDir}", mode: 'copy'
	input:
		tuple val(name), path(concatReads)
	output:
		path "${params.trinityOutDir}/${params.CODE}_Trinity.fasta"
	script:
		"""
		# Sets the Java heap max for the Butterfly step of Trinity
		export _JAVA_OPTIONS="-Xmx${params.bflyHeapSpaceMax}"

		# Runs Trinity
		Trinity --seqType fq --left ${concatReads[0]} --right ${concatReads[1]} --CPU ${params.trinityThreads} --output ${params.trinityOutDir} --max_memory ${params.trinityMemory} --bflyHeapSpaceMax ${params.bflyHeapSpaceMax}
		
		# Manages output files
		python ${params.scriptDir}/strip_returns_fasta.py "${params.trinityOutDir}/Trinity.fasta"
		mv ${params.trinityOutDir}/Trinity.fasta_stripped ./
		python ${params.scriptDir}/rename_trinity.py Trinity.fasta_stripped ${params.CODE}
		mv ${params.CODE}_Trinity.fasta_stripped "${params.trinityOutDir}/${params.CODE}_Trinity.fasta"
		rm Trinity.fasta_stripped
		"""
}

// Assembles a transcriptome from a specified unpaired fastq file
// This is an alternative for the TRINITY process above when working with unparied data
process TRINITY_SINGLE {
	publishDir "${params.outDir}", mode: 'copy'
	input:
		path concatReads 
	output:
		path "${params.trinityOutDir}/${params.CODE}_Trinity.fasta"
	script:
		"""
		# Sets the Java heap max for the Butterfly step of Trinity
		export _JAVA_OPTIONS="-Xmx${params.bflyHeapSpaceMax}"

		# Runs Trinity
		Trinity --seqType fq --single ${concatReads} --CPU ${params.trinityThreads} --output ${params.trinityOutDir} --max_memory ${params.trinityMemory} --bflyHeapSpaceMax ${params.bflyHeapSpaceMax}
		
		# Manages output files
		python ${params.scriptDir}/strip_returns_fasta.py "${params.trinityOutDir}/Trinity.fasta"
		mv ${params.trinityOutDir}/Trinity.fasta_stripped ./
		python ${params.scriptDir}/rename_trinity.py Trinity.fasta_stripped ${params.CODE}
		mv ${params.CODE}_Trinity.fasta_stripped "${params.trinityOutDir}/${params.CODE}_Trinity.fasta"
		rm Trinity.fasta_stripped
		"""
}


// Finds open reading frames within the transcriptome from Trinity and translates them
process TRANSDECODER {
	publishDir "${params.outDir}", mode: 'copy'
	input:
		path transcriptome
	output:
		path "${params.transdecoderOutDir}/${params.CODE}_Trinity.fasta.transdecoder.cds"
		path "${params.transdecoderOutDir}/${params.CODE}_Trinity.fasta.transdecoder.pep"
	script:
		"""
		# Runs TransDecoder
		TransDecoder.LongOrfs -t ${transcriptome} -O ${params.transdecoderTemp}
		TransDecoder.Predict -t ${transcriptome} -O ${params.transdecoderTemp}
		
		# Makes backups of each original file and reformats the main outputs
		mkdir -p ${params.transdecoderOutDir}
		for i in ${params.transdecoderTemp}/${params.CODE}_Trinity.fasta.transdecoder.{cds,pep}
		do
			cp \$i "\${i%/*}/originalFile_\${i##*/}"
			python ${params.scriptDir}/strip_returns_fasta.py \$i
			mv \$i\\_stripped \$i
			mv \$i ${params.transdecoderOutDir}
		done
		"""
}


// Clusters the transcriptome by a specified percent sequence identity. Default value is 99%
// Note that this uses local alignments for anything covering at least 'params.CDHrepAlnLen'% of the larger sequence. Check CD-HIT documentation to ensure this is appropriate for your use.
process CDHIT {
	publishDir "${params.outDir}", mode: 'copy'
	input:
		path transcriptome
	output:
		path "${params.CDHoutput}"
	script:
		"""
		mkdir -p ${params.CDHoutDir}
		echo ${params.CDHoutput}
		cd-hit-est -i ${transcriptome} -o ${params.CDHoutput} -d ${params.CDHdescLen} -c ${params.CDHpercentID} -n ${params.CDHwordSize} -M ${params.CDHmemory} -T ${params.CDHthreads} -G ${params.CDHalnType} -aL ${params.CDHrepAlnLen}
		"""
}


// Filters the translated peptide file by associated sequences which were not filtered out when clustering the transcriptome in the CDHIT process
// This is done to ensure that both the peptide and nucleotide files do not have a mismatch in which sequences where filtered out by percent identity.
process MATCHFASTAS {
	publishDir "${params.outDir}", mode: 'copy'
	input:
		path clusteredTranscriptome
		path translatedPEP
	output:
		path "${params.matchFastasOutput}"
	script:
		"""
		mkdir -p ${params.CDHoutDir}
		python3 ${params.scriptDir}/clustering-fastaMatch.py ${clusteredTranscriptome} ${translatedPEP} ${params.matchFastasOutput}
		"""
}


// Makes an index file out of the clustered transcriptome from CDHIT. This will be used to map reads from individual RNAseq samples.
process SALMON_INDEX {
	publishDir "${params.outDir}", mode: 'copy'
	input:
		path clusteredTranscriptome
	output:
		path "${params.salmonIndexOutput}"
	script:
		"""
		mkdir -p ${params.salmonOutDir}
		salmon index --threads ${params.threads} -t ${clusteredTranscriptome} -i ${params.salmonIndexOutput}
		"""
}

// Maps individual reads from specified fastq file pairs to the index prepared with SALMON_INDEX
process SALMON_MAP {
	publishDir "${params.outDir}", mode: 'copy'
	input:
		path salmonIndex
		tuple val(sampleName), path(individualReads)
	output:
		path "${params.salmonQuantsDir}/quant.${sampleName}"
	script:
		"""
		salmon quant -i ${salmonIndex} -l ${params.salmonLibType} -1 ${individualReads[0]} -2 ${individualReads[1]} -p ${params.salmonThreads} --validateMappings -o "${params.salmonQuantsDir}/${sampleName}.quant"
		mv "${params.salmonQuantsDir}/${sampleName}.quant/quant.sf" "${params.salmonQuantsDir}/quant.${sampleName}"
		"""
}

// Maps individual reads from a specified unpaired fastq file to the index prepared with SALMON_INDEX
// This is an alternative for the SALMON_MAP process above when working with unparied data
process SALMON_MAP_SINGLE {
	publishDir "${params.outDir}", mode: 'copy'
	input:
		path salmonIndex
		path individualReads
	output:
		path "${params.salmonQuantsDir}/quant.${sampleName}"
	script:
		"""
		salmon quant -i ${salmonIndex} -r ${individualReads} -2 ${individualReads[1]} -p ${params.salmonThreads} --validateMappings -o "${params.salmonQuantsDir}/${sampleName}.quant"
		mv "${params.salmonQuantsDir}/${sampleName}.quant/quant.sf" "${params.salmonQuantsDir}/quant.${sampleName}"
		"""
}


// Combines all individual outputs from SALMON_MAP or SALMON_MAP_SINGLE and merges them into a single matrix
// Outputs separate matrices for raw counts and TPM
process SALMON_MERGE {
	publishDir "${params.outDir}", mode: 'copy'
	input:
		path expressionMatrices
	output:
		path "${params.salmonOutDir}/${params.salmonOutTPM}"
		path "${params.salmonOutDir}/${params.salmonOutCounts}"
	script:
		"""
		python ${params.scriptDir}/salmonOut_merge.py ${params.CODE} ${expressionMatrices}
        mkdir ${params.salmonOutDir}
        mv ${params.salmonOutTPM} ${params.salmonOutDir}
        mv ${params.salmonOutCounts} ${params.salmonOutDir}
		"""
}


// Specify workflow here
// Ensure that you run the "_SINGLE" processes for TRINITY and SALMON_MAP when working with unpaired data
workflow {
	
	// --- For an assembly with paired-end reads, use the following for Trinity
	assemblyReads_ch = Channel.fromFilePairs(params.assemblyReads)
	trinity_ch = TRINITY(assemblyReads_ch)
	// --- For an assembly with single-end reads, use the following for Trinity
	//assemblyReads_ch = channel.fromPath(params.assemblyReads)
	//trinity_ch = TRINITY_SINGLE(assemblyReads_ch)
	
	
	transdecoder_ch = TRANSDECODER(trinity_ch)
	
	
	cdhit_ch = CDHIT(trinity_ch)
	
	
	matchFastas_ch = MATCHFASTAS(cdhit_ch, transdecoder_ch[1])
	
	
	salmonIndex_ch = SALMON_INDEX(cdhit_ch)
	
	
	// --- For an assembly with paired-end reads, use the following as input for mapping reads with Salmon
	mappingReads_ch = Channel.fromFilePairs("${params.mappingReads}/*_{1,2}.{fastq,fq}")
	mappedReads_ch = SALMON_MAP(salmonIndex_ch, mappingReads_ch).collect()
	// --- For an assembly with single-end reads, use the following as input for mapping reads with Salmon
	//mappingReads_ch = channel.fromPath("${params.mappingReads}/*.{fastq,fq}")
	//mappedReads_ch = SALMON_MAP_SINGLE(salmonIndex_ch, mappingReads_ch).collect()
	
	
	expMatrices_ch = SALMON_MERGE(mappedReads_ch)
	
}
