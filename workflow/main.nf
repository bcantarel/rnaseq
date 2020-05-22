#!/usr/bin/env nextflow

// Default parameter values to run tests
params.input = "$baseDir"
params.output= "$baseDir/output"
params.design="$params.input/design.txt"
params.genome="/project/shared/bicf_workflow_ref/human/GRCh38/"
params.markdups="picard"
params.stranded="0"
params.pairs="pe"
params.geneset = 'h.all.v5.1.symbols.gmt'
params.align = 'hisat'
params.fusion = 'skip'
params.dea = 'detect'

design_file = file(params.design)
fqs = Channel.fromPath("$params.input/*")
gtf_file = file("$params.genome/gencode.gtf")
genenames = file("$params.genome/genenames.txt")
geneset = file("$params.genome/../gsea_gmt/$params.geneset")
dbsnp="$params.genome/dbSnp.vcf.gz"
indel="$params.genome/GoldIndels.vcf.gz"
knownindel=file(indel)
dbsnp=file(dbsnp)

// params genome is the directory
// base name for the index is always genome
index_path = file(params.genome)

process checkdesignfile {
	queue 'super'
	module 'parallel/20150122:pigz/2.4'
	publishDir "$params.output", mode: 'copy'

	input:
        file design_file name 'design.ori.txt'
	file ("*") from fqs.collect()

	output:
	file("design.valid.txt") into newdesign
	file("*.fastq.gz") into fastqs mode flatten
	stdout spltnames

	script:
	"""
	bash $baseDir/scripts/check_inputfiles.sh 
	perl -p -e 's/\\r\\n*/\\n/g' design.ori.txt > design.fix.txt
	perl $baseDir/scripts/check_designfile.pl ${params.pairs} design.fix.txt
	"""
}

def fileMap = [:]

fastqs
    .subscribe { 
	def fileName = it.getFileName()
	fileMap."$fileName" = it
    }

if (params.pairs == 'pe') {
  spltnames
  .splitCsv()
  .filter { fileMap.get(it[1]) != null & fileMap.get(it[2]) != null }
  .map { it -> tuple(it[0], fileMap.get(it[1]), fileMap.get(it[2])) }
  .set { read }
} else {
  spltnames
  .splitCsv()
  .filter { fileMap.get(it[1]) != null }
  .map { it -> tuple(it[0], fileMap.get(it[1]),'') }
  .set { read }
}
if( ! read ) { error "Didn't match any input files with entries in the design file" }

// Trim raw reads using trimgalore
process trim {
  errorStrategy 'ignore'
  input:
  set pair_id, file(read1), file(read2) from read
  output:
  set pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into trimread
  set pair_id, file("${pair_id}.trim.R1.fastq.gz"),file("${pair_id}.trim.R2.fastq.gz") into fusionfq
  script:
  """
  bash $baseDir/process_scripts/preproc_fastq/trimgalore.sh -p ${pair_id} -a ${read1} -b ${read2}
  """
}

// Align trimmed reads to genome indes with hisat2
// Sort and index with samtools
// QC aligned reads with fastqc
// Alignment stats with samtools

process starfusion {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set pair_id, file(fq1), file(fq2) from fusionfq
  output:
  file("${pair_id}.starfusion.txt") into fusionout
  when:
  params.fusion == 'detect' && params.pairs == 'pe'
  script:
  """
  bash $baseDir/process_scripts/alignment/starfusion.sh -p ${pair_id} -r ${index_path} -a ${fq1} -b ${fq2} -m trinity -f
  """
}

process align {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set pair_id, file(fq1), file(fq2) from trimread
  output:
  set pair_id, file("${pair_id}.bam") into aligned
  set pair_id, file("${pair_id}.bam") into aligned2
  file("${pair_id}.alignerout.txt") into hsatout
  script:
  """
  bash $baseDir/process_scripts/alignment/rnaseqalign.sh -a $params.align -p ${pair_id} -r ${index_path} -x ${fq1} -y ${fq2}
  """
}

process alignqc {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set pair_id, file(bam) from aligned2
  output:
  file("${pair_id}.flagstat.txt") into alignstats
  set file("${pair_id}_fastqc.zip"),file("${pair_id}_fastqc.html") into fastqc
  script:
  """
  bash $baseDir/process_scripts/alignment/bamqc.sh -p ${pair_id} -b ${bam} -y rna
  """
}

// Summarize all flagstat output

process parse_alignstat {
  publishDir "$params.output", mode: 'copy'
  input:
  file(txt) from alignstats.toList()
  file(txt) from  hsatout.toList()
  output:
  file('alignment.summary.txt')
  script:
  """
  perl $baseDir/scripts/parse_flagstat.pl *.flagstat.txt
  """
}

// Identify duplicate reads with Picard

process markdups {
  publishDir "$params.output", mode: 'copy'
  input:
  set pair_id, file(sbam) from aligned
  output:
  set pair_id, file("${pair_id}.dedup.bam") into deduped1
  set pair_id, file("${pair_id}.dedup.bam") into deduped2
  script:
  """
  bash $baseDir/process_scripts/alignment/markdups.sh -a $params.markdups -b $sbam -p $pair_id
  """
}

// Read summarization with subread
// Assemble transcripts with stringtie

process geneabund {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set pair_id, file(sbam) from deduped1
  output:
  file("${pair_id}.cts") into counts
  file("${pair_id}.cts.summary") into ctsum
  file("${pair_id}_stringtie") into strcts
  file("${pair_id}.fpkm.txt") into fpkm
  script:
  """
  bash $baseDir/process_scripts/genect_rnaseq/geneabundance.sh -s $params.stranded -g ${gtf_file} -p ${pair_id} -b ${sbam}
  """
}

process statanal {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  file count_file from counts.toList()
  file count_sum from ctsum.toList()
  file newdesign name 'design.txt'
  file genenames
  file geneset name 'geneset.gmt'
  file fpkm_file from fpkm.toList()
  file stringtie_dir from strcts.toList()
  output:
  file "*.txt" into txtfiles
  file "*.png" into psfiles
  file("*.rda") into rdafiles
  file("geneset.shiny.gmt") into gmtfile
  script:
  """
  bash $baseDir/process_scripts/genect_rnaseq/statanal.sh -d $params.dea
  """
}
