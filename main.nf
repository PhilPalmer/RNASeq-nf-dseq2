params.reads = false
params.transcriptome = false 
params.accession     = false
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"

log.info """\
 R N A S E Q - N F - D S E Q 2  P I P E L I N E    
 ===================================
 transcriptome: ${params.transcriptome}
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 """


transcriptome_file = file(params.transcriptome)
multiqc_file = file(params.multiqc)
accessionID = params.accession
 

Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching" }
    .into { reads_deseq; sra_file; sra_desseq2 } 

int threads = Runtime.getRuntime().availableProcessors()

process preprocess_sra {
    tag "$samples"

	input:
	file samples from sra_file

	output:
	file 'sras.txt' into sraIDs

	script:
	"""
	awk -F, '{print \$1}' $samples | tail -n +2 > sras.txt
	"""	
}

sraIDs.splitText().map { it -> it.trim() }.set { singleSRAId }

process fastqDump {
    tag "$id"
    container 'lifebitai/kallisto-sra'

	//publishDir params.resultdir, mode: 'copy'

	cpus threads

	input:
	val id from singleSRAId

	output:
	file '*.fastq.gz' into reads

	script:
	"""
	parallel-fastq-dump --sra-id $id --threads ${task.cpus} --gzip
	"""	
}

reads
    .map { file -> tuple(file.simpleName, file) }
    .into { reads_fastqc; reads_quant; sra }

process fastqc {
    tag "FASTQC on $sample_id"

    input:
    set val(sample_id), file(reads) from reads_fastqc

    output:
    file("fastqc_${sample_id}_logs") into fastqc_ch


    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}  

process index {
    tag "$transcriptome_file.simpleName"
    
    input:
    file transcriptome from transcriptome_file
     
    output:
    file 'index' into index_ch

    script:       
    """
    salmon index --threads $task.cpus -t $transcriptome -i index
    """
}
 
process quant {
    tag "$name"
    publishDir "${params.outdir}/salmon", mode: 'copy'
    container 'lifebitai/rnaseq-nf-dseq2'
     
    input:
    file index from index_ch
    set val(name), file(reads) from reads_quant
 
    output:
    file(name) into quant
    //val(name) into reads_deseq
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i index -r $reads -o $name
    """
}

quant.into { quant_deseq; quant_multiqc }

process deseq2 {
    tag "deseq2.RData"
	publishDir params.outdir, mode: 'copy'
    container 'lifebitai/rnaseq-nf-dseq2'

    input:
    file(quant) from quant_deseq.collect()
    val(reads) from reads_deseq
    file samples from sra_desseq2

    output:
    file('deseq.RData')

    script:
    """
    #!/usr/bin/env Rscript

    library("tximport")
    library("readr")
    library("tximportData")
    library("DESeq2")

    dir <- system.file("extdata", package="tximportData")
    samples <- read.table("$samples", header=TRUE, sep=",")

    tx2gene <- read_csv(file.path(dir, "tx2gene.gencode.v27.csv"))

    files <- file.path(".",samples\$samples, "quant.sf")
    names(files) <- samples\$samples

    txi <- tximport(files, type="salmon", tx2gene=tx2gene, ignoreAfterBar=TRUE)

    condition <- samples\$condition

    ddsTxi <- DESeqDataSetFromTximport(txi,
                                    colData = samples,
                                    design = ~ condition)

    save(ddsTxi, file="deseq.RData")
    """
}


process multiqc {
    publishDir "${params.outdir}/MultiQC", mode:'copy'
       
    input:
    file('*') from quant_multiqc.mix(fastqc_ch).collect()
    
    output:
    file('*')  
     
    script:
    """
    multiqc . 
    """
}
 
workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
