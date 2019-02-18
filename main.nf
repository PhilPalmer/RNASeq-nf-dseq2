/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'RNASEQ-NF'.
 *
 *   RNASEQ-NF is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   RNASEQ-NF is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with RNASEQ-NF.  If not, see <http://www.gnu.org/licenses/>.
 */
 
 
/* 
 * Proof of concept of a RNAseq pipeline implemented with Nextflow
 * 
 * Authors:
 * - Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * - Emilio Palumbo <emiliopalumbo@gmail.com> 
 * - Evan Floden <evanfloden@gmail.com> 
 */ 

 
/*
 * Default pipeline parameters. They can be overriden on the command line eg. 
 * given `params.foo` specify on the run command line `--foo some_value`.  
 */
 
//params.reads = "$baseDir/data/ggal/*_{1,2}.fq"
params.reads = false
readsChannel = "${params.reads}/*_{1,2}.fq"

params.transcriptome = false //"$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.accession     = false
params.outdir = "results"
params.multiqc = "$baseDir/multiqc"

log.info """\
 R N A S E Q - N F   P I P E L I N E    
 ===================================
 transcriptome: ${params.transcriptome}
 reads        : ${params.reads}
 outdir       : ${params.outdir}
 readChannel  : ${readsChannel}
 """


transcriptome_file = file(params.transcriptome)
multiqc_file = file(params.multiqc)
//projectSRId = params.project
accessionID = params.accession
 

Channel
    .fromPath( params.reads )
    .ifEmpty { error "Cannot find any reads matching" }
    .set { sraIDs } 

// Channel.fromPath(params.reads)
//     .ifEmpty { exit 1, "Text file containing SRA id's not found: ${params.reads}" }
//     .into { sraIDs; rna_sraIDs }

int threads = Runtime.getRuntime().availableProcessors()

// process getSRAIDs {
//     container 'lifebitai/kallisto-sra'
	
// 	cpus 1

// 	input:
// 	val projectID from projectSRId
	
// 	output:
// 	file 'sra.txt' into sraIDs
	
// 	script:
// 	"""
// 	esearch -db sra -query $projectID  | efetch --format runinfo | grep SRR | cut -d ',' -f 1 > sra.txt
// 	"""
// }
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
    .map { file -> tuple(file.baseName, file) }
    .into { reads_fastqc; reads_quant }

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
    publishDir params.outdir, mode: 'copy'
     
    input:
    file index from index_ch
    set val(name), file(reads) from reads_quant
 
    output:
    file(name) into quant
    val(name) into reads_deseq
 
    script:
    """
    salmon quant --threads $task.cpus --libType=U -i index -r $reads -o $name
    """
}

quant.into { quant_deseq; quant_multiqc }

process deseq2 {
	publishDir params.outdir, mode: 'copy'
    container 'lifebitai/rnaseq-nf-salmon-dseq2'

    input:
    file(quant) from quant_deseq.collect()

    output:
    //file('counts.RData') into countData

    script:
    """
    #!/usr/bin/env Rscript

    library("tximport")
    library("readr")
    library("tximportData")

    library("DESeq2")


    #dir <- system.file("extdata", package="tximportData")
    samples <- read.table(file.path(dir,"samples.txt"), header=TRUE)
    samples\$condition <- factor(rep(c("A","B"),each=3))
    rownames(samples) <- samples\$run
    samples[,c("pop","center","run","condition")]

    files <- file.path(dir,"salmon", samples\$run, "quant.sf.gz")
    names(files) <- samples\$run

    tx2gene <- read_csv(file.path(".", ".csv"))

    txi <- tximport(files, type="salmon", tx2gene=tx2gene)


    
    #ddsTxi <- DESeqDataSetFromTximport(txi,
    #                            colData = samples,
    #                            design = ~ condition)
    """
}


// process multiqc {
//     publishDir params.outdir, mode:'copy'
       
//     input:
//     file('*') from quant_multiqc.mix(fastqc_ch).collect()
//     file(config) from multiqc_file
    
//     output:
//     file('multiqc_report.html')  
     
//     script:
//     """
//     cp $config/* .
//     echo "custom_logo: \$PWD/logo.png" >> multiqc_config.yaml
//     multiqc . 
//     """
// }
 
workflow.onComplete { 
	println ( workflow.success ? "\nDone! Open the following report in your browser --> $params.outdir/multiqc_report.html\n" : "Oops .. something went wrong" )
}
