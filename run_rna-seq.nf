#!/usr/bin/env nextflow

params.kallisto_index='/lustre/scratch/projects/rnaseq_nod/indices/kallisto/Araport11_H2b-mCherry_cdna'
params.star_index='/lustre/scratch/projects/rnaseq_nod/indices/star/araport_tetools/'
params.annotation ='/lustre/scratch/projects/rnaseq_nod/indices/star/araport_tetools/'

log.info "RNA-seq NF  ~  version 3.1"
log.info "====================================="
log.info "input paths: ${params.input}"
log.info "output paths: ${params.output}"
log.info "kallisto index: ${params.kallisto_index}"
log.info "star index: ${params.star_index}"
log.info "annotation file: ${params.annotation}"
log.info "\n"

// input channels

//Channel
//    .fromPath(params.input)
//    .splitCsv(header:true)
//    .filter(row.sequencing.type == 'PE')
//    .filter(row.file.type == 'bam')
//    .map{ row-> tuple(row.sampleId, file(row.read1), file(row.read2)) }
//    .set { samples_ch }
// Channel
//     .fromPath(params.genome)
//     .ifEmpty { error "Cannot find genome fasta file: $params.genome." }
//     .set { genome_fasta }



Channel
    .fromFilePairs('reads/*.{1,2}.fq')
    .ifEmpty { error "Cannot find any fq files: $params.input." }
    .map{ file -> tuple(file.name.replaceAll(/.{1,2}.fq$/,''), file) }
    .set { fastq_files }

Channel
    .fromPath(params.kallisto_index)
    .ifEmpty { error "Cannot find kallisto index: $params.kallisto_index." }
    .set { kallisto_idx }


Channel
    .fromPath(params.annotation)
    .ifEmpty { error "Cannot find annotation file: $params.annotation." }
    .into { align_annotation }
    //.into { index_annotation, align_annotation }


process trim_adapters{
    
    tag "quantifiying: $name"

    publishDir "$params.output/$name/fastq/trimmed", mode: 'copy'

    input:
    set name, file("*.fq") from fastq_files.dump(tag: 'tim_input')

    output:
    file "*"
    set name, file("*_val_{1,2}.fq") into fastq_kallisto, fastq_star

    """
    trim_galore --dont_gzip \
      --length 18 \
      --stringency 4 \
      --cores ${task.cpus} \
      --paired ${name}.1.fq ${name}.2.fq
    """

}

process quantify_kallisto{

    tag "quantifiying: $name"
    
    publishDir "$params.output/$name/kallisto", mode: 'copy', pattern: "$name/*"

    input:
    set name, file(fastq) from fastq_quant.dump(tag: 'kallisto_input')
    file index from kallisto_idx

    output:
    file "*"
    
    """
    kallisto quant -i $index -o $name ${name}_val_1.fq ${name}_val_2.fq
    """
}

process align_star{
        
    tag "quantifiying: $name"
    
    publishDir "$params.output/$name/star", mode: 'copy'

    input:
    set name, file(fastq) from fastq_star.dump(tag: 'star_input')
    file index from star_idx
    file anno from annotation

    output: 
    file "*"

    """
    STAR
    --runMode alignReads \
    --alignIntronMax 5000 \
    --alignMatesGapMax 5500 \
    --alignSJDBoverhangMin 1 \
    --outFilterMismatchNmax 5 \
    --outFilterMismatchNoverLmax .05 \
    --outFilterType BySJout \
    --genomeDir ${index} \
    --limitBAMsortRAM "${task.memory}"000000000 \
    --outFilterIntronMotifs RemoveNoncanonicalUnannotated \
    --outFilterMultimapNmax 100 \
    --winAnchorMultimapNmax 100 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMprimaryFlag AllBestScore \
    --outTmpDir ${name} \
    --runThreadN ${task.cpus} \
    --sjdbGTFfile ${anno} \
    --outSJfilterOverhangMin -1 16 -1 -1 \
    --outSJfilterCountTotalMin -1 2 -1 -1 \
    --outSJfilterCountUniqueMin -1 2 -1 -1 \
    --outFileNamePrefix ${name}. \
    --readFilesIn ${name}_val_1.fq ${name}_val_2.fq
    """
}

// process build_star_index{

//     tag "building star index"
//     storeDir '/lustre/scratch/projects/rnaseq_nod/indices/star/'

//     input
//     file fasta from genome_fasta
//     file annotation from index_annotation

//     """
//     star --runMode genomeGenerate
//         --runThreadN ${task.cpus} \
//         --genomeDir  \
//         --genomeFastaFiles $fasta
//         --sjdbGTFfile $annotation
//         --sjdbOverhang ReadLength-1
//     """

// }

// process sort_bam{
//         tag "sorting bam: $name"

//         input:
//         set name, file(bam) from bam_files(tag: 'input')

//         output:
//         set val(name), file("${name}.sorted.bam") into sorted_bam
     
//         """
//         samtools sort -n -@ ${task.cpus} -o "${name}.sorted.bam" ${name}.bam
//         """
// }

// process bam_to_fastq{

//     tag "converting to fastq: $name"

//     publishDir "$params.output/$name/fastq", mode: 'copy'

//     input:
//     set name, file(bam) from sorted_bam.dump(tag: 'bam_to_fq')
    
//     output:

//     set name, file("*.fq") into fastq
 
//     """
//     samtools fastq -1 ${name}.1.fq -2 ${name}.2.fq\
//       -0 /dev/null -s /dev/null -N -F 0x900 -@ ${task.cpus} $bam
//     """
// }


// process build_kallisto_index{   

//     storeDir "/lustre/scratch/projects/rnaseq_nod/indices/kallisto/$name"

//     input:
//     set name, file(fasta) from transcripts_fasta

//     output:
//     file "${name}.kallisto"

//     """
//     kallisto index -i ${name}.kallisto $fasta
//     """
// }
