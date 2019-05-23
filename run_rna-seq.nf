#!/usr/bin/env nextflow

params.kallisto_index='/lustre/scratch/projects/rnaseq_nod/indices/kallisto/araport11_all_H2b-mCherry'
params.annotation ='/lustre/scratch/projects/rnaseq_nod/indices/star/files/Araport11_TES_H2b-mCherry.gtf'
params.star_index='/lustre/scratch/projects/rnaseq_nod/indices/star/araport_tetools/'

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
//    .fromFilePairs(params.input)
//	// { file -> tuple(file.name.replaceAll(/.{1,2}.fq$/,''), file) }
//    //.fromFilePairs(params.input +'*.{1,2}.fq') { file -> tuple(file.name.replaceAll(/.{1,2}.fq$/,''), file) }
//    .ifEmpty { error "Cannot find any fq files: $params.input." }
//    .set { fastq_files }

//Channel
//    .fromFilePairs(params.input)
//	// { file -> tuple(file.name.replaceAll(/.{1,2}.fq$/,''), file) }


Channel
   .fromPath(params.input)
   .ifEmpty { error "Cannot find any bam files: $params.input." }
   .map { file -> tuple(file.name.replaceAll(/.bam$/,''), file) }
   .set { bam_files }

Channel
    .fromPath(params.kallisto_index)
    .ifEmpty { error "Cannot find kallisto index: $params.kallisto_index." }
    .set { kallisto_idx }

Channel
    .fromPath(params.star_index)
    .ifEmpty { error "Cannot find kallisto index: $params.star_index." }
    .set { star_idx }

Channel
    .fromPath(params.annotation)
    .ifEmpty { error "Cannot find annotation file: $params.annotation." }
    .set { align_annotation }
    //.into { index_annotation, align_annotation }


process sort_bam{
        tag "sorting bam: $name"

        input:
        set name, file(bam) from bam_files.dump(tag: 'input')

        output:
        set val(name), file("${name}.sorted.bam") into sorted_bam
     
        """
        samtools sort -n -@ ${task.cpus} -o "${name}.sorted.bam" ${name}.bam
        """
}

process bam_to_fastq{

    tag "converting to fastq: $name"

    input:
    set name, file(bam) from sorted_bam.dump(tag: 'bam_to_fq')
    
    output:
    set name, file("*.fq") into fastq_files
 
    """
    samtools fastq -1 ${name}.1.fq -2 ${name}.2.fq\
      -0 /dev/null -s /dev/null -N -F 0x900 $bam
    """
}

process trim_adapters{
    
    tag "trimming: $name"

    publishDir "$params.output/$name/trimmed", mode: 'copy', pattern: "*_trimming_report.txt"

    input:
    set name, file(fastq) from fastq_files.dump(tag: 'trim_input')

    output:
    file "*"
    set name, file("*_val_{1,2}.fq") into fastq_kallisto, fastq_star

    """
    trim_galore --dont_gzip \
      --length 18 \
      --stringency 4 \
      --paired ${name}.1.fq ${name}.2.fq
    """

}

process quantify_kallisto{

    tag "quantifiying: $name"
    publishDir "$params.output/$name/kallisto", mode: 'copy'
    //publishDir "$params.output/$name/kallisto", mode: 'copy', pattern: "${name}/*", saveAs: { filename -> "${name}_$filename" }

    input:
    set name, file(fastq) from fastq_kallisto.dump(tag: 'kallisto_input')
    file index from kallisto_idx.collect()

    output:
    //file "${name}/*"
    file "*.{tsv,h5,json}"
    
    """
    kallisto quant -i $index -o ${name} ${name}.1_val_1.fq ${name}.2_val_2.fq
    mv ${name}/abundance.h5 ${name}_abundance.h5
    mv ${name}/abundance.tsv ${name}_abundance.tsv
    mv ${name}/run_info.json ${name}_run_info.json
    rmdir ${name}
    """
}

process align_star{
        
    tag "quantifiying: $name"
    
    publishDir "$params.output/$name/star", mode: 'copy'

    input:
    set name, file(fastq) from fastq_star.dump(tag: 'star_input')
    file index from star_idx.collect()
    file anno from align_annotation.collect()

    output: 
    file "*"

    """
    STAR --runMode alignReads \
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
    --readFilesIn ${name}.1_val_1.fq  ${name}.2_val_2.fq
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
