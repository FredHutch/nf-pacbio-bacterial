#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/nf-pacbio-bacterial <ARGUMENTS>
    
    Required Arguments:
      --input_folder        Folder containing PacBio output data
      --output_folder       Folder to place analysis outputs

    Optional QC Arguments:
      --min_length          Minimum read length filter (default: 1000)
      
    For more details on SequelTools, see https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03751-8

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
params.help = false
params.input_folder = null
params.prefix = null
params.output_folder = null
if (params.help || params.input_folder == null || params.output_folder == null || params.prefix == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 1
}

// Default options listed here
params.read_type = "pacbio-corr"
params.min_qscore = 25
params.min_qscore_pct = 90
params.min_length = 1000
params.iterations = 1
params.subset_n_reads = false
params.fastqc_max_reads = 10000

/////////////////////
// DEFINE FUNCTIONS /
/////////////////////

// Run QC with SequelTools
process sequeltools_QC {

  // Docker container to use
  container "quay.io/fhcrc-microbiome/sequeltools:latest"
  label "mem_medium"
  errorStrategy 'finish'

  publishDir "${params.output_folder}/${prefix}/qc/" 
  
  input:
    tuple val(prefix), file(subreads_bam), file(subreads_pbi), file(scraps_bam), file(scraps_pbi)

  output:
    file "SequelToolsResults/*"

"""
#!/bin/bash

set -Eeuo pipefail

echo "${subreads_bam.name}" > subFiles.txt
echo "${scraps_bam.name}" > scrFiles.txt

# Run QC
bash /usr/local/sequeltools/SequelTools/Scripts/SequelTools.sh \
    -v \
    -n ${task.cpus} \
    -t Q \
    -u subFiles.txt \
    -c scrFiles.txt

"""

}

// Run subsampling with SequelTools
process sequeltools_subsampling {

  // Docker container to use
  container "quay.io/fhcrc-microbiome/sequeltools:latest"
  label "mem_medium"
  errorStrategy 'finish'

  publishDir "${params.output_folder}/${prefix}/subsampling/" 
  
  input:
    tuple val(prefix), file(subreads_bam), file(subreads_pbi), file(scraps_bam), file(scraps_pbi)

  output:
    file "SequelToolsResults/*"

"""
#!/bin/bash

set -Eeuo pipefail

echo "${subreads_bam.name}" > subFiles.txt
echo "${scraps_bam.name}" > scrFiles.txt

# Run Subsampling
bash /usr/local/sequeltools/SequelTools/Scripts/SequelTools.sh \
    -v \
    -n ${task.cpus} \
    -t S \
    -T l \
    -u subFiles.txt \
    -c scrFiles.txt

"""

}

// Run read filtering with SequelTools
process sequeltools_filtering {

  // Docker container to use
  container "quay.io/fhcrc-microbiome/sequeltools:latest"
  label "mem_medium"
  errorStrategy 'finish'

  publishDir "${params.output_folder}/${prefix}/filtering/" 
  
  input:
    tuple val(prefix), file(subreads_bam), file(subreads_pbi), file(scraps_bam), file(scraps_pbi)

  output:
    file "SequelToolsResults/*"

"""
#!/bin/bash

set -Eeuo pipefail

echo "${subreads_bam.name}" > subFiles.txt
echo "${scraps_bam.name}" > scrFiles.txt

# Run Read Filtering
bash /usr/local/sequeltools/SequelTools/Scripts/SequelTools.sh \
    -v \
    -n ${task.cpus} \
    -t F \
    -C -P -N -Z ${params.min_length} \
    -u subFiles.txt \
    -c scrFiles.txt

"""

}

// Start the workflow
workflow {

    // Get the input files ending with {subreads | scraps}.bam(.pbi)
    subreads_bam_ch = Channel.fromPath(
        "${params.input_folder}**.subreads.bam"
    ).map { it -> [it.name.replaceAll(/.subreads.bam/, ''), it]}
    
    subreads_pbi_ch = Channel.fromPath(
        "${params.input_folder}**.subreads.bam.pbi"
    ).map { it -> [it.name.replaceAll(/.subreads.bam.pbi/, ''), it]}

    scraps_bam_ch = Channel.fromPath(
        "${params.input_folder}**.scraps.bam"
    ).map { it -> [it.name.replaceAll(/.scraps.bam/, ''), it]}
    
    scraps_pbi_ch = Channel.fromPath(
        "${params.input_folder}**.scraps.bam.pbi"
    ).map { it -> [it.name.replaceAll(/.scraps.bam.pbi/, ''), it]}

    combined_ch = subreads_bam_ch.join(
        subreads_pbi_ch
    ).join(
        scraps_bam_ch
    ).join(
        scraps_pbi_ch
    )

    // Run SequelTools
    sequeltools_QC(
        combined_ch
    )
    sequeltools_filtering(
        combined_ch
    )
    sequeltools_subsampling(
        combined_ch
    )

}
