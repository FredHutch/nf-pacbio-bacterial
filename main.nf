#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/nf-pacbio-bacterial <ARGUMENTS>
    
    Required Arguments:
      --input_folder        Folder containing all PacBio data in *.bam files (including subdirectories)
      --output_folder       Folder to place analysis outputs

    Optional Arguments (passed directly to flye):
      --read_type           Type of PacBio or MinION reads passed in to Flye, either raw, corrected, or HiFi (PacBio-only)
                            Default: pacbio-raw
                            Options: pacbio-raw, pacbio-corr, pacbio-hifi, nano-raw, nano-corr, subassemblies
      --iterations          Number of polishing iterations
                            Default: 1
    
    For more details on Flye, see https://github.com/fenderglass/Flye/blob/flye/docs/USAGE.md

    """.stripIndent()
}

// Show help message if the user specifies the --help flag at runtime
params.help = false
if (params.help || params.input_folder == null || params.output_folder == null){
    // Invoke the function above which prints the help message
    helpMessage()
    // Exit out and do not run anything else
    exit 1
}

// Default options listed here
params.read_type = "pacbio-raw"
params.iterations = 1

/////////////////////
// DEFINE FUNCTIONS /
/////////////////////

// Extract reads from BAM to FASTQ format
process extractBAM {

  // Docker container to use
  container "quay.io/biocontainers/bam2fastx:1.3.1--he1c1bb9_0"
  label "io_limited"
  errorStrategy 'finish'

  input:
    tuple val(name), file(bam)

  // The block below points to the files inside the process working directory which will be retained as outputs (and published to the destination above)
  output:
  tuple val(name), file("${name}/*")

"""
#!/bin/bash

set -Eeuo pipefail

bam2fastq -o ${name} ${bam}

"""

}

// Run Flye
process flye {

  // Docker container to use
  container "quay.io/biocontainers/flye:2.8--py37h8270d21_0"
  label "mem_medium"
  errorStrategy 'finish'

  // The `publishDir` tag points to a folder in the host system where all of the output files from this folder will be placed
  publishDir "${params.output_folder}" 
  
  input:
    tuple val(name), file(reads)

  // The block below points to the files inside the process working directory which will be retained as outputs (and published to the destination above)
  output:
  file "${name}/*"

"""
#!/bin/bash

set -e

df -h
echo ""
ls -lahtr
echo ""
echo "STARTING FLYE"
echo ""

flye \
    --${params.read_type} ${reads} \
    --out-dir ${name} \
    --threads ${task.cpus} \
    --iterations ${params.iterations} \
    --plasmids

"""

}

// Start the workflow
workflow {

    // Get the input files ending with BAM
    bam_ch = Channel.fromPath(
        "${parmas.input_folder}**.bam"
    ).map {
        it -> (it.name.replaceAll(/.bam/, ''), it)
    }

    // Extract the BAM files to FASTQ
    extractBAM(
        bam_ch
    )

}
