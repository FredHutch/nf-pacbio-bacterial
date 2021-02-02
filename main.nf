#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/nf-pacbio-bacterial <ARGUMENTS>
    
    Required Arguments:
      --input_folder        Folder containing all PacBio data in BAM files (including subdirectories)
      --output_folder       Folder to place analysis outputs

    Input Files:
      --suffix              Process all files ending with this string (default: .subreads.bam)
                            Note that files must have a paired .bam.pbi file in the same folder.

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
params.suffix = ".subreads.bam"
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
    tuple file(bam), file(pbi)

  output:
    file "output/*"

"""
#!/bin/bash

set -Eeuo pipefail

mkdir output

bam2fastq -o output ${bam}

"""

}

// Run Flye
process flye {

  // Docker container to use
  container "quay.io/biocontainers/flye:2.8--py37h8270d21_0"
  label "mem_medium"
  errorStrategy 'finish'

//   publishDir "${params.output_folder}" 
  
  input:
    tuple val(name), file(reads)

  output:
  file "${name}/*"

"""
#!/bin/bash

set -e

flye \
    --${params.read_type} ${reads} \
    --out-dir ${name} \
    --threads ${task.cpus} \
    --iterations ${params.iterations} \
    --plasmids

"""

}

// Run FastQC
process fastQC {

  // Docker container to use
  container "quay.io/biocontainers/fastqc:0.11.9--0"
  label "io_limited"
  errorStrategy 'finish'

//   publishDir "${params.output_folder}" 
  
  input:
    tuple val(name), file(reads)

  output:
  file "${name}/*"

"""
#!/bin/bash

set -Eeuo pipefail

fastqc -t ${task.cpus} *.fq

"""

}

// Run CheckM
process checkM {

  container "quay.io/biocontainers/checkm-genome:1.1.3--py_1"
  label "mem_medium"
  errorStrategy 'finish'

//   publishDir "${params.output_folder}" 
  
  input:
    file "input/*"

  output:
    file "output/*"

"""
#!/bin/bash

set -Eeuo pipefail

mkdir output

checkm lineage_wf input/ output/

"""

}

// Start the workflow
workflow {

    // Get the input files ending with BAM
    bam_ch = Channel.fromPath(
        "${params.input_folder}**{${params.suffix},${params.suffix}.pbi}"
    )

    // Extract the BAM files to FASTQ
    extractBAM(
        bam_ch
    )

    // Run FastQC on the reads
    fastQC(
        extractBAM.out
    )

    // Also run the assembler
    flye(
        extractBAM.out
    )

    // Check the quality of the assemblies
    checkM(
        flye.out
    )

}
