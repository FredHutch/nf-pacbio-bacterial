# Process Raw PacBio Data from Bacterial Isolates

The purpose of this workflow is to process the raw PacBio data
generated from bacterial isolates, including:

- Extracting reads from BAM format (if necessary)
- Generating quality metrics for all specimens
- Performing _de novo_ assembly
- Measuring the completeness and contamination of each assembly

## Usage

To run this workflow, you will need to first:

- Install [Nextflow](nextflow.io) and configure it for your compute resources
- Identify the folder which contains all PacBio results
- Identify the folder in which all output files should be placed

NOTE: While reading input files, this script will traverse all subdirectories
and analyze any file ending with ".bam". The output for every individual input
file will be placed in a subdirectory within the output folder, with that
subdirectory named for the input file. Please ensure that all of these input
files have unique names.

### Basic Execution

The example immediately below provides the simplest example of how to run
this workflow, without including any of the advanced commands for controlling
the compute resources used to perform the task.

```#!/bin/bash

set -Eeuo pipefail

nextflow \
    run \
    FredHutch/nf-pacbio-bacterial \
    --input-folder <INPUT FOLDER> \
    --output-folder <OUTPUT FOLDER> \
    -resume \
    -with-report \
    -with-trace
```

### Advanced Execution

The example below contains a number of additional optional flags
or parameters which can be used to customize the workflow execution.

```#!/bin/bash

set -Eeuo pipefail

AWS_PROFILE=$AWS_PROFILE \
NXF_VER=20.10.0 \
nextflow \
    run \
    -c <NEXTFLOW CONFIG FILG> \
    FredHutch/nf-pacbio-bacterial \
    --input-folder <INPUT FOLDER> \
    --output-folder <OUTPUT FOLDER> \
    -resume \
    -latest \
    -with-report \
    -with-trace \
    -process.queue <COMPUTE QUEUE>
```
