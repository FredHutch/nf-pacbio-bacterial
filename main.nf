#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Containers
container__unicycler = "quay.io/biocontainers/unicycler:0.4.7--py37hdbcaa40_1"
container__pandas = "quay.io/fhcrc-microbiome/python-pandas:latest"
container__fastqc = "quay.io/biocontainers/fastqc:0.11.9--0"
container__multiqc = "quay.io/biocontainers/multiqc:1.10--py_1"
container__checkm = "quay.io/biocontainers/checkm-genome:1.1.3--py_1"
container__prokka = "quay.io/fhcrc-microbiome/prokka:latest"
container__lima = "quay.io/biocontainers/lima:2.2.0--h9ee0642_0"
container__busco = "ezlabgva/busco:v5.2.2_cv1"

// Function which prints help message text
def helpMessage() {
    log.info"""
    Usage:

    nextflow run FredHutch/nf-pacbio-bacterial <ARGUMENTS>
    
    Required Arguments:
      --input_folder        Folder containing all PacBio data in FASTQ files (including subdirectories)
      --output_folder       Folder to place analysis outputs

    Input Files:
      --suffix              Process all files ending with this string (default: .fastq.gz)
      --barcodes            Optional FASTA file with per-sample barcode sequences. 
                            If specified, input FASTQ files will be demultiplexed as described:
                            https://lima.how/get-started.html

    Optional Assembly Arguments (passed UniCycler):
      --mode                Assembly mode
                            Default: normal
                            Options: normal | bold | conservative
    
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
params.suffix = ".fastq.gz"
params.barcodes = false
params.demux_flags = "--same --split-named"
params.mode = "normal"


/////////////////////
// DEFINE FUNCTIONS /
/////////////////////

// Run lima
process demultiplex {
    container "${container__lima}"
    label "mem_medium"
    publishDir "${params.output_folder}/demux_reads/", mode: "copy", overwrite: true, pattern: "${genome_name}.demux.*.fq.gz"
    publishDir "${params.output_folder}/demux_report/", mode: "copy", overwrite: true, pattern: "${genome_name}.demux.lima.*"

    input:
    tuple val(genome_name), file(input_fastq)
    file barcodes_fasta

    output:
    path "${genome_name}.demux.*.fq.gz", emit: reads
    path "${genome_name}.demux.lima.*", emit: reports

"""
set -Eeuo pipefail

lima \
    "${input_fastq}" \
    "${barcodes_fasta}" \
    "${genome_name}.demux.fq.gz" \
    ${params.demux_flags}

"""
}

// Run UniCycler
process unicycler {
    container "${container__unicycler}"
    label "mem_veryhigh"
    publishDir "${params.output_folder}/${genome_name}/${params.mode}/", mode: "copy", overwrite: true

    input:
    tuple val(genome_name), file(long_reads)

    output:
    tuple val(genome_name), file("${genome_name}/${genome_name}.fasta.gz")
    path "${genome_name}/${genome_name}.gfa"
    path "${genome_name}/${genome_name}.log"

"""
set -Eeuo pipefail

mkdir ${genome_name}

unicycler \
    -l ${long_reads} \
    -o ${genome_name} \
    --keep 0 \
    --mode ${params.mode} \
    -t ${task.cpus}

mv ${genome_name}/assembly.gfa ${genome_name}/${genome_name}.gfa
mv ${genome_name}/assembly.fasta ${genome_name}/${genome_name}.fasta
mv ${genome_name}/unicycler.log ${genome_name}/${genome_name}.log
gzip ${genome_name}/${genome_name}.fasta
"""
}


// Summarize the assembly
process summarizeAssemblies {
    container "${container__pandas}"
    label "io_limited"
    publishDir "${params.output_folder}/${genome_name}/${params.mode}/", mode: "copy", overwrite: true

    input:
    tuple val(genome_name), file(contigs_fasta_gz)

    output:
    path "${genome_name}.${params.mode}.json"

"""#!/usr/bin/env python3
import gzip
import json
import pandas as pd

# Count the size, depth, and circularity of each contig
contig_info = []
contig_lengths = dict()
length_buffer = 0
contig_name = None
for l in gzip.open("${contigs_fasta_gz}", "rt"):
    if l.startswith(">"):
        if contig_name is not None:
            contig_lengths[contig_name] = length_buffer
            length_buffer = 0
        if " " in l:
            contig_name, contig_dict = l.lstrip(">").rstrip("\\n").split(" ", 1)
            contig_dict = dict([
                (i.split("=")[0], i.split("=")[1])
                for i in contig_dict.split(" ")
            ])
            contig_dict["name"] = contig_name
            contig_dict["circular"] = contig_dict.get("circular", "false")
            contig_info.append(contig_dict)
        else:
            contig_name = l.lstrip(">").rstrip("\\n")
    else:
        length_buffer += len(l.rstrip("\\n"))
# Add the final contig
contig_lengths[contig_name] = length_buffer

# Make into a DataFrame
if len(contig_info) > 0:
    contig_info = pd.DataFrame(contig_info)
    contig_info["length"] = contig_info["length"].apply(int)
    contig_info["depth"] = contig_info["depth"].apply(lambda s: float(s.rstrip("x")))
    contig_info["circular"] = contig_info["circular"].fillna("false") == "true"
else:
    contig_info = pd.DataFrame(dict([("length", contig_lengths)])).reset_index()
    contig_info["depth"] = 1
    contig_info["circular"] = False
contig_info.sort_values(by="length", ascending=False, inplace=True)

# Calculate N50
running_total = 0
n50 = None
for nbases in contig_info["length"].values:
    running_total += nbases
    if running_total >= contig_info["length"].sum() / 2.:
        n50 = int(nbases)
        break
assert n50 is not None

# Summarize these contigs
output = {
    "mode": "${params.mode}",
    "genome_name": "${genome_name}",
    "num_contigs": int(contig_info.shape[0]),
    "num_circular_contigs": int(contig_info["circular"].sum()),
    "circular_contig_lengths": ", ".join(contig_info.query("circular")["length"].apply(str).tolist()),
    "linear_contig_lengths": ", ".join(contig_info.query("circular == False")["length"].apply(str).tolist()),
    "num_bases": int(contig_info["length"].sum()),
    "longest_contig": int(contig_info["length"].max()),
    "num_over_1Mb": int((contig_info["length"] >= 1000000).sum()),
    "num_100kb_to_1Mb": int(((contig_info["length"] < 1000000) & (contig_info["length"] >= 100000)).sum()),
    "num_10kb_to_100kb": int(((contig_info["length"] < 100000) & (contig_info["length"] >= 10000)).sum()),
    "num_1kb_to_10kb":    int(((contig_info["length"] < 10000) & (contig_info["length"] >= 1000)).sum()),
    "num_under_1kb":        int((contig_info["length"] < 1000).sum()),
    "N50": n50,
}

json.dump(output, open("${genome_name}.${params.mode}.json", "wt"), indent=4)

"""
}


// Combine the genome summaries
process combineSummaries {
    container "${container__pandas}"
    label "io_limited"
    publishDir "${params.output_folder}/", mode: "copy", overwrite: true

    input:
    file "*"

    output:
    file "assembly_summary.${params.mode}.csv"

"""#!/usr/bin/env python3

import json
import os
import pandas as pd

# Make a DataFrame
df = pd.DataFrame(
    [
        # Containing data read from the file in JSON format
        json.load(open(fp, 'r'))
        # Reading over each file in the folder
        for fp in os.listdir(".")
        # If the file ends with the expected suffix
        if fp.endswith(".${params.mode}.json")
    ]
)

# Write it out to a file
df.to_csv("assembly_summary.${params.mode}.csv", index=None)

"""
}

// Run FastQC
process fastQC {

  // Docker container to use
  container "${container__fastqc}"
  label "io_limited"

  publishDir "${params.output_folder}/${name}/${params.mode}/", mode: "copy", overwrite: true
  
  input:
    tuple val(name), file(reads)

  output:
  path "${name}/*"
  path "${name}/*_fastqc.zip", emit: zip

"""
#!/bin/bash

set -Eeuo pipefail

mkdir ${name}

fastqc -t ${task.cpus} -o ${name} ${reads}

ls -lahtr

"""

}

// Run MultiQC
process multiQC {

  container "${container__multiqc}"

  publishDir "${params.output}", mode: 'copy', overwrite: true
  
  input:
    file "*_fastqc.zip"

  output:
    file "multiqc_report.html"

"""
#!/bin/bash

set -Eeuo pipefail

multiqc .

ls -lahtr
"""

}


// Run CheckM
process checkM {

  container "${container__checkm}"
  label "mem_medium"

  publishDir "${params.output_folder}/${name}/${params.mode}/checkm/", mode: "copy", overwrite: true 
  
  input:
    tuple val(name), file(fasta_gz)

  output:
    file "*"

"""
#!/bin/bash

set -Eeuo pipefail

mkdir input
mkdir output
gunzip -c ${fasta_gz} > input/${name}.fa
rm ${fasta_gz}

checkm lineage_wf -t ${task.cpus} -x fa input/ output/

mv output/* ./
rmdir output
rm -r input

"""

}

// Run BUSCO
process busco {

  container "${container__busco}"
  label "mem_medium"

  publishDir "${params.output_folder}/${name}/${params.mode}/busco/", mode: "copy", overwrite: true 
  
  input:
    file faa_gz

  output:
    file "*"

"""
#!/bin/bash

set -Eeuo pipefail

# Get the base file name
BASE_NAME="${faa_gz.name.replaceAll(/.faa.gz/, "")}"

# Decompress the input
gunzip -c "${faa_gz}" > \$BASE_NAME.faa

# Run BUSCO
busco \
    -m protein \
    -i  "\$BASE_NAME.faa" \
    -o "\$BASE_NAME" \
    --auto-lineage-prok

# Remove the temporary decompressed copy of the input files
rm \$BASE_NAME.faa

"""

}

process prokka {
    container "${container__prokka}"
    label "mem_medium"
    publishDir "${params.output_folder}/${genome_name}/${params.mode}/", mode: "copy", overwrite: true

    input:
    tuple val(genome_name), file(contigs_fasta_gz)
    
    output:
    path "prokka/${genome_name}.${params.mode}.faa.gz", emit: faa
    path "prokka/${genome_name}.${params.mode}.gbk.gz"
    path "prokka/${genome_name}.${params.mode}.gff.gz"
    path "prokka/${genome_name}.${params.mode}.tsv.gz"
    
    """#!/bin/bash

set -Eeuo pipefail

# Decompress the assembly
gunzip -c ${contigs_fasta_gz} > ${contigs_fasta_gz.simpleName}.fasta

prokka \
    --outdir prokka/ \
    --prefix ${genome_name}.${params.mode} \
    --cpus ${task.cpus} \
    --compliant \
    ${contigs_fasta_gz.simpleName}.fasta


gzip prokka/${genome_name}.${params.mode}.faa
gzip prokka/${genome_name}.${params.mode}.gbk
gzip prokka/${genome_name}.${params.mode}.gff
gzip prokka/${genome_name}.${params.mode}.tsv
    """
}

// Start the workflow
workflow {

    // Get the input files ending with BAM
    input_ch = Channel.fromPath(
        "${params.input_folder}**${params.suffix}"
    ).map {
        it -> [it.name.replaceAll(/${params.suffix}/, ''), it]
    }

    // If barcodes were provided
    if (params.barcodes) {

        // Demultiplex the reads
        demultiplex(
            input_ch,
            file(params.barcodes)
        )

        // Set up a channel with the demultiplexed reads
        fastq_ch = demultiplex
            .out
            .reads
            .flatten()
            .map {
                it -> [it.name.replaceAll(/.fq.gz/, ''), it]
            }


    } else {

        // Run the downstream analysis on the input reads
        fastq_ch = input_ch

    }

    // Run the assembler on the input reads
    unicycler(
        fastq_ch
    )

    // Run FastQC on the input reads
    fastQC(
        fastq_ch
    )

    // Compute metrics on the assembly
    summarizeAssemblies(
        unicycler.out[0]
    )

    // Join together the assembly metrics
    combineSummaries(
        summarizeAssemblies.out.toSortedList()
    )

    // Check the quality of the assemblies
    checkM(
        unicycler.out[0]
    )

    // Annotate the assemblies
    prokka(
        unicycler.out[0]
    )

    // Score the completeness of the coding gene content
    busco(
        prokka.out.faa
    )

}