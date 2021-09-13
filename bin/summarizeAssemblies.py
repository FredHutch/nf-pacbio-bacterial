#!/usr/bin/env python3

import gzip
import json
import pandas as pd
import sys

# Parse the input arguments
genome_name = sys.argv[1]
assembly_mode = sys.argv[2]
contigs_fasta_gz = sys.argv[3]
assembly_log = sys.argv[4]

# USEFUL FUNCTIONS

def parse_assembly_log(assembly_log):
    """Read through the UniCycler logs to get a list of circular contigs."""
    
    circular_contigs = []

    # Set a flag which we will use as we read through the file
    reading_circular_contigs = False

    with open(assembly_log, 'r') as handle:

        for line in handle:

            # If we are already in a list of circular contigs
            if reading_circular_contigs:

                # Try to read the contig from the line
                contig_name = parse_line(line)

                # If there is no contig name parsed
                if contig_name is None:

                    # Then we have ended the section
                    reading_circular_contigs = False

                # Otherwise, if the name is valid
                else:

                    # Add the contig name to the list
                    circular_contigs.append(contig_name)

            # Otherwise
            else:

                # If the line starts with "Segment"
                if line.startswith("Segment"):

                    # Then we have entered the indicated portion of the output
                    reading_circular_contigs = True

    return circular_contigs

def parse_line(l):
    """Read a single line of the UniCycler logs."""

    # Strip all of the whitespaces from both ends
    l = l.strip()

    # The line is empty
    if len(l) <= 1:
        return

    # The line has no spaces
    elif " " not in l:
        return

    # Otherwise
    else:
        # Return the string up until the first space
        return l.split(" ", 1)[0]

def parse_contig_info(contigs_fasta_gz, circular_contigs):
    """Count the size, depth, and circularity of each contig."""
    
    contig_info = []
    contig_lengths = dict()
    length_buffer = 0
    contig_name = None
    for l in gzip.open(contigs_fasta_gz, "rt"):
        if l.startswith(">"):
            if contig_name is not None:
                contig_lengths[contig_name] = length_buffer
                length_buffer = 0
            
            contig_name = l.lstrip(">").rstrip("\n").split()[0]
            contig_dict = dict(
                name=contig_name,
                circular=contig_name in circular_contigs
            )
            contig_info.append(contig_dict)

        else:
            length_buffer += len(l.rstrip("\n"))

    # Add the final contig
    contig_lengths[contig_name] = length_buffer

    # Make into a DataFrame
    contig_info = pd.DataFrame(
        contig_info
    )
    
    # Add the contig length
    contig_info = contig_info.assign(
        length=contig_info["name"].apply(contig_lengths.get)
    )

    # Sort it
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

    return contig_info, n50

# Get the list of circular contigs from the logs
circular_contigs = parse_assembly_log(assembly_log)

# Parse the FASTA itself
contig_info, n50 = parse_contig_info(contigs_fasta_gz, circular_contigs)

# Summarize these contigs
output = {
    "mode": assembly_mode,
    "genome_name": genome_name,
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

json.dump(output, open(f"{genome_name}.{assembly_mode}_mqc.json", "wt"), indent=4)
