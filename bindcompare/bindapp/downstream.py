###############
# BindCompare: downstream.py
# Author: Pranav Mahableshwarkar
# Date: 04/12/2023
# Note: Create a formatted CSV for looking through results.
#       Extracts sequences from overlapping binding sites.
#       Performs MEME/STREME motif analysis.
# Output: MEME/STREME Folder of Results, Overlaps.CSV, sequences.fasta.
###############

import os
import pandas as pd
from Bio import SeqIO
import sys
import subprocess


def get_chrom2seq(FASTA_FILE, capitalize=True):
    """
    Load in the genome fasta file to extract sequences from the BED Files.
    """
    chrom2seq = {}
    for seq in SeqIO.parse(FASTA_FILE, "fasta"):
        chrom2seq[seq.description.split()[0]] = (
            seq.seq.upper() if capitalize else seq.seq
        )
    return chrom2seq


def downstream(input_genes: str, fasta: str, outdir: str):
    if fasta == "None":
        os.write(2, b"Skipping Sequence Extraction...\n")
    else:
        os.write(2, b"Completed BED Merge... starting sequence extraction!\n")
        df = pd.read_csv(input_genes)
        bed = df[["Chrom", "Begin Ref Site", "End Ref Site"]]

        chrom2seq = get_chrom2seq(fasta)

        df["Sequences"] = df.apply(
            lambda row: str(
                chrom2seq[row.Chrom][row["Begin Ref Site"] : row["End Ref Site"]]
            ),
            axis=1,
        )
        df.to_csv(input_genes)

        fastapath = outdir + "/sequences.fasta"
        # Create the FASTA file of extracted sequences
        fp = open(fastapath, "w")
        for _, row in df.iterrows():
            chrom_pos = (
                row["Chrom"]
                + ":"
                + str(row["Begin Ref Site"])
                + "-"
                + str(row["End Ref Site"])
            )
            name = (
                ">"
                + row["GeneIDs"].replace('"', "").replace(";", "")
                + "; "
                + chrom_pos
                + "\n"
            )
            fp.write(name)
            sequence = row["Sequences"] + "\n"
            fp.write(sequence)
        fp.close()

        # Need to use STREME if num_sequences > 50:
        if len(df) > 50:
            cmd_str = (
                "streme -p "
                + fastapath
                + " --dna --nmotifs 10 --oc "
                + outdir
                + "/motifanalysis 2> /dev/null"
            )
        else:
            cmd_str = (
                "meme -p"
                + fastapath
                + " --dna --nmotifs 5 --maxw 50 --oc "
                + outdir
                + "/motifanalysis 2> /tmp/out.txt"
            )
        # Run STREME or MEME and produce the motif analysis
        # subprocess.run(cmd_str, shell=True, stdout=subprocess.DEVNULL)

        # os.write(2, b"Finished Sequence Extraction and MEME/STREME... Beginning GO!\n")


if __name__ == "__main__":
    downstream()
