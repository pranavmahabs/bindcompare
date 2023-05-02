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
import pybedtools as bd
import sys
import subprocess

os.write(2, b"Completed BED Merge... starting downstream analysis!\n")

input_genes = sys.argv[1]
df = pd.read_csv(input_genes)
bed = df[['Chrom', 'Begin Bind Site', 'End Bind Site']]

fasta = bd.example_filename(sys.argv[2])
a = bd.BedTool.from_dataframe(bed)
a = a.sequence(fi=fasta)

data = open(a.seqfn).read().splitlines()
sequences = list(filter(lambda seq: seq[0] != '>', data))
df["Sequences"] = sequences
df = df.drop(columns=['Unnamed: 0'])
df.to_csv(input_genes)

outdir = sys.argv[3]
fastapath = outdir + '/sequences.fasta'
# Create the FASTA file of extracted sequences
fp = open(fastapath, 'w')
for _, row in df.iterrows():
    chrom_pos = row["Chrom"] + str(row['Begin Bind Site']) + ":" + str(row['End Bind Site'])
    name = ">" + row['Fbgns'].replace("\"", "").replace(";","") + ":" + chrom_pos + "\n"
    fp.write(name)
    sequence = row['Sequences'] + "\n"
    fp.write(sequence)
fp.close()

# Need to use STREME if num_sequences > 50:
if len(df) > 50:
    cmd_str = "streme -p " + fastapath + " --dna --nmotifs 10 --oc " + outdir + "/motifanalysis 2> /dev/null"
else:
    cmd_str = "meme -p" + fastapath +  " --dna --nmotifs 5 --maxw 50 --oc " + outdir + "/motifanalysis 2> /tmp/out.txt"
# Run STREME or MEME and produce the motif analysis
subprocess.run(cmd_str, shell=True, stdout=subprocess.DEVNULL)

os.write(2, b"Finished Sequence Extraction and MEME/STREME... Beginning GO!\n")