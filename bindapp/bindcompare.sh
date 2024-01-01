#!/bin/bash

###############
# BindCompare: bindcompare.sh
# Author: Pranav Mahableshwarkar
# Date: 04/12/2023
# Note: Runs the BindCompare analysis pipeline. 
###############

if [ $# -ne 7 ]; then
	echo $0: "Usage: 'bash preprocess.sh'
            1) DNA Bed File Path
            2) RNA Bed File Path
            3) Scope (How many nucleotides to look in each direction from the middle of the DNA ref peak)
            4) Sample Name
            5) Output Directory
            6) All Genes GTF File
            7) All Genome fasta File
            "
        exit 1
fi

# eval "$(conda shell.bash hook)"
# conda activate bindcompare

DNA=$1
RNA=$2
SCOPE=$3
SNAME=$4
OUT=$5
GTF=$6
FASTA=$7

# echo "Beginning BindCompare! Any errors will be printed in the Summary File or Command Line."
python3 bindapp/merge.py ${DNA} ${RNA} ${SCOPE} ${SNAME} ${OUT}/ ${GTF} > ${OUT}/${SNAME}_summary.txt

# # echo "Completed Merge, Beginning Downstream Analysis!"
python3 bindapp/downstream.py ${OUT}/${SNAME}_overlaps.csv ${FASTA} ${OUT} >> ${OUT}/${SNAME}_summary.txt

# # echo "Finished Motif Analysis, Beginning Gene Ontology!"
# Rscript bindapp/geneont.R --outdir ${OUT} --expcondition ${SNAME}

echo "Completed BindCompare! Time Stamp:"
date

rm Rplots.pdf 2> /dev/null

# EXAMPLE TERMINAL CODE
# ./bindcompare.sh ~/CNR_bedfiles/KC_consensusPeaks.bed ~/iCLIP_bedfiles/Kc_ChF.bed 1000 KC ~/bindcompare/KCSampleOut ~/reference/dmel-all-r6.42.gtf ~/reference/dm6.fa
