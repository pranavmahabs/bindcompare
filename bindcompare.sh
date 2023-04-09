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

DNA=$1
RNA=$2
SCOPE=$3
SNAME=$4
OUT=$5
GTF=$6
FASTA=$7

python3 merge.py ${DNA} ${RNA} ${SCOPE} ${SNAME} ${OUT} ${GTF} > ${OUT}/${SNAME}_summary.txt

python3 downstream.py ${OUT}/${SNAME}_genes.csv ${FASTA} ${OUT} >> ${OUT}/${SNAME}_summary.txt

Rscript geneont.R --outdir ${OUT} --expcondition ${SNAME}

# EXAMPLE TERMINAL CODE
# ./bindcompare.sh /Users/pranavmahableshwarkar/BrownUniversity/LarschanLab/Analytics/PM_NewScripts/CNR_bedfiles/KC_consensusPeaks.bed /Users/pranavmahableshwarkar/BrownUniversity/LarschanLab/Analytics/PM_NewScripts/iCLIP_bedfiles/Kc_ChF.bed 750 KC SampleOut /Users/pranavmahableshwarkar/BrownUniversity/LarschanLab/Analytics/reference/dmel-all-r6.42.gtf /Users/pranavmahableshwarkar/BrownUniversity/LarschanLab/Analytics/reference/dm6.fa
