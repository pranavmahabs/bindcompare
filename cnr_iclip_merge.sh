#!/bin/bash

if [ $# -ne 5 ]; then
	echo $0: "Usage: 'sbatch preprocess.sh' or 'bash preprocess.sh'
            1) DNA Bed File Path
            2) RNA Bed File Path
            3) Scope (How many nucleotides to look in each direction from the middle of the DNA ref peak)
            4) Sample Name
            5) Output Directory
            "
        exit 1
fi

DNA=$1
RNA=$2
SCOPE=$3
SNAME=$4
OUT=$5

python3 cnr_iclip_merge.py ${DNA} ${RNA} ${SCOPE} ${SNAME} ${OUT} > ${OUT}/${SNAME}_summary.txt

