###############
# BindCompare: merge.py
# Author: Pranav Mahableshwarkar
# Date: 04/12/2023
# Note: Merge/Prep the Two BED Files for downstream analysis.
#       Create overlap profile and initial visualizations.
# Output: Bar Totals, Pie Chart, Summary.txt, Overlaps PNG
#         List of Identified Genes
###############

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# from utils import process_bed, process_gtf, average_peak_size, within
from merge_class import Bed, GTF
from exp_class import BindCompare

os.write(2, b"Beginning BindCompare!\n")

base_bed = sys.argv[1]
overlay_bed = sys.argv[2]
scope = sys.argv[3]
sample_name = sys.argv[4]
out_name = sys.argv[5]


# Initialize the BED Files
base_bed = Bed(base_bed)
overlay_bed = Bed(overlay_bed)

# Process the BED Files
base_bed.process_bed(True, int(scope))
overlay_bed.process_bed(False, int(scope))

# Set up the BindCompare Experiment
exp = BindCompare(base_bed, overlay_bed, int(scope))
exp.compare_binds()

# Get the Chromosomes
b_chroms = base_bed.get_chroms()
e_chroms = overlay_bed.get_chroms()


print(f"\nAverage Peak Sizes for {sample_name}:")
print(f"Base BED File:      {base_bed.average_peak_size(b_chroms)} base pairs.")
print(f"Overlayed BED File: {overlay_bed.average_peak_size(e_chroms)} base pairs.")

# # Gene Coordinate Section
gtf = sys.argv[6]
gtf = GTF(gtf)
gtf.process_gtf()

# Get the BC Dictionary for All Chromosomes
exp.compare_binds()
bc_it = exp.get_experiments_overlaps_it(b_chroms)

# Perform all Plotting
exp.generate_all(bc_it, out_name, sample_name, gtf)
