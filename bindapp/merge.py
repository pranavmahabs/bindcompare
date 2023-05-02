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
from utils import process_bed, process_gtf, average_peak_size, within

os.write(2, b"Beginning BindCompare!\n")

NUM_OVERLAPS = 0
unique_ref_overlaps = set()
unique_ola_overlaps = set()
# If reference is "RNA", then rna_loci is the reference loci and vice versa.
def compare_chrom_binds(dna_loci, rna_loci, reference, scope):
    """
    Returns the integer locations on a -1000 to 1000 range for any overlaps for overall tabulation.
    Also returns a list of all chrom locations that had the overlap. 
    """

    # Set the correct bed file for the reference.
    if reference == "RNA":
        ref_loci = rna_loci
        exp_loci = dna_loci
    else:
        ref_loci = dna_loci
        exp_loci = rna_loci

    # Set Up Return Structures
    chroms = ref_loci.keys()
    overlap_full = []
    overlap_front = []
    overlap_end = []
    overlap_locs = []
    overlap_ext = []

    # Make all the Binding Comparisons
    for chromosome in chroms:
        if chromosome not in exp_loci:
            continue
        for ref_bind in ref_loci[chromosome]:
            for exp_bind in exp_loci[chromosome]:
                overlap = within(exp_bind, ref_bind, scope)
                if len(overlap) == 0:
                    # No overlap in the binding sites... continue. 
                    continue
                else:
                    if exp_bind[0] >= ref_bind[0] and exp_bind[1] <= ref_bind[1]:
                        # The RNA peak is fully contained by the DNA peak
                        overlap_full.extend(overlap)
                        overlap_type = "Fully Contained"
                        ot = "OF"
                    elif exp_bind[1] > ref_bind[1] and exp_bind[0] <= ref_bind[1]:
                        # The RNA Peak overlaps the end of the DNA peak.
                        overlap_end.extend(overlap)
                        overlap_type = "Overlaps with End of Ref Peak"
                        ot = "OE"
                    elif exp_bind[0] < ref_bind[0] and exp_bind[1] >= ref_bind[0]:
                        # The RNA Peak overlaps the front of the DNA peak.
                        overlap_front.extend(overlap)
                        overlap_type = "Overlaps with Beg. of Ref Peak"
                        ot = "OB"
                    else:
                        overlap_ext.extend(overlap)
                        overlap_type = "Overlaps Outside of the Ref. Peak"
                        ot = "OX"
                    unique_ref_overlaps.add(ref_bind)
                    unique_ola_overlaps.add(exp_bind)
                    global NUM_OVERLAPS
                    NUM_OVERLAPS += 1
                    overlap_locs.append((chromosome, ref_bind[0], ref_bind[1], overlap_type, ot))

    return overlap_full, overlap_front, overlap_end, overlap_locs, overlap_ext

base_bed = sys.argv[1]
overlay_bed = sys.argv[2]
scope = sys.argv[3]
sample_name = sys.argv[4]
out_name = sys.argv[5]
processed_base_bed = process_bed(base_bed)
processed_overlay_bed = process_bed(overlay_bed)

full_o, front_o, end_o, overlap_coords, ext_o = compare_chrom_binds(processed_base_bed, processed_overlay_bed, "DNA", int(scope))

print(f'\nAverage Peak Sizes for {sample_name}:')
print(f'Base BED File:      {average_peak_size(processed_base_bed)} base pairs.')
print(f'Overlayed BED File: {average_peak_size(processed_overlay_bed)} base pairs.')

#setting up the array in numpy
overlap_full = np.asarray(full_o,dtype='int')
overlap_front = np.asarray(front_o,dtype='int')
overlap_end = np.asarray(end_o,dtype='int')
overlap_ext = np.asarray(ext_o,dtype='int')

hfont = {'fontname':'Sans Serif'}

fig = plt.figure()
ax1 = fig.add_subplot(111)

x ,y  = np.unique(overlap_full, return_counts=True) # counting occurrence of each loan
ax1.scatter(x, y, s=10, c='m', marker="s", label='Complete Peak Overlap')

x ,y  = np.unique(overlap_front, return_counts=True) # counting occurrence of each loan
ax1.scatter(x, y, s=5, c='r', marker="o", label='Overlap Ref Front')

x ,y  = np.unique(overlap_end, return_counts=True) # counting occurrence of each loan
ax1.scatter(x, y, s=5, c='b', marker="o", label='Overlaps Ref End')

x ,y  = np.unique(overlap_ext, return_counts=True) # counting occurrence of each loan
ax1.scatter(x, y, s=5, c='y', marker="o", label='External Overlaps')

plt.legend(loc='upper left')
plt.ylabel("Frequency", **hfont)
plt.xlabel("Overlap of Binding Sites", **hfont)
plt.title(f"Frequency of Binding Overlaps Over {2 * int(scope)} Base Pair Range", **hfont)
outpath = out_name + "/" + sample_name + "_" + "overlaps.png"
plt.savefig(outpath)
plt.close()

categories = ["Complete Overlap", "Partial Overlap Front", "Partial Overlap End", "External Overlap"]

abbr = ["OF", "OE", "OB", "OX"]
counts = [0, 0, 0, 0]
for olap in overlap_coords:
    index = abbr.index(olap[4])
    counts[index] = counts[index] + 1
fig1, ax1 = plt.subplots()

def pie_fmt(x):
    return '{:.0f}'.format((NUM_OVERLAPS)*x/100)

ax1.pie(counts, labels=categories, autopct=pie_fmt, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.title(f'Categorization of {NUM_OVERLAPS} Found Overlaps', **hfont)
outpath = out_name + "/" + sample_name + "_" + "pie.png"
plt.savefig(outpath)
plt.close()

categories = ["Total\nBinding Peaks", "Unique Overlaps", "Total Number\nof Overlaps", "Reference\nPeaks Overlapped"]
num_olay = 0
for rna_chrom in processed_overlay_bed.keys():
    num_olay = num_olay + len(processed_overlay_bed[rna_chrom])
vals = [num_olay, len(unique_ola_overlaps), NUM_OVERLAPS, len(unique_ref_overlaps)]
bars = plt.bar(categories, vals, color=['mediumblue', 'sandybrown', 'lightgrey', 'cornflowerblue'])

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + 0.3, yval + .3, yval, wrap=True)

plt.ylabel("")
plt.xlabel("")
plt.title('Total Number of Overlaps and Binding Peaks in Overlayed Bed', **hfont)
outpath = out_name + "/" + sample_name + "_" + "bartotals.png"
plt.savefig(outpath)
plt.close()

print("All plots in merge.py saved successfully.\n")

# # Gene Coordinate Section
gtf = sys.argv[6]
gene_coords = process_gtf(gtf)
all_geneids = set()
def find_fbgn(gtf, chrom, begin, end):
    overlaps = set()
    begin = int(begin) - int(scope)
    end = int(end) + int(scope)
    found = "Not in GTF"
    if chrom in gtf:
        for genes in gtf[chrom]:
            gene_beg = int(genes[1])
            gene_end = int(genes[2])
            if (begin >= gene_beg and begin < gene_end) or (end > gene_beg and end <= gene_end): 
                to_add = genes[0].replace("\"", "").replace(";", "")
                overlaps.add(to_add)
                all_geneids.add(to_add)
        if len(overlaps) != 0:
            found = " ".join(list(overlaps))
        return found
    else:
        return found

# overlap_locs.append((chromosome, ref_bind[0], ref_bind[1], overlap_type, ot))
df = pd.DataFrame(overlap_coords, columns=["Chrom", "Begin Bind Site", "End Bind Site", "Full Overlap Type", "Overlay_Type"])

df = df.groupby(df.columns.tolist(),as_index=False).size()
df.rename(columns={'size': 'Occurrences'}, inplace=True)
df.drop(columns=['Full Overlap Type'])
df['Occurrences'] = df['Occurrences'].apply(lambda x: str(x) + ",")
df['Overlay_Type'] = df['Overlay_Type'].apply(lambda x: str(x) + ",")

columns = ["Chrom", "Begin Bind Site", "End Bind Site", "Overlay Type", "Occurrences"]
df = df.groupby(["Chrom", "Begin Bind Site", "End Bind Site"]).agg(
     Overlay_Type = ('Overlay_Type','sum'),
     Occurrences = ('Occurrences','sum'),
     ).reset_index()

df['Occurrences'] = df['Occurrences'].apply(lambda x: x[:-1])
df['Overlay_Type'] = df['Overlay_Type'].apply(lambda x: x[:-1])
df['Fbgns'] = df.apply(lambda row: find_fbgn(gene_coords, row['Chrom'], row['Begin Bind Site'], row['End Bind Site']), axis=1)

outpath = out_name + "/" + sample_name + "_" + "overlaps.csv"
df.to_csv(outpath)

print("Here is a list of all Gene IDs in the overlaps CSV:")
for gene in all_geneids:
    print(gene, end=" ")
gene_arr = np.array(list(all_geneids))
df = pd.DataFrame(gene_arr.reshape(len(gene_arr), -1), columns=["Gene ID"])
print("\n")
outpath = "/tmp/gene_list.csv"
df.to_csv(outpath)