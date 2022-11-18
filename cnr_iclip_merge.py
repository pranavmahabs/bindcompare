import sys
import numpy as np
import matplotlib.pyplot as plt

def process_bed(bedfile):
    # Read in table using dask library.
    reference_peaks = {}
    with open(bedfile) as table:
        for line in table:
            bed_row = line.split()
            if (bed_row[0])[0:3] != "chr" and (bed_row[0])[0:3] != "Chr":
                bed_row[0] = "chr" + bed_row[0]
            if (bed_row[0])[0:5] == "chrUn" or (len(bed_row[0]) > 6) or (bed_row[0] == "chr"):
                continue
            else:
                this_chrom = bed_row[0]
                if not(str.isnumeric(bed_row[1])) or not(str.isnumeric(bed_row[2])):
                    continue
                elif this_chrom in reference_peaks:
                    reference_peaks[this_chrom].append((bed_row[1], bed_row[2], bed_row[4]))
                else:
                    reference_peaks[this_chrom] = [(bed_row[1], bed_row[2], bed_row[4])]
    return reference_peaks

def within(exp_bind, ref_bind, scope):
    midpoint = (int)((int(ref_bind[0]) + int(ref_bind[1])) / 2)
    # Create the scope. 
    overlap = range(max(int(exp_bind[0]), midpoint - scope), min(int(exp_bind[-2]), midpoint + scope)+1)
    if len(overlap) == 0:
        return overlap
    else:
        # Adjust to -1000 to 1000 base pair range. 
        overlap = [x - midpoint for x in overlap]
        return overlap
# Track the number of within/overlap calls made by the below function. 
within.counter = 0

# If reference is "RNA", then rna_loci is the reference loci and vice versa.
def compare_chrom_binds(dna_loci, rna_loci, reference, scope, count_more):
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
                    # Insignificant p-value... continue if p-value is considered.
                    # No overlap in the binding sites... continue. 
                    continue
                else:
                    overlap_type = "None"
                    ot = "None"
                    # overlap_all.extend(overlap)
                    if exp_bind[0] >= ref_bind[0] and exp_bind[1] <= ref_bind[1]:
                        # The RNA peak is fully contained by the DNA peak
                        within.counter = within.counter + 1
                        overlap_full.extend(overlap)
                        overlap_type = "Fully Contained"
                        ot = "OF"
                    elif exp_bind[1] > ref_bind[1] and exp_bind[0] <= ref_bind[1]:
                        # The RNA Peak overlaps the end of the DNA peak.
                        within.counter = within.counter + 1
                        overlap_end.extend(overlap)
                        overlap_type = "Overlaps with End of Ref Peak"
                        ot = "OE"
                    elif exp_bind[0] < ref_bind[0] and exp_bind[1] >= ref_bind[0]:
                        # The RNA Peak overlaps the front of the DNA peak.
                        within.counter = within.counter + 1
                        overlap_front.extend(overlap)
                        overlap_type = "Overlaps with Beg. of Ref Peak"
                        ot = "OB"
                    else:

                        within.counter = within.counter + 1
                        overlap_ext.extend(overlap)
                        if exp_bind[1] < ref_bind[0]:
                            overlap_type = "Before Ref Peak"
                            ot = "BR"
                        elif exp_bind[0] > ref_bind[1]:
                            overlap_type = "After Ref Peak"
                            ot = "AR"
                    overlap_locs.append((chromosome, ref_bind[0], ref_bind[1], overlap_type, ot))
                    if not count_more:
                        break

    return overlap_full, overlap_front, overlap_end, overlap_locs, overlap_ext

dna = sys.argv[1]
rna = sys.argv[2]
scope = sys.argv[3]
sample_name = sys.argv[4]
out_name = sys.argv[5]

processed_dna_bed = process_bed(dna)
processed_rna_bed = process_bed(rna)

overlap_full, overlap_front, overlap_end, overlap_locs, overlap_ext = compare_chrom_binds(processed_dna_bed, processed_rna_bed, "DNA", int(scope), True)

#setting up the array in numpy
overlap_full = np.asarray(overlap_full,dtype='int')
overlap_front = np.asarray(overlap_front,dtype='int')
overlap_end = np.asarray(overlap_end,dtype='int')
overlap_ext = np.asarray(overlap_ext,dtype='int')


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
plt.ylabel("Frequency")
plt.xlabel("Overlap of Binding Sites")
plt.title(f"Frequency of Binding Overlaps Over {2 * int(scope)} Base Pair Range")
outpath = out_name + "/" + sample_name + "_" + "overlaps.png"
plt.savefig(outpath)
plt.close()

categories = ["Complete Overlap", "Partial Overlap Front", "Partial Overlap End", "Before Reference Peak", "After Reference Peak"]
abbr = ["OF", "OE", "OB", "BR", "AR"]
counts = [0, 0, 0, 0, 0]
for olap in overlap_locs:
    index = abbr.index(olap[4])
    counts[index] = counts[index] + 1
fig1, ax1 = plt.subplots()

def pie_fmt(x):
    return '{:.0f}'.format((within.counter)*x/100)

ax1.pie(counts, labels=categories, autopct=pie_fmt, startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
plt.title(f'Categorization of {within.counter} Found Overlaps')
outpath = out_name + "/" + sample_name + "_" + "pie.png"
plt.savefig(outpath)
plt.close()

categories = ["Overlaps", "DNA Binding Sites", "RNA Binding Sites"]
num_rna = 0
num_dna = 0
for rna_chrom in processed_rna_bed.keys():
    num_rna = num_rna + len(processed_rna_bed[rna_chrom])
for dna_chrom in processed_dna_bed.keys():
    num_dna = num_dna + len(processed_dna_bed[dna_chrom])
vals = [within.counter, num_dna, num_rna]
bars = plt.bar(categories, vals)

for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + 0.3, yval + .25, yval)

plt.ylabel("")
plt.xlabel("")
plt.title('Number of Overlaps and Total Binding Sites in DNA/RNA')
outpath = out_name + "/" + sample_name + "_" + "bartotals.png"
plt.savefig(outpath)
plt.close()










