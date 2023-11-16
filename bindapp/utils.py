import textwrap


def process_gtf(gene_file):
    all_genes = {}
    with open(gene_file) as table:
        for line in table:
            gtf_row = line.split()
            if (gtf_row[0])[0:3] == "chr" or (gtf_row[0])[0:3] == "Chr":
                gtf_row[0] = (gtf_row[0])[3:]
            if (
                (gtf_row[0])[0:2] == "Un"
                or (len(gtf_row[0]) > 6)
                or (gtf_row[0] == "chr")
            ):
                continue
            else:
                this_chrom = gtf_row[0]
                # if not(str.isnumeric(bed_row[1])) or not(str.isnumeric(bed_row[2])):
                #     continue
                if this_chrom in all_genes:
                    all_genes[this_chrom].append((gtf_row[9], gtf_row[3], gtf_row[4]))
                else:
                    all_genes[this_chrom] = [(gtf_row[9], gtf_row[3], gtf_row[4])]
    return all_genes


def average_peak_size(processed_bed):
    total_sum, total_peaks = 0, 0
    for chrom in processed_bed:
        values = processed_bed[chrom]
        total_sum = total_sum + sum((int(j) - int(i)) for (i, j) in values)
        total_peaks = total_peaks + len(values)
    return total_sum / total_peaks


# Track the number of within/overlap calls made by the below function.
def within(exp_bind, ref_bind, scope):
    midpoint = (int)((int(ref_bind[0]) + int(ref_bind[1])) / 2)
    # Create the scope.
    overlap = range(
        max(int(exp_bind[0]), midpoint - scope),
        min(int(exp_bind[1]), midpoint + scope) + 1,
    )
    if len(overlap) == 0:
        return overlap
    else:
        # Adjust to -1000 to 1000 base pair range.
        overlap = [x - midpoint for x in overlap]
        return overlap


# Function from https://www.dunderdata.com/blog/Automatically%20Wrap%20Graph%20Labels%20in%20Matplotlib%20and%20Seaborn
def wrap_labels(ax, width, break_long_words=False):
    labels = []
    for label in ax.get_xticklabels():
        text = label.get_text()
        labels.append(
            textwrap.fill(text, width=width, break_long_words=break_long_words)
        )
    ax.set_xticklabels(labels, rotation=0)
