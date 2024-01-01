from intervaltree import Interval, IntervalTree


class Bed:
    def __init__(self, bedfile: str):
        self.bedfile = bedfile
        self.chroms = []
        self.loci = {}
        self.loci_it = {}
        self.num_peaks = 0

    def update_chroms(self, chrom: str, start: str, end: str, is_ref: bool, scope: int):
        """Update the chroms and loci attributes. Cast start and end to int."""
        start, end = int(start), int(end)
        self.num_peaks = self.num_peaks + 1
        if chrom not in self.chroms:
            self.chroms.append(chrom)
            self.loci[chrom] = [(start, end)]
            self.loci_it[chrom] = IntervalTree()
        else:
            self.loci[chrom].append((start, end))

        midpoint = int((start + end) / 2)
        if is_ref:
            self.loci_it[chrom][midpoint - scope : midpoint + scope + 1] = (start, end)
            # "{}:{}-{}".format(chrom, start, end)
        else:
            self.loci_it[chrom][start:end] = (start, end)
            # "{}:{}-{}".format(chrom, start, end)

    def process_bed(self, is_ref: bool, scope: int):
        """Process the bedfile and update the chroms and loci attributes."""
        reference_peaks = {}
        with open(self.bedfile) as table:
            for line in table:
                bed_row = line.split()
                if (bed_row[0])[0:3] == "chr" or (bed_row[0])[0:3] == "Chr":
                    bed_row[0] = (bed_row[0])[3:]
                if (
                    (bed_row[0])[0:2] == "Un"
                    or (len(bed_row[0]) > 5)
                    or (bed_row[0] == "chr")
                ):
                    continue
                else:
                    this_chrom = bed_row[0]
                    if str.isnumeric(bed_row[1]) and str.isnumeric(bed_row[2]):
                        self.update_chroms(
                            this_chrom, bed_row[1], bed_row[2], is_ref, scope
                        )

    def get_chroms(self):
        """Returns a list of chromosomes."""
        return self.chroms

    def get_loci(self, chroms: list):
        """Returns a dictionary of loci for the given chromosomes."""
        loci = {}
        for chrom in chroms:
            loci[chrom] = self.loci[chrom]
        return loci

    def get_chrom(self, chrom: str) -> list:
        """Returns a list of loci for the given chromosome."""
        return self.loci[chrom]

    def get_chrom_it(self, chrom: str) -> IntervalTree:
        """Returns a interval tree of loci for the given chromosome."""
        return self.loci_it[chrom]

    def average_peak_size(self, chroms: list):
        """Returns the average peak size for the given chromosomes."""
        total_sum, total_peaks = 0, 0
        for chrom in chroms:
            values = self.loci[chrom]
            total_sum = total_sum + sum((int(j) - int(i)) for (i, j) in values)
            total_peaks = total_peaks + len(values)
        return total_sum / total_peaks


class GTF:
    def __init__(self, gtf_file: str):
        self.gtf_file = gtf_file
        self.chroms = []
        self.loci = {}

    def update_chroms(self, chrom: str, gene: str, start: str, end: str):
        """Update the chroms and loci attributes. Cast start and end to int."""
        start, end = int(start), int(end)
        if chrom not in self.chroms:
            self.chroms.append(chrom)
            # self.loci[chrom] = [(gene, start, end)]
            self.loci[chrom] = IntervalTree()
        # self.loci[chrom].append((gene, start, end))
        if start == end:
            end += 1
        self.loci[chrom][start:end] = gene

    def process_gtf(self):
        """Process the gtf file and update the chroms and loci attributes."""
        with open(self.gtf_file) as table:
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
                    if str.isnumeric(gtf_row[3]) and str.isnumeric(gtf_row[4]):
                        try:
                            gene_index = gtf_row.index("gene_id")
                            if gene_index == len(gtf_row) - 1:
                                continue
                        except ValueError:
                            print('"gene_id" not found in GTF attributes column.')
                        self.update_chroms(
                            this_chrom, gtf_row[gene_index + 1], gtf_row[3], gtf_row[4]
                        )

    def get_chroms(self):
        """Returns a list of chromosomes."""
        return self.chroms

    def get_loci(self, chroms: list):
        """Returns a dictionary of loci for the given chromosomes."""
        loci = {}
        for chrom in chroms:
            loci[chrom] = self.loci[chrom]
        return loci

    def find_fbgn(self, chrom: str, begin: str, end: str, scope: int, all_genes: set):
        """Returns the gene ID for the given chromosome, begin, and end."""
        # print(f"Finding FBGN for {chrom}:{begin}-{end}")
        genes = set()
        begin = int(begin) - scope
        end = int(end) + scope
        found = "Not in GTF"
        if chrom in self.loci:
            overlaps = self.loci[chrom][begin:end]
            if len(overlaps) == 0:
                return found
            for gene_interval in overlaps:
                _, _, gene_id = gene_interval
                to_add = gene_id.replace('"', "").replace(";", "")
                genes.add(to_add)
                all_genes.add(to_add)
            if len(genes) != 0:
                found = " ".join(list(genes))
            return found
        else:
            return found
