class Bed:
    def __init__(self, bedfile: str):
        self.bedfile = bedfile
        self.chroms = []
        self.loci = {}
        self.num_peaks = 0

    def update_chroms(self, chrom: str, start: str, end: str):
        """Update the chroms and loci attributes. Cast start and end to int."""
        start, end = int(start), int(end)
        self.num_peaks = self.num_peaks + 1
        if chrom not in self.chroms:
            self.chroms.append(chrom)
            self.loci[chrom] = [(start, end)]
        else:
            self.loci[chrom].append((start, end))

    def process_bed(self):
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
                        self.update_chroms(this_chrom, bed_row[1], bed_row[2])

    def get_chroms(self):
        """Returns a list of chromosomes."""
        return self.chroms

    def get_loci(self, chroms: list):
        """Returns a dictionary of loci for the given chromosomes."""
        loci = {}
        for chrom in chroms:
            loci[chrom] = self.loci[chrom]
        return loci

    def get_chrom(self, chrom: str):
        """Returns a list of loci for the given chromosome."""
        return self.loci[chrom]

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
            self.loci[chrom] = [(gene, start, end)]
        else:
            self.loci[chrom].append((gene, start, end))

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
                        self.update_chroms(
                            this_chrom, gtf_row[9], gtf_row[3], gtf_row[4]
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
        overlaps = set()
        begin = int(begin) - scope
        end = int(end) + scope
        found = "Not in GTF"
        if chrom in self.loci:
            for genes in self.loci[chrom]:
                gene_beg = int(genes[1])
                gene_end = int(genes[2])
                if (begin >= gene_beg and begin < gene_end) or (
                    end > gene_beg and end <= gene_end
                ):
                    to_add = genes[0].replace('"', "").replace(";", "")
                    overlaps.add(to_add)
                    all_genes.add(to_add)
            if len(overlaps) != 0:
                found = " ".join(list(overlaps))
            return found
        else:
            return found
