import numpy as np

import matplotlib.pyplot as plt

from multiprocessing import Process, Manager, Pool, cpu_count
import sys


class BedProcessor:
    def __init__(self, bed_file):
        self.bed_file = bed_file
        self.name = self.get_name_from_bed()
        self.global_dict = {}
        self.counts = {}
        self.total_sites = 0

    def get_name_from_bed(self):
        form = self.bed_file.find(".bed")
        if form == -1:
            print(f"Improper input {self.bed_file}. Please provide a BED (.bed) file.")
            quit()
        beg = self.bed_file.rfind("/")
        if beg != -1:
            name = self.bed_file[beg + 1 : form]
        else:
            name = self.bed_file[0:form]
        return name

    def process_bed(self, bin_size: int):
        with open(self.bed_file) as bed_obj:
            for line in bed_obj:
                bed_row = line.split()
                chrom = bed_row[0]

                if chrom[0:3] == "chr" or chrom[0:3] == "Chr":
                    chrom = (bed_row[0])[3:]
                if chrom[0:2] == "Un" or (len(chrom) > 5) or (chrom == "chr"):
                    continue
                if str.isnumeric(bed_row[1]) and str.isnumeric(bed_row[2]):
                    start = int(bed_row[1])
                else:
                    continue

                bin_id = int(start / bin_size)
                if chrom not in self.global_dict:
                    self.global_dict[chrom] = []
                    self.counts[chrom] = 0
                self.global_dict[chrom].append(bin_id)
                self.counts[chrom] += 1

        total_sites = 0
        for chrom in self.global_dict:
            self.counts[chrom] = len(self.global_dict[chrom])
            total_sites += self.counts[chrom]
            self.global_dict[chrom] = np.array(self.global_dict[chrom])

        self.total_sites = total_sites

    def __reduce__(self):
        return (self.__class__, (self.bed_file,))


class ProcessManager:
    def __init__(self, bed_files: list) -> None:
        self.manager = Manager()
        self.global_dict = self.manager.dict()
        self.processors = [BedProcessor(bf) for bf in bed_files]
        self.names = [bf.name for bf in self.processors]
        self.N = len(self.names)
        self.correl = np.ones((self.N, self.N))

    def process_bed_files(self, bin_size: int):
        for processor in self.processors:
            processor.process_bed(bin_size=bin_size)

    def overlap(self, args):
        i, j = args
        a = self.processors[i].global_dict
        b = self.processors[j].global_dict
        self.correl[i, j] = (
            sum(
                [
                    np.count_nonzero(np.isin(b[chrom], a[chrom])) if chrom in a else 0
                    for chrom in b
                ]
            )
            / self.processors[j].total_sites
        )
        self.correl[i, j]

    def plot_correlation_matrix(self):
        indices = [(i, j) for i in range(self.N) for j in range(self.N)]
        for index_pair in indices:
            self.overlap(index_pair)

        plt.figure(figsize=(10, 8))

        plt.imshow(self.correl, cmap="viridis", interpolation="nearest")
        plt.xticks(np.arange(self.N), self.names)
        plt.yticks(np.arange(self.N), self.names, rotation=90)
        plt.colorbar()

        # Create the Correlation Matrix
        plt.title(
            "Pair-wise Binding Overlap Frequencies",
            fontsize=16,
        )
        plt.xlabel("Reference Transcription Factor", fontsize=13)
        plt.ylabel("Overlapped Transcription Factors", fontsize=13)
        plt.savefig("explore.png")


def main():
    usage = "Usage: bindexplore <scope: int> <bed file 1> <bed file 2> ...\n"
    if len(sys.argv) == 1:
        print("\nError: No scope or BED file provided.\n" + usage)
        sys.exit(1)
    if not sys.argv[1].isnumeric():
        print("\nError: Integer scope not provided.\n" + usage)
        sys.exit(1)
    if len(sys.argv) < 4:
        print("\nError: Only 1 or less BED files provided.\n" + usage)
        sys.exit(1)

    scope = int(sys.argv[1])
    beds = sys.argv[2:]

    pm = ProcessManager(bed_files=beds)
    pm.process_bed_files(bin_size=scope)
    pm.plot_correlation_matrix()