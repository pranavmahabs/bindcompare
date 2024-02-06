import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from intervaltree import Interval, IntervalTree
import pandas as pd
from .merge_class import Bed, GTF
import os

import time


def moving_average(data, window_size):
    smoothed_data = np.zeros_like(data)
    half_window = window_size // 2

    for i in range(half_window, len(data) - half_window):
        smoothed_data[i] = np.mean(data[i - half_window : i + half_window + 1])

    # Pad the edges of the smoothed data with initial and final values
    for i in range(half_window):
        smoothed_data[i] = np.mean(data[: i + half_window + 1])
        smoothed_data[len(data) - 1 - i] = np.mean(
            data[len(data) - i - half_window - 1 :]
        )

    return smoothed_data


class Chromosome:
    def __init__(
        self,
        name: str,
        scope: int,
        ref_it: IntervalTree,
        exp_it: IntervalTree,
    ):
        self.chr = name
        self.ref_binds_it = ref_it
        self.exp_binds_it = exp_it
        self.scope = scope

        self.overlap_counts_it = {"Total": 0, "CRO": 0, "ORE": 0, "ORF": 0, "PXP": 0}
        self.overlaps_it = []
        self.overlap_full_it = []
        self.overlap_front_it = []
        self.overlap_end_it = []
        self.overlap_ext_it = []
        self.gene_ids_it = set()

        self.unique_ref_overlaps_it = set()
        self.unique_ola_overlaps_it = set()

        self.unique_ref_proxpeak_it = set()
        self.unique_exp_proxpeak_it = set()

    def within(self, exp_bind: tuple, ref_bind: tuple):
        """Returns the overlap between exp_bind and ref_bind."""
        midpoint = (int)((int(ref_bind[0]) + int(ref_bind[1])) / 2)
        # Create the scope.
        overlap = range(
            max(int(exp_bind[0]), midpoint - self.scope),
            min(int(exp_bind[1]), midpoint + self.scope) + 1,
        )
        if len(overlap) == 0:
            return overlap, None, None
        else:
            # Adjust to -1000 to 1000 base pair range.
            lower = overlap[0]
            upper = overlap[-1]
            overlap = [x - midpoint for x in overlap]
            return overlap, lower, upper

    def coord_format(self, coord: tuple, chrom: str):
        """Returns the coordinate in the format 'chr:start-end'."""
        return f"{chrom}:{coord[0]}-{coord[1]}"

    def compare_chrom_bind_it(self):
        """Compares the binding sites for the given chromosome."""
        # Make all the Binding Comparisons
        for exp_bind in self.exp_binds_it:
            e_beg, e_end, e_bind = exp_bind
            overlaps = self.ref_binds_it[e_beg:e_end]
            # overlap, upper, lower = self.within(exp_bind, ref_bind)
            if len(overlaps) == 0:
                # No overlap in the binding sites... continue.
                continue
            else:
                for r_interval in overlaps:
                    r_beg, r_end, ref_bind = r_interval
                    overlap, upper, lower = self.within((e_beg, e_end), ref_bind)
                    if exp_bind[0] >= ref_bind[0] and exp_bind[1] <= ref_bind[1]:
                        # The EXP peak is fully contained by the REF peak
                        self.overlap_full_it.extend(overlap)
                        ot = "CRO"
                    elif exp_bind[1] > ref_bind[1] and exp_bind[0] <= ref_bind[1]:
                        # The EXP Peak overlaps the end of the REF peak.
                        self.overlap_end_it.extend(overlap)
                        ot = "ORE"
                    elif exp_bind[0] < ref_bind[0] and exp_bind[1] >= ref_bind[0]:
                        # The EXP Peak overlaps the front of the REF peak.
                        self.overlap_front_it.extend(overlap)
                        ot = "ORF"
                    else:
                        # The EXP Peak overlaps the REF peak externally.
                        self.overlap_ext_it.extend(overlap)
                        ot = "PXP"

                    if ot == "PXP":
                        self.unique_ref_proxpeak_it.add(ref_bind)
                        self.unique_exp_proxpeak_it.add(exp_bind)
                    else:
                        self.unique_ref_overlaps_it.add(ref_bind)
                        self.unique_ola_overlaps_it.add(exp_bind)
                    self.overlap_counts_it["Total"] += 1
                    self.overlap_counts_it[ot] += 1
                    self.overlaps_it.append(
                        (
                            self.chr,
                            ref_bind[0],
                            ref_bind[1],
                            self.coord_format(exp_bind, self.chr),
                            ot,
                        )
                    )

    def get_num_overlaps(self):
        """Returns the number of overlaps."""
        return self.num_overlaps

    def get_all_overlaps(self):
        """Returns all the overlaps."""
        return (
            self.overlaps,
            np.asarray(self.overlap_full),
            np.asarray(self.overlap_front),
            np.asarray(self.overlap_end),
            np.asarray(self.overlap_ext),
        )


class BindCompare:
    def __init__(self, ref_bed: Bed, exp_bed: Bed, scope: int):
        self.ref_bed = ref_bed
        self.exp_bed = exp_bed
        self.scope = scope
        self.experiments = {}
        self.font = {"fontname": "Sans Serif"}

    def compare_binds(self):
        """Compares the binding sites for each chromosome."""
        # Make all the Binding Comparisons
        for chromosome in self.ref_bed.get_chroms():
            if chromosome not in self.exp_bed.get_chroms():
                continue
            self.experiments[chromosome] = Chromosome(
                chromosome,
                self.scope,
                self.ref_bed.get_chrom_it(chromosome),
                self.exp_bed.get_chrom_it(chromosome),
            )
            # begin = time.time()
            self.experiments[chromosome].compare_chrom_bind_it()
            # end = time.time()
            # print("interval-tree took {}".format(end - begin))

    def get_experiment(self, chromosome: str):
        """Returns the experiment for the given chromosome."""
        return self.experiments[chromosome]

    def get_experiments_overlaps_it(self, chromosomes: list):
        """Returns the experiments for the given chromosomes."""
        olaps = []
        full_o = []
        front_o = []
        end_o = []
        ext_o = []
        overlap_counts = {"Total": 0, "CRO": 0, "ORE": 0, "ORF": 0, "PXP": 0}
        unique_ref_overlaps = set()
        unique_ola_overlaps = set()
        unique_ref_proxpeak = set()
        unique_exp_proxpeak = set()
        for chromosome in chromosomes:
            if chromosome not in self.experiments:
                continue
            olaps.extend(self.experiments[chromosome].overlaps_it)
            full_o.extend(self.experiments[chromosome].overlap_full_it)
            front_o.extend(self.experiments[chromosome].overlap_front_it)
            end_o.extend(self.experiments[chromosome].overlap_end_it)
            ext_o.extend(self.experiments[chromosome].overlap_ext_it)
            for key in self.experiments[chromosome].overlap_counts_it:
                overlap_counts[key] += self.experiments[chromosome].overlap_counts_it[
                    key
                ]
            unique_ref_overlaps.update(
                self.experiments[chromosome].unique_ref_overlaps_it
            )
            unique_ola_overlaps.update(
                self.experiments[chromosome].unique_ola_overlaps_it
            )
            unique_ref_proxpeak.update(
                self.experiments[chromosome].unique_ref_proxpeak_it
            )
            unique_exp_proxpeak.update(
                self.experiments[chromosome].unique_exp_proxpeak_it
            )

        return {
            "overlaps": olaps,
            "full": np.asarray(full_o),
            "front": np.asarray(front_o),
            "end": np.asarray(end_o),
            "ext": np.asarray(ext_o),
            "overlap_counts": overlap_counts,
            "unique_ref_overlaps": unique_ref_overlaps,
            "unique_ola_overlaps": unique_ola_overlaps,
            "unique_ref_proxpeak": unique_ref_proxpeak,
            "unique_exp_proxpeak": unique_exp_proxpeak,
        }

    def scatter_overlap_freq(self, bc_dict: dict, filepath: str):
        fig = plt.figure()
        ax1 = fig.add_subplot(111)

        x, y = np.unique(
            bc_dict["full"], return_counts=True
        )  # counting occurrence of each loan
        ax1.scatter(x, y, s=10, c="m", marker="s", label="Complete Peak Overlap (CRO)")

        x, y = np.unique(
            bc_dict["front"], return_counts=True
        )  # counting occurrence of each loan
        ax1.scatter(x, y, s=5, c="r", marker="o", label="Overlap Ref Front (ORF)")

        x, y = np.unique(
            bc_dict["end"], return_counts=True
        )  # counting occurrence of each loan
        ax1.scatter(x, y, s=5, c="b", marker="o", label="Overlaps Ref End (ORE)")

        x, y = np.unique(
            bc_dict["ext"], return_counts=True
        )  # counting occurrence of each loan
        ax1.scatter(x, y, s=5, c="y", marker="o", label="Proximal Peaks (PXP)")

        ax1.set_xlabel("Overlap of Binding Sites", **self.font)
        ax1.set_ylabel("Frequency", **self.font)

        plt.legend(loc="upper left")
        plt.title(
            f"Frequency of Binding Overlaps Over {2 * self.scope} Base Pair Range"
        )

        plt.savefig(filepath + "_overlaps.png", dpi=300)
        plt.close()

    def overlap_distribution_barplot(self, bc_dict: dict, filepath: str):
        """Takes the overlap counts and generates a single bar split by overlap category."""
        # Sample data generation
        overlap_types = [
            "Full Overlap (CRO)",
            "Front Overlap (ORF)",
            "End Overlap (ORE)",
            "Proximal Peak (PXP)",
        ]
        counts = bc_dict["overlap_counts"]
        values = [counts["CRO"], counts["ORF"], counts["ORE"], counts["PXP"]]
        total = sum(values)

        if total == 0:
            print("Overlap Distribution Plot Terminated... no counts found.")
            return

        values = [val / total for val in values]

        # Your code for plotting
        fig, ax = plt.subplots(figsize=(10, 2))

        left = 0
        colors = ["purple", "red", "blue", "yellow"]
        for i, value in enumerate(values):
            ax.barh(
                0,
                value,
                left=left,
                label=overlap_types[i],
                edgecolor="black",
                linewidth=1,
                color=colors[i],
            )
            left += value

        ax.set_xlim(0, 1)
        ax.set_ylim(-0.5, 0.5)
        ax.set_xlabel("Distribution of Overlap Type Frequencies")

        # Hide y-axis labels
        ax.get_yaxis().set_visible(False)

        # Set a title
        ax.set_title("Stacked Bar Plot of Overlap Categories")

        # Add a custom legend item for the total without a color
        legend_total = mpatches.Patch(label=f"Total: {total:.0f}", color="none")
        handles, labels = ax.get_legend_handles_labels()
        handles.append(legend_total)

        # Create the legend
        ax.legend(handles=handles, loc="center left", bbox_to_anchor=(1, 0.5))

        plt.subplots_adjust(
            left=0.05, right=0.85, top=0.9, bottom=0.1
        )  # Adjust the layout

        plt.savefig(filepath + "_bardist.png", bbox_inches="tight", dpi=300)
        plt.close()

    def overlap_distribution_piechart(self, bc_dict: dict, filepath: str):
        """Takes the overlap counts and generates a pie chart split by overlap category."""
        categories = [
            "Complete Overlap",
            "Overlap Ref. Front",
            "Overlap Ref. End",
            "Proximal Peak",
        ]

        counts = [
            bc_dict["overlap_counts"]["CRO"],
            bc_dict["overlap_counts"]["ORF"],
            bc_dict["overlap_counts"]["ORE"],
            bc_dict["overlap_counts"]["PXP"],
        ]
        colors = ["purple", "red", "blue", "yellow"]
        total = sum(counts)
        fig1, ax1 = plt.subplots()

        def pie_fmt(x):
            return "{:.0f}".format((total) * x / 100)

        ax1.pie(
            counts, labels=categories, autopct=pie_fmt, startangle=90, colors=colors
        )
        ax1.axis("equal")

        plt.title(f"Categorization of {total} Found Overlaps", **self.font)
        plt.savefig(filepath + "_pie.png", dpi=300)
        plt.close()

    def overlap_bar_totals(self, bc_dict: dict, filepath: str):
        categories = [
            "Exp. Binding\nPeaks",
            "Unique Overlaps",
            "Total No.\nof Overlaps",
            "Unique\nProx. Peaks",
            "Total No. of\nProx. Peaks",
            "Reference\nPeaks Identified",
        ]

        overlap_total = (
            bc_dict["overlap_counts"]["CRO"]
            + bc_dict["overlap_counts"]["ORF"]
            + bc_dict["overlap_counts"]["ORE"]
        )

        complete = set()
        complete.update(bc_dict["unique_ref_overlaps"])
        complete.update(bc_dict["unique_ref_proxpeak"])

        vals = [
            self.exp_bed.num_peaks,
            len(bc_dict["unique_ola_overlaps"]),
            overlap_total,
            len(bc_dict["unique_exp_proxpeak"]),
            bc_dict["overlap_counts"]["PXP"],
            len(complete),
        ]
        bars = plt.bar(
            categories,
            vals,
            color=[
                "steelblue",
                "sandybrown",
                "lightgrey",
                "peru",
                "slategrey",
                "cornflowerblue",
            ],
        )

        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + 0.25, yval + 0.35, yval, wrap=True)

        plt.xticks(fontsize=7)

        plt.ylabel("")
        plt.xlabel("")
        plt.title("Total Number of Overlaps and Binding Peaks in Overlayed Bed")
        plt.savefig(filepath + "_barsummary.png", dpi=300)
        plt.close()

    def generate_csv(
        self, bc_dict: dict, filepath: str, outdir: str, name: str, gtf: GTF = None
    ):
        df = pd.DataFrame(
            bc_dict["overlaps"],
            columns=[
                "Chrom",
                "Begin Ref Site",
                "End Ref Site",
                "Experimental_Peaks",
                "Overlay_Type",
            ],
        )

        group_columns = ["Chrom", "Begin Ref Site", "End Ref Site", "Overlay_Type"]
        df = df.groupby(group_columns, as_index=False).agg(
            Experimental_Peaks=("Experimental_Peaks", ",".join),
            Occurrences=("Overlay_Type", "count"),
        )

        df["Occurrences"] = df["Occurrences"].apply(lambda x: str(x) + ",")
        df["Overlay_Type"] = df["Overlay_Type"].apply(lambda x: str(x) + ",")
        columns = [
            "Chrom",
            "Begin Ref Site",
            "End Ref Site",
            "Experimental_Peaks",
            "Overlay_Type",
            "Occurrences",
        ]
        df = (
            df.groupby(["Chrom", "Begin Ref Site", "End Ref Site"])
            .agg(
                Overlay_Type=("Overlay_Type", "sum"),
                Occurrences=("Occurrences", "sum"),
                Experimental_Peaks=("Experimental_Peaks", "sum"),
            )
            .reset_index()
        )

        df["Occurrences"] = df["Occurrences"].apply(lambda x: x[:-1])
        df["Overlay_Type"] = df["Overlay_Type"].apply(lambda x: x[:-1])
        if gtf is not None:
            all_genes = set()
            df["GeneIDs"] = df.apply(
                lambda row: gtf.find_fbgn(
                    row["Chrom"],
                    row["Begin Ref Site"],
                    row["End Ref Site"],
                    self.scope,
                    all_genes,
                ),
                axis=1,
            )
            bc_dict["all_genes"] = all_genes
            gene_arr = np.array(list(all_genes))
            gdf = pd.DataFrame(gene_arr.reshape(len(gene_arr), -1), columns=["Gene ID"])
            outpath = "/tmp/gene_list.csv"
            gdf.to_csv(outpath, index=False)
        df.to_csv(filepath + "_overlaps.csv", index=False)

        ## Create Separate CSVs for Each Overlap Type
        path = os.path.join(outdir, "CategorizedCSVs")
        os.makedirs(path, exist_ok=True)
        for type_o in ["CRO", "ORF", "ORE", "PXP"]:
            result = df[df["Overlay_Type"].apply(lambda x: type_o in x)]
            path2 = os.path.join(path, type_o + f"_{name}_overlaps.csv")
            result.to_csv(path2, index=False)

    def plot_perchrom_ref_peak(self, filepath: str):
        """The same plot as plot_average_ref_peak but creates N subplots for each chromosome on one panel.
        There will be at most 3 plots on each row and there will be however many rows needed to fit all the chromosomes.
        """
        # Determine Number of Chromosomes
        chroms = sorted(
            list(self.experiments.keys()),
            key=lambda x: (x.isdigit(), int(x) if x.isdigit() else x),
        )
        Nchroms = len(chroms)

        num_rows = Nchroms // 3
        num_rows = num_rows if num_rows * 3 == Nchroms else num_rows + 1
        fig, axs = plt.subplots(
            num_rows, 3, sharex=True, sharey="all", figsize=(13, Nchroms // 2 + 3)
        )
        if Nchroms == 1:
            axs = [axs]
        else:
            axs = axs.flatten()

        for i in range(Nchroms):
            chrom, ax = chroms[i], axs[i]
            chrom_t: Chromosome = self.experiments[chrom]
            # num_peaks = len(chrom_t.ref_binds)
            num_peaks = len(chrom_t.ref_binds_it)
            x = np.arange(-self.scope, self.scope + 1, 1)
            y = np.zeros(2 * self.scope + 1)

            chrom_dict = self.ref_bed.get_loci([chrom])
            for binding_site in chrom_dict[chrom]:
                midpoint = (binding_site[1] + binding_site[0]) / 2
                start_index = int(binding_site[0] - midpoint + self.scope)
                end_index = int(binding_site[1] - midpoint + self.scope)
                y[start_index:end_index] += 1
            y = y / num_peaks
            ax.plot(x, y, label="Average Ref. Peak", c="k", alpha=0.7)
            ax2 = ax.twinx()

            data_arrays = [
                chrom_t.overlap_full_it,
                chrom_t.overlap_front_it,
                chrom_t.overlap_end_it,
                chrom_t.overlap_ext_it,
            ]
            colors = ["m", "r", "b", "y"]

            for arr, color in zip(data_arrays, colors):
                x, z = np.unique(arr, return_counts=True)
                # if len(arr) != 0:
                # z = moving_average(z, 20)
                ax2.plot(x, z, c=color, alpha=0.7)

            ax.set_title(f"Chromosome {chrom}")

        # Create a custom legend outside the function
        custom_legend = [
            plt.Line2D([0], [0], color="k", lw=2, label="Average Ref. Peak"),
            plt.Line2D([0], [0], color="m", lw=2, label="CRO"),
            plt.Line2D([0], [0], color="r", lw=2, label="ORF"),
            plt.Line2D([0], [0], color="b", lw=2, label="ORE"),
            plt.Line2D([0], [0], color="y", lw=2, label="PXP"),
        ]
        fig.legend(
            handles=custom_legend,
            loc="lower center",
            fancybox=True,
            shadow=True,
            fontsize=11,
            ncol=5,
        ).set_bbox_to_anchor((0.5, 0.0))
        # Set x-axis label at the center of the left
        fig.text(
            0.5,
            0.07,
            "Distance from Reference Peak Midpoint",
            ha="center",
            va="center",
            fontsize=12,
        )

        # Set y-axis label at the center of the bottom
        fig.text(
            0.02,
            0.5,
            "Frequency of Reference Peaks",
            ha="center",
            va="center",
            rotation="vertical",
            fontsize=12,
        )

        fig.text(
            0.98,
            0.5,
            "Overlaps Counts",
            ha="center",
            va="center",
            rotation="vertical",
            fontsize=12,
        )

        # plt.subplots_adjust(bottom=0.4)
        fig.suptitle(
            "Per Chromosome Counts of Binding Overlaps Across Reference Binding Peak",
        )
        fig.tight_layout()
        fig.subplots_adjust(left=0.08, right=0.92, top=0.9, bottom=0.12)
        plt.savefig(filepath + "_chrom_ref_freq.png", dpi=300)
        plt.close()

    def plot_average_ref_peak(self, bc_dict: dict, filepath: str):
        """With the x-axis being -scope to scope, plot all the frequency of the reference peaks over this domain.
        normalize the y-axis and make the plot a line curve."""

        num_peaks = self.ref_bed.num_peaks
        x = np.arange(-self.scope, self.scope + 1, 1)
        y = np.zeros(2 * self.scope + 1)
        chrom_dict = self.ref_bed.get_loci(self.ref_bed.get_chroms())
        for chr in chrom_dict:
            for binding_site in chrom_dict[chr]:
                midpoint = (binding_site[1] + binding_site[0]) / 2
                start_index = int(binding_site[0] - midpoint + self.scope)
                end_index = int(binding_site[1] - midpoint + self.scope)
                y[start_index:end_index] += 1
        y = y / num_peaks
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ln1 = ax1.plot(x, y, label="Average Ref. Peak", c="k", alpha=0.7)

        ax2 = ax1.twinx()
        x, z = np.unique(bc_dict["full"], return_counts=True)
        # if len(z) != 0:
        #     z = moving_average(z, 20)
        ln2 = ax2.plot(x, z, c="m", label="Complete Ref Overlap (CRO)", alpha=0.7)

        x, z = np.unique(bc_dict["front"], return_counts=True)
        # if len(z) != 0:
        #     z = moving_average(z, 20)
        ln3 = ax2.plot(x, z, c="r", label="Overlap Ref Front (ORF)", alpha=0.7)

        x, z = np.unique(bc_dict["end"], return_counts=True)
        # if len(z) != 0:
        #     z = moving_average(z, 20)
        ln4 = ax2.plot(x, z, c="b", label="Overlaps Ref End (ORE)", alpha=0.7)

        x, z = np.unique(bc_dict["ext"], return_counts=True)
        # if len(z) != 0:
        #     z = moving_average(z, 20)
        ln5 = ax2.plot(x, z, c="y", label="Proximal Peaks (PXP)", alpha=0.7)

        lns = ln1 + ln2 + ln3 + ln4 + ln5
        labs = [l.get_label() for l in lns]

        plt.title("Counts of Binding Overlaps Over Average Reference Peak Profile")
        plt.xlabel("Distance from Reference Peak Midpoint")
        ax1.set_ylabel("Frequency of Reference Peaks")
        ax2.set_ylabel("Overlaps Counts")

        plt.legend(
            lns,
            labs,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.1),
            fancybox=True,
            shadow=True,
            fontsize=8,
            ncol=3,
        )
        plt.tight_layout()
        plt.savefig(filepath + "_ref_freq.png", dpi=300)
        plt.close()

    def generate_summary(self, bc_dict: dict, filepath: str, gtf: GTF = None):
        with open(filepath + "_summary.txt", "a") as summary:
            summary.write(
                f"Total Number of Reference Peaks: {self.ref_bed.num_peaks}\n"
            )
            summary.write(
                f"Total Number of Experimental Peaks: {self.exp_bed.num_peaks}\n"
            )
            summary.write(
                f"Total Number of Unique Overlaps: {len(bc_dict['unique_ola_overlaps'])}\n"
            )
            summary.write(
                f"Total Number of Overlaps: {bc_dict['overlap_counts']['Total']}\n"
            )
            summary.write(
                f"Total Number of Reference Peaks Overlapped: {len(bc_dict['unique_ref_overlaps'])}\n"
            )
            if gtf:
                summary.write(
                    f"Total Number of Genes Overlapped: {len(bc_dict['all_genes'])}\n"
                )
                all_genes = ""
                for gene in list(bc_dict["all_genes"]):
                    all_genes += gene + " "
                summary.write(f"List of All Genes:\n{all_genes}")

    def generate_all(self, bc_dict: dict, outpath: str, name: str, gtf: GTF = None):
        """Generates all the visualizations and csv files."""
        # self.scatter_overlap_freq(bc_dict, filepath)
        filepath = outpath + name
        self.overlap_distribution_barplot(bc_dict, filepath)
        self.overlap_distribution_piechart(bc_dict, filepath)
        self.plot_average_ref_peak(bc_dict, filepath)
        self.plot_perchrom_ref_peak(filepath)
        self.overlap_bar_totals(bc_dict, filepath)
        self.generate_csv(bc_dict, filepath, outpath, name, gtf)
        self.generate_summary(bc_dict, filepath, gtf)
