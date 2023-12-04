import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from merge_class import Bed, GTF


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
    def __init__(self, name: str, ref_binds: list, exp_binds: list, scope: int):
        self.chr = name
        self.ref_binds = ref_binds
        self.exp_binds = exp_binds
        self.scope = scope
        self.overlap_counts = {"Total": 0, "CRO": 0, "ORE": 0, "ORF": 0, "PXP": 0}
        self.overlaps = []
        self.overlap_full = []
        self.overlap_front = []
        self.overlap_end = []
        self.overlap_ext = []
        self.gene_ids = set()

        self.unique_ref_overlaps = set()
        self.unique_ola_overlaps = set()

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

    def compare_chrom_bind(self):
        """Compares the binding sites for the given chromosome."""
        # Make all the Binding Comparisons
        for ref_bind in self.ref_binds:
            for exp_bind in self.exp_binds:
                overlap, upper, lower = self.within(exp_bind, ref_bind)
                if len(overlap) == 0:
                    # No overlap in the binding sites... continue.
                    continue
                else:
                    if exp_bind[0] >= ref_bind[0] and exp_bind[1] <= ref_bind[1]:
                        # The EXP peak is fully contained by the REF peak
                        self.overlap_full.extend(overlap)
                        ot = "CRO"
                    elif exp_bind[1] > ref_bind[1] and exp_bind[0] <= ref_bind[1]:
                        # The EXP Peak overlaps the end of the REF peak.
                        self.overlap_end.extend(overlap)
                        ot = "ORE"
                    elif exp_bind[0] < ref_bind[0] and exp_bind[1] >= ref_bind[0]:
                        # The EXP Peak overlaps the front of the REF peak.
                        self.overlap_front.extend(overlap)
                        ot = "ORF"
                    else:
                        # The EXP Peak overlaps the REF peak externally.
                        self.overlap_ext.extend(overlap)
                        ot = "PXP"
                    self.unique_ref_overlaps.add(ref_bind)
                    self.unique_ola_overlaps.add(exp_bind)
                    self.overlap_counts["Total"] += 1
                    self.overlap_counts[ot] += 1
                    self.overlaps.append(
                        (
                            self.chr,
                            ref_bind[0],
                            ref_bind[1],
                            self.coord_format(exp_bind, self.chr),
                            # self.coord_format((lower, upper), self.chr),
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
                self.ref_bed.get_chrom(chromosome),
                self.exp_bed.get_chrom(chromosome),
                self.scope,
            )
            self.experiments[chromosome].compare_chrom_bind()

    def get_experiment(self, chromosome: str):
        """Returns the experiment for the given chromosome."""
        return self.experiments[chromosome]

    def get_experiments_overlaps(self, chromosomes: list):
        """Returns the experiments for the given chromosomes."""
        olaps = []
        full_o = []
        front_o = []
        end_o = []
        ext_o = []
        overlap_counts = {"Total": 0, "CRO": 0, "ORE": 0, "ORF": 0, "PXP": 0}
        unique_ref_overlaps = set()
        unique_ola_overlaps = set()
        for chromosome in chromosomes:
            if chromosome not in self.experiments:
                continue
            olaps.extend(self.experiments[chromosome].overlaps)
            full_o.extend(self.experiments[chromosome].overlap_full)
            front_o.extend(self.experiments[chromosome].overlap_front)
            end_o.extend(self.experiments[chromosome].overlap_end)
            ext_o.extend(self.experiments[chromosome].overlap_ext)
            for key in self.experiments[chromosome].overlap_counts:
                overlap_counts[key] += self.experiments[chromosome].overlap_counts[key]
            unique_ref_overlaps.update(self.experiments[chromosome].unique_ref_overlaps)
            unique_ola_overlaps.update(self.experiments[chromosome].unique_ola_overlaps)

        return {
            "overlaps": olaps,
            "full": np.asarray(full_o),
            "front": np.asarray(front_o),
            "end": np.asarray(end_o),
            "ext": np.asarray(ext_o),
            "overlap_counts": overlap_counts,
            "unique_ref_overlaps": unique_ref_overlaps,
            "unique_ola_overlaps": unique_ola_overlaps,
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
        total = sum(counts)
        fig1, ax1 = plt.subplots()

        def pie_fmt(x):
            return "{:.0f}".format((total) * x / 100)

        ax1.pie(counts, labels=categories, autopct=pie_fmt, startangle=90)
        ax1.axis("equal")

        plt.title(f"Categorization of {total} Found Overlaps", **self.font)
        plt.savefig(filepath + "_pie.png", dpi=300)
        plt.close()

    def overlap_bar_totals(self, bc_dict: dict, filepath: str):
        categories = [
            "Total Exp.\nBinding Peaks",
            "Unique Overlaps",
            "Total Number\nof Overlaps",
            "Reference\nPeaks Overlapped",
        ]

        vals = [
            self.exp_bed.num_peaks,
            len(bc_dict["unique_ola_overlaps"]),
            bc_dict["overlap_counts"]["Total"],
            len(bc_dict["unique_ref_overlaps"]),
        ]
        bars = plt.bar(
            categories,
            vals,
            color=["mediumblue", "sandybrown", "lightgrey", "cornflowerblue"],
        )

        for bar in bars:
            yval = bar.get_height()
            plt.text(bar.get_x() + 0.3, yval + 0.3, yval, wrap=True)

        plt.ylabel("")
        plt.xlabel("")
        plt.title("Total Number of Overlaps and Binding Peaks in Overlayed Bed")
        plt.savefig(filepath + "_barsummary.png", dpi=300)
        plt.close()

    def generate_csv(self, bc_dict: dict, filepath: str, gtf: GTF = None):
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

    def plot_perchrom_ref_peak(self, filepath: str):
        """The same plot as plot_average_ref_peak but creates N subplots for each chromosome on one panel.
        There will be at most 3 plots on each row and there will be however many rows needed to fit all the chromosomes.
        """
        # Determine Number of Chromosomes
        chroms = list(self.experiments.keys())
        Nchroms = len(chroms)

        fig, axs = plt.subplots(
            Nchroms, 1, sharex=True, sharey="all", figsize=(4, Nchroms + 3)
        )
        if Nchroms == 1:
            axs = [axs]
        else:
            axs = axs.flatten()

        for i in range(Nchroms):
            chrom, ax = chroms[i], axs[i]
            chrom_t: Chromosome = self.experiments[chrom]
            num_peaks = len(chrom_t.ref_binds)
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
                chrom_t.overlap_full,
                chrom_t.overlap_front,
                chrom_t.overlap_end,
                chrom_t.overlap_ext,
            ]
            colors = ["m", "r", "b", "y"]

            for arr, color in zip(data_arrays, colors):
                x, z = np.unique(arr, return_counts=True)
                if len(arr) != 0:
                    z = moving_average(z, 20)
                ax2.plot(x, z, c=color, alpha=0.7)

            if (i + 1) % Nchroms == 0:
                ax.set_xlabel("Distance from Reference Peak Midpoint")
            if i == int(Nchroms / 2):
                ax.set_ylabel("Frequency of Reference Peaks")
                ax2.set_ylabel("Overlaps Counts")
            ax.set_title(f"Chromosome {chrom}")

        # Create a custom legend outside the function
        custom_legend = [
            plt.Line2D([0], [0], color="k", lw=2, label="Average Ref. Peak"),
            plt.Line2D([0], [0], color="m", lw=2, label="CRO"),
            plt.Line2D([0], [0], color="r", lw=2, label="ORF"),
            plt.Line2D([0], [0], color="b", lw=2, label="ORE"),
            plt.Line2D([0], [0], color="y", lw=2, label="PXP"),
        ]
        plt.legend(
            handles=custom_legend,
            loc="upper center",
            bbox_to_anchor=(0.55, -0.5),
            fancybox=True,
            shadow=True,
            fontsize=7,
            ncol=3,
        )
        plt.subplots_adjust(bottom=0.4)
        fig.suptitle(
            "Per Chromosome Counts of Binding Overlaps\n Across Reference Binding Peak"
        )
        plt.tight_layout(rect=[0, 0, 1, 1])
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
        if len(z) != 0:
            z = moving_average(z, 20)
        ln2 = ax2.plot(x, z, c="m", label="Complete Ref Overlap (CRO)", alpha=0.7)

        x, z = np.unique(bc_dict["front"], return_counts=True)
        if len(z) != 0:
            z = moving_average(z, 20)
        ln3 = ax2.plot(x, z, c="r", label="Overlap Ref Front (ORF)", alpha=0.7)

        x, z = np.unique(bc_dict["end"], return_counts=True)
        if len(z) != 0:
            z = moving_average(z, 20)
        ln4 = ax2.plot(x, z, c="b", label="Overlaps Ref End (ORE)", alpha=0.7)

        x, z = np.unique(bc_dict["ext"], return_counts=True)
        if len(z) != 0:
            z = moving_average(z, 20)
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
        with open(filepath + "_summary.txt", "w") as summary:
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
            if gtf is not None:
                summary.write(
                    f"Total Number of Genes Overlapped: {len(bc_dict['all_genes'])}\n"
                )

    def generate_all(self, bc_dict: dict, filepath: str, gtf: GTF = None):
        """Generates all the visualizations and csv files."""
        # self.scatter_overlap_freq(bc_dict, filepath)
        self.overlap_distribution_barplot(bc_dict, filepath)
        self.overlap_distribution_piechart(bc_dict, filepath)
        self.plot_average_ref_peak(bc_dict, filepath)
        self.plot_perchrom_ref_peak(filepath)
        self.overlap_bar_totals(bc_dict, filepath)
        self.generate_csv(bc_dict, filepath, gtf)
        self.generate_summary(bc_dict, filepath)
