import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
from merge_class import Bed, GTF


class Chromosome:
    def __init__(self, name: str, ref_binds: list, exp_binds: list, scope: int):
        self.chr = name
        self.ref_binds = ref_binds
        self.exp_binds = exp_binds
        self.scope = scope
        self.overlap_counts = {"Total": 0, "OF": 0, "OE": 0, "OB": 0, "OX": 0}
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
                        ot = "OF"
                    elif exp_bind[1] > ref_bind[1] and exp_bind[0] <= ref_bind[1]:
                        # The EXP Peak overlaps the end of the REF peak.
                        self.overlap_end.extend(overlap)
                        ot = "OE"
                    elif exp_bind[0] < ref_bind[0] and exp_bind[1] >= ref_bind[0]:
                        # The EXP Peak overlaps the front of the REF peak.
                        self.overlap_front.extend(overlap)
                        ot = "OB"
                    else:
                        # The EXP Peak overlaps the REF peak externally.
                        self.overlap_ext.extend(overlap)
                        ot = "OX"
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
        overlap_counts = {"Total": 0, "OF": 0, "OE": 0, "OB": 0, "OX": 0}
        unique_ref_overlaps = set()
        unique_ola_overlaps = set()
        for chromosome in chromosomes:
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
        ax1.scatter(x, y, s=10, c="m", marker="s", label="Complete Peak Overlap")

        x, y = np.unique(
            bc_dict["front"], return_counts=True
        )  # counting occurrence of each loan
        ax1.scatter(x, y, s=5, c="r", marker="o", label="Overlap Ref Front")

        x, y = np.unique(
            bc_dict["end"], return_counts=True
        )  # counting occurrence of each loan
        ax1.scatter(x, y, s=5, c="b", marker="o", label="Overlaps Ref End")

        x, y = np.unique(
            bc_dict["ext"], return_counts=True
        )  # counting occurrence of each loan
        ax1.scatter(x, y, s=5, c="y", marker="o", label="External Overlaps")

        ax1.set_xlabel("Overlap of Binding Sites", **self.font)
        ax1.set_ylabel("Frequency", **self.font)

        plt.legend(loc="upper left")
        plt.title(
            f"Frequency of Binding Overlaps Over {2 * self.scope} Base Pair Range"
        )

        plt.savefig(filepath + "_overlaps.png")
        plt.close()

    def overlap_distribution_barplot(self, bc_dict: dict, filepath: str):
        """Takes the overlap counts and generates a single bar split by overlap category."""
        # Sample data generation
        overlap_types = ["Full Overlap", "Front Overlap", "End Overlap", "Ext. Overlap"]
        counts = bc_dict["overlap_counts"]
        values = [counts["OF"], counts["OB"], counts["OE"], counts["OX"]]
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

        plt.savefig(filepath + "_bardist.png", bbox_inches="tight")
        plt.close()

    def overlap_distribution_piechart(self, bc_dict: dict, filepath: str):
        """Takes the overlap counts and generates a pie chart split by overlap category."""
        categories = [
            "Complete Overlap",
            "Partial Overlap Front",
            "Partial Overlap End",
            "External Overlap",
        ]

        abbr = ["OF", "OE", "OB", "OX"]
        counts = [
            bc_dict["overlap_counts"]["OF"],
            bc_dict["overlap_counts"]["OE"],
            bc_dict["overlap_counts"]["OB"],
            bc_dict["overlap_counts"]["OX"],
        ]
        total = sum(counts)
        fig1, ax1 = plt.subplots()

        def pie_fmt(x):
            return "{:.0f}".format((total) * x / 100)

        ax1.pie(counts, labels=categories, autopct=pie_fmt, startangle=90)
        ax1.axis("equal")

        plt.title(f"Categorization of {total} Found Overlaps", **self.font)
        plt.savefig(filepath + "_pie.png")
        plt.close()

    def overlap_bar_totals(self, bc_dict: dict, filepath: str):
        categories = [
            "Total\nBinding Peaks",
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
        plt.savefig(filepath + "_barsummary.png")
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
            outpath = "gene_list.csv"
            gdf.to_csv(outpath, index=False)
        df.to_csv(filepath + "_overlaps.csv", index=False)

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
        self.scatter_overlap_freq(bc_dict, filepath)
        self.overlap_distribution_barplot(bc_dict, filepath)
        self.overlap_distribution_piechart(bc_dict, filepath)
        self.overlap_bar_totals(bc_dict, filepath)
        self.generate_csv(bc_dict, filepath, gtf)
        self.generate_summary(bc_dict, filepath)
