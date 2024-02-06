import os
import re
import sys
import argparse
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from collections import defaultdict


def extract_genes_from_summary(file_path):
    genes = set()
    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("List of All Genes:"):
                genes.update(re.findall(r"\b\w+\b", next(file)))
    return genes


def create_summary_file(genes_folder1, genes_folder2, prefix1, prefix2):
    # Create Venn diagram
    venn = venn2([genes_folder1, genes_folder2], set_labels=(prefix1, prefix2))

    # Prepare genes in each section
    venn_labels = defaultdict(list)
    for idx, label in venn_labels:
        if label != "":
            genes = venn.get_label_by_id(idx).get_text().split("\n")
            venn_labels[label].extend(genes)

    plt.title(f"Comparing Co-Regulatory Gene Regions for \n{prefix1} and {prefix2}")
    plt.tight_layout()

    venn_file_name = f"{prefix1}_v_{prefix2}_venn.png"
    # Show Venn diagram
    plt.savefig(venn_file_name, dpi=300)

    # Calculate Jaccard similarity score
    intersection = len(genes_folder1.intersection(genes_folder2))
    union = len(genes_folder1.union(genes_folder2))
    jaccard_similarity = intersection / union

    # Categorize genes
    exclusive_to_folder1 = list(genes_folder1.difference(genes_folder2))
    common_genes = genes_folder1.intersection(genes_folder2)
    exclusive_to_folder2 = genes_folder2.difference(genes_folder1)

    # Create comparison summary file
    summary_file_name = f"{prefix1}_v_{prefix2}_summary.txt"
    with open(summary_file_name, "w") as summary_file:
        summary_file.write(f"Jaccard Similarity: {jaccard_similarity:.4f}\n")
        summary_file.write(
            f"Genes Exclusive to {prefix1}: {' '.join(exclusive_to_folder1)}\n"
        )
        summary_file.write(f"Genes in Both Samples: {' '.join(common_genes)}\n")
        summary_file.write(
            f"Genes Exclusive to {prefix2}: {' '.join(exclusive_to_folder2)}\n"
        )

    return summary_file_name


def verify_summary_file(folder_path):
    summary_files = []
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith("_summary.txt"):
                summary_file_path = os.path.join(root, file)
                with open(summary_file_path, "r") as summary_file:
                    if any(
                        line.startswith("List of All Genes:") for line in summary_file
                    ):
                        summary_files.append(summary_file_path)
    return summary_files


def main():
    parser = argparse.ArgumentParser(
        description="comparexp: Compare two bindcompare experiments."
    )

    parser.add_argument(
        "-a",
        "--bindpath_1",
        help="Path to the first bindcompare Output Directory.",
        required=True,
    )
    parser.add_argument(
        "-b",
        "--bindpath_2",
        help="Path to the second bindcompare Output Directory.",
        required=True,
    )

    args = parser.parse_args()

    folder1_path = args.bindpath_1
    folder2_path = args.bindpath_2

    # Verify summary files in both folders
    summary_files_folder1 = verify_summary_file(folder1_path)
    summary_files_folder2 = verify_summary_file(folder2_path)

    if summary_files_folder1 and summary_files_folder2:
        # Extract gene lists from summary files
        if len(summary_files_folder1) > 1 or len(summary_files_folder2) > 1:
            print("Error: multiple summary files detected in one of the folders.")
            sys.exit(1)

        genes_folder1 = extract_genes_from_summary(summary_files_folder1[0])
        genes_folder2 = extract_genes_from_summary(summary_files_folder2[0])

        # Get prefixes for summary files
        folder1_prefix = os.path.basename(summary_files_folder1[0]).split(
            "_summary.txt"
        )[0]
        folder2_prefix = os.path.basename(summary_files_folder2[0]).split(
            "_summary.txt"
        )[0]

        # Create Venn diagram and comparison summary file
        summary_file_name = create_summary_file(
            genes_folder1, genes_folder2, folder1_prefix, folder2_prefix
        )
        print(f"Comparison summary file '{summary_file_name}' created successfully.")
    else:
        print("No suitable summary files found in both folders.")
