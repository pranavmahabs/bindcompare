import argparse
import sys
import os

from .bindapp import merge
from .bindapp import downstream


def main():
    parser = argparse.ArgumentParser(
        description="bindcompare: Reveal and analyze co-regulatory sites across two protein-binding datasets."
    )

    parser.add_argument(
        "-r",
        "--ref",
        required=True,
        help="The file path for your reference BED file. Should be the DNA path if comparing DNA and RNA.",
    )
    parser.add_argument(
        "-e", "--exp", required=True, help="The file path for your overlayed BED file."
    )
    parser.add_argument(
        "-s",
        "--scope",
        required=True,
        help="Nucleotides upstream and downstream from the ref peak's center that BC will search for peaks.",
    )
    parser.add_argument(
        "-n", "--name", required=True, help="Str name for this experiment."
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output Directory. All outputs for bindcompare will be found here. Will be created if DNE.",
    )
    parser.add_argument(
        "-g",
        "--gtf",
        default="None",
        help="Gene GTF file in proper format.",
    )
    parser.add_argument(
        "-f",
        "--fasta",
        default="None",
        help="Genome file corresponding to your BED Files.",
    )

    args = parser.parse_args()

    summary_output = os.path.join(args.out, f"{args.name}_summary.txt")

    if not os.path.exists(args.out):
        os.makedirs(args.out, exist_ok=True)

    with open(summary_output, "w") as summary_file:
        merge.main(args.ref, args.exp, args.scope, args.name, f"{args.out}/", args.gtf)

    with open(summary_output, "a") as summary_file:
        downstream.downstream(
            os.path.join(args.out, f"{args.name}_overlaps.csv"), args.fasta, args.out
        )

    print("Completed BindCompare! Time Stamp:")
    os.system("date")
