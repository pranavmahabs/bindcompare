import argparse
from datetime import datetime
import os

from .bindapp import merge
from .bindapp import downstream


def is_valid_file(parser, arg, name):
    if not os.path.isfile(arg):
        parser.error(f"The {name} {arg} does not exist!")
    else:
        return arg


def is_valid_directory(parser, arg, name):
    if not os.path.isdir(arg):
        parser.error(f"The {name} directory {arg} does not exist!")
    else:
        return arg


def is_valid_scope(parser, arg):
    try:
        return int(arg)
    except ValueError:
        parser.error(f"The value {arg} is not a valid scope!")


def main():
    parser = argparse.ArgumentParser(
        description="bindcompare: Reveal and analyze co-regulatory sites across two protein-binding datasets."
    )

    parser.add_argument(
        "-r",
        "--ref",
        required=True,
        type=lambda x: is_valid_file(parser, x, "reference BED"),
        help="The file path for your reference BED file. Should be the DNA path if comparing DNA and RNA.",
    )
    parser.add_argument(
        "-e",
        "--exp",
        required=True,
        type=lambda x: is_valid_file(parser, x, "overlayed BED"),
        help="The file path for your overlayed BED file.",
    )
    parser.add_argument(
        "-s",
        "--scope",
        required=True,
        type=lambda x: is_valid_scope(parser, x),
        help="Nucleotides upstream and downstream from the ref peak's center that BC will search for peaks.",
    )
    parser.add_argument(
        "-n", "--name", required=True, help="Str name for this experiment."
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        type=lambda x: is_valid_directory(parser, x, "output"),
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

    if args.gtf != "None" and not os.path.isfile(args.gtf):
        parser.error(f"Valid GTF file not provided.")

    if args.fasta != "None" and not os.path.isfile(args.fasta):
        parser.error(f"Valid FASTA file not provided.")

    summary_output = os.path.join(args.out, f"{args.name}_summary.txt")

    if not os.path.exists(args.out):
        os.makedirs(args.out, exist_ok=True)

    with open(summary_output, "w") as summary_file:
        merge.main(args.ref, args.exp, args.scope, args.name, f"{args.out}/", args.gtf)

    with open(summary_output, "a") as summary_file:
        downstream.downstream(
            os.path.join(args.out, f"{args.name}_overlaps.csv"), args.fasta, args.out
        )

    date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"Completed BindCompare! Time Stamp: {date}")
