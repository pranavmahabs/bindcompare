import sys
import os

from .bindapp import merge
from .bindapp import downstream


def main():
    if len(sys.argv) != 8:
        print(
            f"\nUsage: bindcompare <Ref BED Path> <Exp BED Path> <Scope> <Sample Name> <Output Directory> <GTF File> <Genome FA>\n\n"
            + "Ref BED Path: The file path for your reference BED file. Should be the DNA path if comparing DNA and RNA.\n"
            + "Exp BED Path: The file path for your overlayed BED file.\n"
            + "Scope: Nucleotides upstream and downstream from the ref peak's center that BC will search for peaks.\n"
            + "Sample Name: str name for this experiment.\n"
            + "Output Directory: All outputs for bindcompare will be found here. Will be created if DNE.\n"
            + "GTF File: Gene GTF file in proper format. Enter None if not available.\n"
            + "Genome FA File: Genome file corresponding to your BED Files.\n\n"
            + "If using D. Melanogaster, run retrievedm6 to get GTF and Genome file. Further instructions at https://github.com/pranavmahabs/bindcompare.\n"
        )
        sys.exit(1)

    DNA = sys.argv[1]
    RNA = sys.argv[2]
    SCOPE = sys.argv[3]
    SNAME = sys.argv[4]
    OUT = sys.argv[5]
    GTF = sys.argv[6]
    FASTA = sys.argv[7]

    summary_output = os.path.join(OUT, "{SNAME}_summary.txt")

    if not os.path.exists(OUT):
        os.makedirs(OUT, exist_ok=True)

    with open(summary_output, "w") as summary_file:
        merge.main(DNA, RNA, SCOPE, SNAME, f"{OUT}/", GTF)

    with open(summary_output, "a") as summary_file:
        downstream.downstream(os.path.join(OUT, "_overlaps.csv"), FASTA, OUT)

    print("Completed BindCompare! Time Stamp:")
    os.system("date")
