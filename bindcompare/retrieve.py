import gzip
import shutil
import os


def copy_file(source_path, destination_directory="."):
    file_name = os.path.basename(source_path)
    destination_path = os.path.join(destination_directory, file_name)
    shutil.copyfile(source_path, destination_path)


def main():
    with gzip.open(
        os.path.dirname(__file__) + "/reference_files/dm6.fa.gz", "rb"
    ) as file_in:
        with open("dm6.fa", "wb") as file_out:
            shutil.copyfileobj(file_in, file_out)
            print("dm6 fasta file created")

    with gzip.open(
        os.path.dirname(__file__) + "/reference_files/dmel-all-r6.46.gtf.gz", "rb"
    ) as file_in:
        with open("dmel-all-r6.46.gtf", "wb") as file_out:
            shutil.copyfileobj(file_in, file_out)
            print("dm6 GTF file created")

    copy_file(os.path.dirname(__file__) + "/reference_files/dm6.fa.fai")


if __name__ == "__main__":
    print(os.path.dirname(__file__))
