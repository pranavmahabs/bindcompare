.. raw:: html

   <p align="center">

.. raw:: html

   </p>

This is a novel DNA/RNA Integration tool meant to analyze the overlap
between protein binding sites. The input data is in the form of
peak-called BED files from programs such as MACS2 and SEACR. If you are
interested in selecting candidate co-regulators from a set of BED files,
please go down to the `BindExplore <#BindExplore>`__ section.

Oftentimes, co-regulation from factors occurs within a larger locus
surrounding the marked binding site. The script searches for overlaps
between these binding sites in a chosen scope. Because BindCompare
utilizes BED files, it enables the comparison between RNA and DNA
binding sites, aiding the study of system wide co-transcriptional
regulation.

.. raw:: html

   <!-- ![Illustration of Binding Site Overlapping in Scoped Region](https://github.com/pranavmahabs/bindcompare/blob/main/BindCompareDemo1.png) -->

.. raw:: html

   <p align="center">

.. raw:: html

   </p>

As we can see in this diagram, there are two main categories of overlaps
that can be found from two BED Files. When you run the program, you will
choose one BED file to the reference, or base, BED file and one BED File
to be overlayed. The first time are complete overlaps in peak sites
which can be seen in the right half of the drawing. BindCompare also
looks for overlaps upstream and downstream of the reference peak. On the
left half of the drawing, we see that there is an overlap outside of the
direct binding site. These categories are further broken down and are
explained in the Results category.

Afterwards, Gene IDs and sequences are extracted from these sites. Then,
Gene Ontology analysis is conducted on those genes and MEME/STREME motif
analysis is conducted on the sequences. Finally, you can take gene lists
from two separate comparisons and see if there are overlaps across
different runs of BindCompare!

Quick Start: 4 Steps
--------------------

1. Clone the Git Repository for BindCompare (if you do not have this on
   your machine, go
   `here <https://github.com/git-guides/install-git>`__).

::

   ~/home $ git clone https://github.com/pranavmahabs/bindcompare.git

2. Install the conda environment using the following command in your
   BindCompare directory:

::

   ~/home $ cd bindcompare
   ~/home/bindcompare $ conda env create -f environment.yml

If you are on an M1 device, you will have to use a Rosetta enabled
terminal as many of these packages rely on an x86_64 architecture!

3. Activate BindCompare environment using the following command:

::

   ~/bindcompare $ conda activate bindcompare

4. Then, to launch the BindCompare app, run the following command:

::

   ~/bindcompare $ ./run.sh

This should launch the application in a separate window. The status
results from your experiments will appear in the **terminal** as you go.

Note: If you get a permission denied error when running ``./run.sh``,
running ``chmod +wx run.sh`` should fix this issue!

Reference Files
~~~~~~~~~~~~~~~

To run the script, you are required to provide a Genes GTF file and
optionally a whole Genome FA file. If you are using the DM6 system,
these files are provided - zipped - in the ``reference_files`` folder.
To unzip them and use them, run the following command:

::

   ~/bindcompare/reference_files $ gunzip dm6.fa.gz; gunzip dmel-all-r6.46.gtf.gz

Navigating the Tkinter Application
----------------------------------

Comparing Two Bed Files
~~~~~~~~~~~~~~~~~~~~~~~

In the first frame of the app, you will be able to run the core
functionality of BindCompare. There are seven input boxes on the left
that you will have to fill out before running the tool. 1. **Base Bed
File Path:** Enter the file path for your reference BED file. If
comparing DNA and RNA, then this should be the filepath for the DNA BED
file or more generally, the BED file with the larger peak size. 2.
**Overlayed Bed File Path:** Enter the file path for your overlayed BED
file. Conversely, when applicable, this would be the BED file with the
smaller peak size. 3. **Scope:** The scope is how many nucleotides
upstream and downstream from the reference peak’s center that
BindCompare will search for an overlap. Making this value smaller will
decrease the number of overlaps and vice versa. 4. **Sample Name:** A
short phrase to label the experiment (i.e. CLAMP) 5. **Output Folder:**
An existing folder’s file path where all of the outputs will be
generated. 6. **Genes GTF File:** This file details the chrom location
of every gene in your organism. The GTF file for D. Melanogaster is
provided (gzipped) in the reference category. Enter ``None`` to skip
this feature! 7. **Genome FA File Path:** A FA file with a corresponding
fa.fai (index file) for BedTools to extract sequences of binding sites
and perform motif analysis. Enter ``None`` to skip this feature!

Comparing Two BindCompare Experiments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you run BindCompare with two BED Files, you will get a list of
genes where there was a binding overlap. If you were to run BindCompare
on say male samples and female samples, you would then have two lists of
genes. You can copy and paste these two lists in to the input categories
in the **Comparison Tab** in the RShiny interface and then click submit.
These are the outputs from this analysis:

1. Using the equation :math:`J(A,B) = \frac{|A \cap B|}{|A \cup B|}`,
   the Jaccard Similarity Index is calculated and printed.
2. A size-biased venn diagram is also generated using the R-Eulerr
   package.
3. Gene lists are also printed from each of the following categories:
   Only List 1, Only List 2, Both List 1 and List 2.

Note that this is not saved to your computer and you would need to take
a Screenshot to save this result! Additionally, make sure to copy the
list exactly as it is printed from the first tab.

Understanding the Results
-------------------------

Overlap Profile
~~~~~~~~~~~~~~~

Below is a sample Overlap Profile. The overlaps are categorized into
four main categories based upon the location of the overlap: 1.
Completely overlapping (purple lines in frequency plot). 2. Partially
overlapping at the DNA peak start site (red lines in frequency plot). 3.
Partially overlapping at the DNA peak end site (blue lines in frequency
plot) 4. Non-overlapping, i.e. when there is an overlap in a region
outside the DNA binding site (yellow lines in frequency plot).

This extended region is defined by the scope variable in the script,
allowing the overlap to look for binding sites in the proximity of the
DNA binding site (this scope is 2 kb including the DNA binding site). It
should be noted that multiple RNA peaks can be found on one DNA peak.
All of these overlaps are placed onto a [-scope, scope] region. Then,
each type of overlap shown with a different color is overlaid and
plotted onto a frequency plot. So, if the frequency at a given base pair
is 5, then there are five overlaps that contained that base pair within
the region defined by the scope.

.. raw:: html

   <p align="center">

.. raw:: html

   </p>

Oftentimes, it can be valuable to see **where** this split is occurring!
The values derived for the above plot can be `split over all
chromsomes <https://github.com/pranavmahabs/bindcompare/blob/main/SampleOut/ClampKC_chrom_ref_freq.png>`__.

Bar Graph and Pie Chart
~~~~~~~~~~~~~~~~~~~~~~~

Total Binding Peaks references the number of peaks or rows that are in
the overlayed bed file. Unique overlaps references the number of unique
peaks in the overlayed BED file that were found to overlap with a peak
in the base/reference BED file. The total number of overlaps simply
references how many times an RNA peak overlapped with a DNA peak. Note
that there can be repeats here! Finally, the last column is the number
of unique reference/base peaks that were found within an overlap.

.. raw:: html

   <p align="center">

.. raw:: html

   </p>

Summary File and CSV Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The CSV file contains one row for every reference peak that was involved
in an overlap. This includes the Chromosome, Beginning/Ending Coordinate
of the peak, the corresponding nucleotide sequence, the type of overlap
(as described above), and the Gene IDs that correspond to that region.

The summary file contains the average peak size for both of the BED
files. Additionally, it prints all of the found Gene IDs that are in the
CSV file so that they can be easily converted to gene names.

Other Outputs!
~~~~~~~~~~~~~~

Gene Ontology results from GProfiler2 and motif analysis from either
STREME or MEME or also included in this directory. Please see the `MEME
Suite <https://meme-suite.org/meme/doc/streme.html>`__ page for more
information on MEME/STREME. Please see the `GProfiler2
Manual <https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html>`__
for more information on the Gene Ontology analysis! Feel free to copy
the gene list into your GO tool of choice as well!

BindExplore
===========

If you are interested in taking N BED files for N different binding
experiments, you can visualize pair-wise binding overlaps across all
experiments to select candidate pairs for BindCompare. This script,
right now, can only be run from the terminal:

::

   $ ./bindexplore.sh
   $ Enter BED File Paths (Space-Separated): CLAMP_KC_DNA.bed CLAMP_S2_DNA.bed gaf_chip.bed MLE_DNA.bed
   Provided BED File Paths:
   CLAMP_KC_DNA.bed
   CLAMP_S2_DNA.bed
   gaf_chip.bed
   MLE_DNA.bed
   $ Enter the scope: 5000
   Scope: 5000
   $ Is everything okay? Enter 'yes' to continue or 'no' to cancel: yes

The ``scope`` value essentially bins the genome into bins of size
``scope``. Then, it uses this size to search for overlaps within each
bin. Then a heatma is generated to visualize binding overlaps and can be
seen below. The math for each cell is as follows:

:math:`\frac{\text{Num Ref Binds found in Overlayed Sites}}{\text{Num Ref Binds}}`

.. raw:: html

   <p align="center">

.. raw:: html

   </p>

In this example, we see that we are comparing CLAMP binding in KC and S2
Cells, GAF Binding, and MLE Binding. Understandably, CLAMP KC and S2 has
a significant overlap!

Credits
-------

This was script was written at Brown University in the `Larschan
Lab <https://www.larschanlab.com>`__ by Pranav Mahableshwarkar.