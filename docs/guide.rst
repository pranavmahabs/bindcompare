BindCompare User/Interpretation Guide
======================

This guide will walk you through the BindCompare workflow, how to choose a scope, how to interpret BindExplore results, and how to understand the overlap profiles.

Workflow
########

The BindCompare workflow consists of three main steps. Make sure that you have properly peak-called your data prior to using the BindCompare platform. If you are interested in gene/sequence information, make you have the appropriate respective genome GTF file or FASTA file.

To start, optionally run ``bindexplore`` to find candidate co-regulators. Then, you
can choose two candidate coregulators and run ``bindcompare`` to explore co-regulatory
activities between the protein-binding datasets given. Finally, you can look into 
downstream analysis using ``comparexp`` or gene ontology/motif analysis. 

Make sure to look at the relevant docs page (accessible from the home page) to determine how to run these methods. 

Interpreting and Using BindCompare
##################################

Interpreting ``bindexplore``
----------------------------

BindExplore generates visualizations and data that help you understand the binding interactions:

Choosing a Bin Size
+++++++++++++++++++
Choosing a bin size makes a large impact on the quality of your results. Generally, it is recommended to try a bin size of around 1000bp and determine the quality of your results. 

If the yielded results seem overly diluted (i.e. you are capturing too many overlaps) decreasing the bin size will decrease the number of captured overlaps. Conversely, if the bin size is too small, correlation scores will likely be very low and will not reveal combinatorial patterns.

Normalization and Binning
+++++++++++++++++++++++++

``bindexplore`` normalizes the interaction scores by scaling counts based on the number of reference peaks. This prevents skewing results towards proteins with low specificity. It also allows us to compare all pairs within the same 0-1 range. 

If two peaks from different datasets fall within the same bin (e.g., 5 kb) but do not strictly overlap, they are still considered in the intersection. Bins are used for intersecting lists rather than individual peaks, which helps in summarizing the binding interactions more effectively.

Asymmetry in Correlation Matrix Scores Between Protein Pairs
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

An asymmetric correlation matrix scores represents differential interaction between multiple proteins. Multiple TFs and RBPs can form different hubs with differential interacting partners. To identify different hubs of proteins binding to common nucleic acid targets, the user can leverage the correlation matrix scores in `bindexplore`. 

The difference in score is due to normalization specific to each protein. When comparing A (reference) vs. B (overlapped), the number of co-occurrences across bins are normalized against the number of peaks in A, and vice versa. This allows for accurate proportional comparisons, highlighting differences in binding specificities between the two datasets.

.. image:: ./images/explore.png
   :align: center
   :width: 400

For example, in the figure above, we have run ``bindexplore`` on 17 nucleic-acid binding proteins. When Protein C is the overlapped binding protein, correlation matrix values are quite high (see row 3). However, when Protein C is the reference (column 3), correlation matrix values are much lower. This indicates that HNRNPA1 has low specificity and bins along with many different RBPs but few RBPs, such as Protein C, appear alongside all HNRNPA1 sites with high frequency. 

Interpreting ``bindcompare``
----------------------------

The overlap profile plot is a key output of ``bindcompare``, providing detailed insights into binding interactions. A comprehensive interpretation guide is provided in the ``bindcompare`` `docs page <https://github.com/pranavmahabs/bindcompare/blob/main/docs/bindcompare.rst>`__. 

Choosing a Scope
++++++++++++++++
We have used a default of 1000bp scope region when comparing DNA peaks. While comparing RNA peaks (smaller size), we used a smaller scope area of 250bp. We suggest that the user choose a scope/bin size reflecting the use case. Selecting a larger bin size would pick up more distal co-binding, while a smaller bin size would be more stringent. In our experience, the default 1000bp scope region works for DNA-DNA and DNA-RNA comparisons. The scope size of 250bp works for RNA-RNA comparisons.

Overlap Profiles
++++++++++++++++

Here we will present two overlap profiles produced from `bindcompare` for Protein A and Protein B. In the first image, Protein A is the reference, whose average DNA binding peaks are denoted by the black line. The colored lines indicate the varying overlap categories described in the docs page. 

.. image:: ./images/taf15_fus_dna.png
   :align: center
   :width: 220

Here, we present the DNA binding peaks as a reference for Protein A. Overlaid, we have the RNA binding peaks for Protein A. 

.. image:: ./images/taf15_dnarna.png
   :align: center
   :width: 220

Looking at the first plot, we can see that a large number of Protein B binding sites bind towards the beginning or end of a Protein A binding site with much fewer occurring in the extended scoped region. This indicates that Protein B positions itself near Protein A binding. On the other hand, when comparing Protein A DNA binding to its RNA binding, we see that while there are 100s of sites where DNA and RNA binding overlap, there is no clear skew along the scoped region.
