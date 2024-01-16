Downstream & ``comparexp``
==========================

``comparexp``
-------------
This is a simple way to compare and visualize two ``bindcompare`` runs. 
This is only possible if you provided a GTF file. Make sure that your
original summary file is in the ``bindcompare`` output directory you provide.

Given two directories for ``bindcompare`` outputs, ``comparexp`` will produce
a weighted Venn Diagram based on the genes present in each of the categories.
In a summary file, the program will also provide the jaccard similarity score
in addition to the names of the genes in each category (comma separated so that
it can be easily input into Gene Ontology platforms).

The output names for these files will have the prefix "samplename1_v_samplename2"
where each sample name is extracted from the name of the summary file. They will
be saved in the folder from which you launched the GUI or where you run ``comparexp``.

Command Line
^^^^^^^^^^^^
The same functionality can also be achieved through the command line. 

.. code-block:: bash

   comparexp [-a outputpath_1] [-b outpath_2]
   
Using the GUI
^^^^^^^^^^^^^
When you launch the app using ``bindlaunch``, you may click the `comparexp`
button. This will launch another mini window. 

.. image:: ./images/empty.png
   :align: center
   :width: 200

There, you can enter in the file paths for two folders that contain BindCompare summary 
files with gene lists and click run. The output and GUI should look like this:

.. image:: ./images/comparexp.png
   :align: center
   :width: 200

As you can see, the venn diagram is rendered and the entire summary file
will also be dumped below in the text box. 

Gene Ontology and Motif Analysis
--------------------------------
If you provided a GTF file then your summary file would contain a gene list. You can find
overlap type specific gene lists as well in the produced CSV files. This comma separated
gene list can be fed into popular gene ontology packages including ShinyGO or GProfiler2.
You would also get similar gene lists from ``comparexp`` outputs!

If you provided a genome FA file, then ``bindcompare`` produces a sequences.fa that can be
plugged into various different programs within the MEMESuite such as MEME or STREME or FIMO.