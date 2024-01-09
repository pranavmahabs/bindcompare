``bindexplore``
===============

If you are interested in taking N BED files for N different binding
experiments, you can visualize pair-wise binding overlaps across all
experiments to select candidate pairs for BindCompare. This script,
right now, can only be run from the terminal:

.. code-block:: bash

   bindexplore <scope> <bed_file_1> <bed_file_2> ... <bed_file_n>

The ``scope`` value essentially bins the genome into bins of size
``scope``. Then, it uses this size to search for overlaps within each
bin. Then a heatma is generated to visualize binding overlaps and can be
seen below. The math for each cell is as follows:

.. math:: 

   \frac{ \text{Num Ref Binds found in Overlayed Sites}}{\text{Num Ref Binds}}

.. image:: ./images/explore.png
   :align: center
   :width: 250

In this example, we see that we are comparing CLAMP binding in KC and S2
Cells, GAF Binding, and MLE Binding. Understandably, CLAMP KC and S2 has
a significant overlap!