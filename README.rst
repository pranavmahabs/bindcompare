

.. image:: https://pranavmahabs.github.io/bindcompare/docs/images/bclogo.png
   :align: center
   :width: 300

=================================================================
BindCompare: A novel integrated protein-binding analysis platform
=================================================================


.. _Installation Guide: https://github.com/pranavmahabs/bindcompare/blob/main/docs/installation.rst
.. _bindcompare: https://github.com/pranavmahabs/bindcompare/blob/main/docs/bindcompare.rst
.. _bindexplore: https://github.com/pranavmahabs/bindcompare/blob/main/docs/bindexplore.rst
.. _Downstream & comparexp: https://github.com/pranavmahabs/bindcompare/blob/main/docs/comparexp.rst

`Introduction`_ 
| `Installation Guide`_ 
| `Using BindCompare`_ 
| `Quick Start`_ 
| `Credits`_

BindCompare Function Manuals: 
`bindcompare`_
| `bindexplore`_
| `Downstream & comparexp`_


Introduction
============

BindCompare is a novel integrated protein-binding analysis platform
designed to be user-accessible and interpretable. Given protein-binding
sites on DNA and/or RNA, BindCompare determines and visualizes domains
of co-regulatory activity at the single-nucleotide level.

At the core of BindCompare is defining overlapping binding domains. 
Oftentimes, co-regulation with factors occurs across a larger locus 
surrounding the marked binding site. BindCompare searches and categorizes
overlaps across a scoped regions containing the reference binding sites. 
Because BindCompare utilizes BED files, it enables the comparison between 
RNA and DNA binding sites, aiding the study of system wide 
co-transcriptional regulation.

To support user-accessible software, BindCompare can be launched in a GUI
interface that allows for easy application of bindcompare and comparexp. To
learn how to utilize the tool, please read below and do not hesitate to 
reach out if any bugs or issues are revealed.

Quick Start
===========

Install bindcompare and its dependencies from PyPI using pip::

   pip install bindcompare

Then, to run the core function ``bindcompare``::

   bindcompare -r REF -e EXP -s SCOPE -n NAME -o OUT [-g GTF] [-f FASTA]

Alternatively, to launch the core GUI application, see below. For 
a more detailed walk-through, read the full documentation below and function specific
manuals. Complete installation instructions can also be found at 

Using BindCompare
=================

Here is a general overview of BindCompare usage in the following schematic.

.. image:: https://pranavmahabs.github.io/bindcompare/docs/images/schematic.png
   :align: center
   :width: 400

To start, optionally run `bindexplore`_ to find candidate co-regulators. Then, you
can choose two candidate coregulators and run `bindcompare`_ to explore co-regulatory
activities between the protein-binding datasets given. Finally, you can look into 
downstream analysis using ``comparexp`` or gene ontology/motif analysis. 

As aforementioned, ``bindcompare`` and ``comparexp`` can be run through a tkinter
GUI interface. All of the commands can be run from the command line. This includes
``retrievedm6`` and ``bindexplore``. How to use BindCompare is presented for both 
the GUI and command line approaches. To launch the GUI:

.. code-block:: bash

   bindlaunch

That should launch a platform that looks like this:

.. figure:: https://pranavmahabs.github.io/bindcompare/docs/images/bindlaunch.png
   :align: center
   :width: 350

Please visit the command specific pages for each of the above commands to learn more
about how to use BindCompare. 

Credits
=======

This script was written at Brown University in the `Larschan
Lab <https://www.larschanlab.com>`__ by Pranav Mahableshwarkar under
the guidance of Mukulika Ray, PhD and Erica Larschan, PhD. 

If you want to pull the source-code, this can be done via github. 

.. code-block:: bash

   git pull https://github.com/pranavmahabs/bindcompare.git

Please leave any messages here regarding errors or issues found in using the platform. 
