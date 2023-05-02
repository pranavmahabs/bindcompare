# BindCompare

This is a novel DNA/RNA Integration tool meant to analyze the overlap between protein binding sites. The input data is in the form of peak-called BED files from programs such as MACS2 and SEACR.



Oftentimes, regulation from factors occur within a larger locus surrounding the marked binding site. The script searches for overlaps between these binding sites in a chosen scope. 

![Illustration of Binding Site Overlapping in Scoped Region](https://github.com/pranavmahabs/bindcompare/blob/main/BindCompareDemo1.png)

As we can see in this diagram, there are two main categories of overlaps that can be found from two BED Files. When you run the program, you will choose one BED file to the reference, or base, BED file and one BED File to be overlayed. The first time are complete overlaps in peak sites which can be seen in the right half of the drawing. BindCompare also looks for overlaps upstream and downstream of the reference peak. On the left half of the drawing, we see that there is an overlap outside of the direct binding site. These categories are further broken down and are explained in the Results category. 

Afterwards, Gene IDs and sequences are extracted from these sites. Then, Gene Ontology analysis is conducted on those genes and MEME/STREME motif analysis is conducted on the sequences. Finally, you can take gene lists from two separate comparisons and see if there are overlaps across different runs of BindCompare!

## Quick Start: 3 Steps
  1. Install the conda environment using the following command in your BindCompare directory:
  ```
  ~/bindcompare $ conda env create -f environment.yml
  ```

  2. Activate BindCompare environment using the following command:
  ```
  ~/bindcompare $ conda activate bindcompare
  ```

  3. Then, to launch the RShiny Web Application, run the following command: 
  ```
  ~/bindcompare $ R -e "shiny::runApp('bindapp’)”
  ```
  This will print out an http:// link that you can copy and paste into your browser to run the appliation! This will look something like: "Listening on http://127.0.0.1:3868". 
  
  Do not terminate this process in your terminal until you are done running BindCompare. Status results from your experiments will appear here as you go. 

## Navigating the RShiny

### Comparing Two Bed File
In the first tab of the RShiny, you will be able to run the core functionality of BindCompare. There are seven input boxes on the left that you will have to fill out before running the tool. 
1. **BindCompare Directory Filepath:** Enter the directory that contains your copy of BindCompare. You can get this by running `pwd` in your Terminal/CL while in the `bindcompare` folder. 
2. **Base Bed File Path:** Enter the file path for your reference BED file. If comparing DNA and RNA, then this should be the filepath for the DNA BED file or more generally, the BED file with the larger peak size. 
3. **Overlayed Bed File Path:** Enter the file path for your overlayed BED file. Conversely, when applicable, this would be the BED file with the smaller peak size. 
4. **Scope:** The scope is how many nucleotides upstream and downstream from the reference peak's center that BindCompare will search for an overlap. Making this value smaller will decrease the number of overlaps and vice versa. 
5. **Sample Name:** A short phrase to label the experiment (i.e. CLAMP)
6. **Output Folder:** An existing folder's file path where all of the outputs will be generated. 
7. **Genes GTF File:** This file details the chrom location of every gene in your organism. The GTF file for D. Melanogaster is provided (gzipped) in the reference category. 
8. **Genome FA File Path:** A FA file with a corresponding fa.fai (index file) for BedTools to extract sequences of binding sites and perform motif analysis. Enter `None` to skip this feature!

### Comparing Two BindCompare Experiments

## Understanding the Results
### Overlap Profile 
![Overlap Profile for S2](https://github.com/pranavmahabs/bindcompare/blob/main/SampleOut/S2_overlaps.png)
### Bar Graph and Pie Chart
![Bar Graph and Pie Chart for KC](https://github.com/pranavmahabs/bindcompare/blob/main/SampleOut/S2_bartotals.png)
### Summary File and CSV Output
The CSV file contains one row for every reference peak that was involved in an overlap. This includes the Chromosome, Beginning/Ending Coordinate of the peak, the corresponding nucleotide sequence, the type of overlap (as described above), and the Gene IDs that correspond to that region. 

The summary file contains the average peak size for both of the BED files. Additionally, it prints all of the found Gene IDs that are in the CSV file so that they can be easily converted to gene names. 

### Other Outputs!
Gene Ontology results from GProfiler2 and motif analysis from either STREME or MEME or also included in this directory. Please see the [MEME Suite](https://meme-suite.org/meme/doc/streme.html) page for more information on MEME/STREME. Please see the [GProfiler2 Manual](https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html) for more information on the Gene Ontology analysis! Feel free to copy the gene list into your GO tool of choice as well!

## Credits
This was script was written at Brown University in the Larschan Lab: larschanlab.com.

