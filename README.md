# BindCompare

This is a novel DNA/RNA Integration tool meant to analyze the overlap between protein binding sites. The input data is in the form of peak-called BED files from programs such as MACS2 and SEACR.

### What Does BindCompare Do?

Oftentimes, regulation from factors occur within a larger locus surrounding the marked binding site. The script searches for overlaps between these binding sites in a chosen scope. Using the Genome File and the Genes GTF File, corresponding Gene IDs are extracted 


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


### Navigating the RShiny

## Understanding the Results


This was script was written at Brown University in the Larschan Lab: larschanlab.com.

