# dnarna_bindcompare

This is a novel DNA/RNA Integration tool meant to analyze the overlap between RNA and DNA protein binding sites. The input data is in the form of BED files from either MACS2 or SEACR. 

Oftentimes, regulation from factors occur within a larger locus surrounding the marked binding site. The script searches for overlaps between these binding sites in a chosen scope. The script will soon support a statistical significance metric to compare found overlaps between two different sample types. Additionally, the script will soon contain Gene Ontology analysis and motif analysis for the genes and regions found in the overlaps. 

The package includes a shell script, cnr_iclip_merge.sh, to run the python file in addition to sample outputs. 

This was script was written at Brown University in the Larschan Lab: larschanlab.com.
