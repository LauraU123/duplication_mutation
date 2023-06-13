# Duplication mutation analysis

This workflow analyses the evolution of the  duplication in the RSV G gene.

Inputs:
 -  RSV sequences (fasta)
 -  root sequence (json)
 - tree file (nwk)
 
 These inputs can be obtained by running the without-G workflow. [RSV workflow without G](https://github.com/LauraU123/without-G/tree/main)
 This workflow carries out RSV phylogenetic analyses which have trees built on all genes excluding the G gene. 
 Outputs:
	- graph of each mutation copy location - synonymous (png)
	- graph of each mutation copy location - nonsynonymous (png)
	- cumulative distribution for each mutation copy (png)
	- KS goodness of fit test for each distribution (csv)
	- mutation rate of each copy (csv)

## Running the workflow

To run the workflow, run snakemake --cores all from the relevant folder.

To specify A or B, add A or B in the rule all in the Snakefile



