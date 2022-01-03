# # Support code for Ncube-Kanyika *et.al.* (2022) CASSAVA PROTEINS INTERACTING WITH THE MOVEMENT PROTEIN OF SOUTH AFRICAN CASSAVA MOSAIC VIRUS REVEAL ROLES IN PATHOGENICITY.  

By Valeria Velásquez-Zapata


**Overview**

This repository contains code for performing the bioinformatic analysis of the manuscript. It is separated in two R scripts with instruction on how to format the files to run the SAINTexpress and Y2H-SCORES software in analysis of MS Co-IP data and for posterior network analysis . 


**Citations**

If you use part of this code please cite  

* Ncube-Kanyika B, Velásquez-Zapata V, Wise RP, Rey MEC (2022). Cassava proteins interacting with the movement protein of south African cassava mosaic virus reveal roles in pathogenicity.


**Software requirements**

* R
* R packages: UniprotR, stringr, reshape2, tidyverse, psych, DESeq2, mass, igraph, RCy3
* Unix

**Functions in this repository**

Two R files were created with the main functions used to perform the analysis. A brief description of each file is provided: 

1. `SAINTexpress and Y2H-SCORES analyses.R` contains the formatting of files to run SAINTexpress and Y2H-SCORES.

2. `BC1_subnetwork_analysis.R` contains network analysis of the BC1 interactors in a compiled Arabidopsis interactome. 

