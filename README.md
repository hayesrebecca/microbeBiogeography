
## Description of manuscript  
This repository contains all of the necessary data products and code to reproduce analyses and results from the manuscript titled "Drivers of distance decay of host-associated microbiomes vary by host sociality and microbial association strength". We analyzed wild bee gut microbe communities across mountain meadows in the Sky Islands region of Arizona and New Mexico to understand how host sociality and microbial association strength influence the distance decay of similarity relationship. We used Bray-Curtis dissimilarity of individual host gut communities and interaction dissimilarity approaches to understand these relationships at the individual and whole host community levels, and then decomposed interaction dissimilarity into its components (rewiring, host-driven turnover, microbe-driven turnover, and complete turnover) to explore the mechanisms of spatial dissimilarity.

## Repository organization  

packages.sh  
- This shell script includes code to install all of the required R packages to run analysis and figure-generating scripts. 

microbeBiogeographyData.Rdata  
- This Rdata file includes all of the necessary pieces of data to run the scripts in this repository. These include:
  - spec.net  
    - This data frame includes all of the specimen-level data and relevant explanatory variables for all analyses. Variables used for the SEM analysis are transformed in this data frame (see analysis/SEM/1_microbeSEM.R for a list of transformed variables)
  - spec.orig  
    - This data frame is the same as spec.net, except that it stores the untransformed values to facilitate plotting on the original scale.
  - tree.16s  
    - This phylo object is the full microbial phylogenetic tree for all ASVs included in this manuscript.
  - site.sum  
    - This data frame stores the site-level variables.
  - abund_csv  
    - This data frame stores the abundance of each species of bee across year, site, and sample round.
  - indivNet_micro  
    - This is a list of individual-level bee-microbe networks at each site, with each item in the list representing that site's individual bee-microbe network.
  - spNet_micro  
    - This is a list of species-level bee-microbe networks at each site, with each item in the list representing that site's species-level bee-microbe network.
  - indiv.comm.16sR0  
    - This matrix array represents the community matrix of individual bee hosts as rows, microbial ASVs as columns, and total read abundance as values.
  - finalASVtable  
    - This matrix array represents the community matrix of individual bee hosts as rows, microbial ASVs as columns, and relative abundance of those ASVs in each individual as values.
  - physeq16sR0  
    - This phyloseq object is the is the full microbial phylogenetic tree for all ASVs included in this manuscript.
  - phylo_matrix  
    - This matrix array is the phylogenetic covariance matrix of bee species included in this manuscript. The original published supermatrix was downloaded and trimmed using our species list. (Henriquez Piskulich, Patricia Andrea; Hugall, Andrew F.; Stuart-Fox, Devi (2023). A supermatrix phylogeny of the world’s bees (Hymenoptera: Anthophila) [Dataset]. Dryad. https://doi.org/10.5061/dryad.80gb5mkw1
    )  
  - sites_sf  
    - This sf spatial object includes spatial data from our sites used for generating the map figure.  
  - feature.2.tax.16s  
    - This data frame includes the taxonomic assignment for each 16s feature to facilitate generation of phylogenetic tree figures.

analysis  
  - This folder includes the scripts to run the two main analyses for this manuscript, including:  
    - NetworkTurnover  
      - 1_NetworkTurnover.R  
        - This script runs the network turnover components and dissimilarity by geographic distance analyses and generates the results figures.
      - 2_brayCurtisDissimilarity.R
        - This script runs the microbial community dissimilarity and interaction dissimilarity by geographic distance analyses and generates the results figures.
      - 3_beeCommunityDissimilarity.R
        - This script runs the bee community dissimilarity by geographic distance analysis and generates the results figure from the supplement.
      - 4_plantCommunityDissimilarity.R
        - This script runs the plant community dissimilarity by geographic distance analysis and generates the results figure from the supplement. 
  
figures  
  - This folder includes the scripts to generate the site map and microbial phylogenetic tree figures.  

## To use this repository  
  1. Run packages.sh to install all necessary packages.  
  1. Scripts are self-contained and can be run out of order.
