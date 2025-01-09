#!/usr/bin/env bash

## SEM analysis packages
Rscript -e 'install.packages("picante", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("plyr", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("bayesplot", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("pscl", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("brms", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("performance", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("R2admb", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("shinystan", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("ggplot2", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("dplyr", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("tidyr", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("tidybayes", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("ggthemes", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("bayestestR", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("gtools", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("viridis", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("gridExtra", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("grid", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("RColorBrewer", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("rstantools", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("phyloseq", repos="http://cran.r-project.org")'

## Network Turnover analysis packages
Rscript -e 'install.packages("lme4", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("lmerTest", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("igraph", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("ggpubr", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("emmeans", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("bipartite", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("fields", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("betapart", repos="http://cran.r-project.org")'

## Map figure packages
Rscript -e 'install.packages("sf", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("ggrepel", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("terra", repos="http://cran.r-project.org")'
Rscript -e 'install.packages("ggspatial", repos="http://cran.r-project.org")'
Rscript -e 'devtools::install_github("yutannihilation/ggsflabel")'


