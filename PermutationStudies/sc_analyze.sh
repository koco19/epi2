#!/bin/bash

# Generate figures in R

# Number of permutations to perform
n_perm=10000

mkdir -p out

# this is not really used, but maybe good as comparison
Rscript -e "rmarkdown::render('Perm-analysis.Rmd', output_format='html_document', params=list(n_perm='${n_perm}'), output_file='out/Perm-analysis-${n_perm}.html')"
# the relevant data
Rscript -e "rmarkdown::render('Perm-analysis2.Rmd', output_format='html_document', params=list(n_perm='${n_perm}'), output_file='out/Perm-analysis2-${n_perm}.html')"

# Package
rm -f out-Perm.zip
zip -r out-Perm.zip out
