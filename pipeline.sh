#!/bin/bash 

# analysis is meant to be run in the unzipped repository
R_version="/opt/R/4.0.3/bin/Rscript"
script_dir="./R"

# setup R virtual enviornmment for reproducibility
${R_version} "${script_dir}/setup_renv.R"

# read in array data and perform normmalization
${R_version} "${script_dir}/preprocessing.R" \

# coonduct differential gene expressioon for all contrasts of interest
${R_version} "${script_dir}/differential_gene_expression.R" \

