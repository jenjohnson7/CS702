#!/bin/bash
echo "Requires AD/AR/XLR_clinvar_results.txt"
Rscript get_allele_id_to_query.R
python get_ac_af.py
Rscript writing_obs.R
