
# MRLap
# Code is from: https://github.com/n-mounier/MRlap
# Directly install the package from github
# install.packages("remotes")
# remotes::install_github("n-mounier/MRlap")
library(MRlap)

# read in multimorbidity GWAS data
outcome <- read.delim("multimorbidity_Nikpay.Scott.Daner_noUKBB_non-het.txt.gz")
# add missing info (or rename it)
outcome$Neff <- 562507
outcome$POS <- outcome$BP

# read in exposure GWAS
exposure <- read.delim("SESA_neuro_clus_sumstats.txt.gz")

# rename columns
exposure$N <- exposure$N_analyzed

# run MRlap
A = MRlap(exposure = exposure,
          exposure_name = "SESA", # replace with other traits
          outcome = outcome,
          outcome_name = "multimorbidity",
          ld = "eur_w_ld_chr",
          hm3 = "w_hm3.snplist")

# explore results
str(A)

# # observed effect
A[["MRcorrection"]]$observed_effect

# corrected effect
A[["MRcorrection"]]$corrected_effect

# difference p-value
A[["MRcorrection"]]$p_difference


### repeat above steps with the remaining traits ###


