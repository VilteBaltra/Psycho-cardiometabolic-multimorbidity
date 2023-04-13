
# MR MULTIMORBIDITY PAPER
# MR code is taken from: https://mrcieu.github.io/TwoSampleMR/articles/perform_mr.html

# install.packages("remotes")
# remotes::install_github("MRCIEU/TwoSampleMR")
# remotes::install_github("MRCIEU/MRInstruments")
library(TwoSampleMR)
library(MRInstruments)
library(tidyverse)

# MULTIPLE EXPOSURES AT ONCE

# extract SNPs for 7 risk factors with ukb IDs
exp_dat <- extract_instruments(outcomes = c("ukb-a-225", "ukb-b-19953", "ukb-b-8909", "ieu-b-38", "ieu-b-39", "ukb-d-30710_irnt", "ieu-b-116" )) 
table(exp_dat$exposure)
# Body fat percentage || id:ukb-b-8909 
# 395 
# Body mass index (BMI) || id:ukb-b-19953 
# 458 
# C-reactive protein || id:ukb-d-30710_irnt 
# 210 
# diastolic blood pressure || id:ieu-b-39 
# 460 
# Fasting insulin || id:ieu-b-116 
# 14 
# Smoking status: Current || id:ukb-a-225 
# 15 
# systolic blood pressure || id:ieu-b-38 
# 461 

# extract exposure SNPs from outcome dataset
outcome_dat <- read_outcome_data(
  snps = exp_dat$SNP,
  filename = "multimorbidity_Nikpay.Scott.Daner_noUKBB_non-het.txt.gz",
  sep = " ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  id_col = "MM",
  chr_col = "CHR",
  pos_col = "BP"
)

dat <- harmonise_data(
  exposure_dat = exp_dat, 
  outcome_dat = outcome_dat
)
res <- mr(dat)


# Combine all results
het <- mr_heterogeneity(dat)
plt <- mr_pleiotropy_test(dat)
sin <- mr_singlesnp(dat)
all_res <- combine_all_mrresults(res,
                                 het,
                                 plt,
                                 sin,
                                 ao_slc = FALSE,
                                 Exp = TRUE,
                                 split.exposure = FALSE,
                                 split.outcome = FALSE)

head(all_res[, c("Method", "outcome", "exposure", "nsnp",
                 "b", "se", "pval", "intercept", "intercept_se",
                 "intercept_pval", "Q", "Q_df", "Q_pval",
                 # "consortium","ncase","ncontrol","pmid","population"
                 "or", "or_lci95","or_uci95")])

# save it in a csv file 
results_tidy <- all_res[, c("Method", "outcome", "exposure", "nsnp",
                            "b", "se", "pval", "intercept", "intercept_se",
                            "intercept_pval", "Q", "Q_df", "Q_pval",
                            # "consortium","ncase","ncontrol","pmid","population"
                            "or", "or_lci95","or_uci95")] 

write.csv(results_tidy, file = "two-sample-MR-output.csv", quote = F)


##################### PMID ##################### 

# using PMID for insomnia
insomnia_gwas <- subset(gwas_catalog, grepl("Jansen", Author) & Phenotype == "Insomnia"  & Year == 2019 )

exp_data <- format_data(insomnia_gwas) 
exp_data <- clump_data(exp_data) 

# extract exposure SNPs from outcome dataset
out <- read_outcome_data(
  snps = exp_data$SNP,
  filename = "multimorbidity_Nikpay.Scott.Daner_noUKBB_non-het.txtz",
  sep = " ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  id_col = "MM",
  chr_col = "CHR",
  pos_col = "BP")


dat <- harmonise_data(  
    exposure_dat = exp_data, 
    outcome_dat = out)

res <- mr(dat)

# add MR Egger intercept and heterogeneity test 
source("my_MR_tests.R")
res_all <- my_MR_tests(res, dat)

write.csv(res_all, file = "MR-insomnia-output.csv", quote = F)


##################### MR FOR REMAINING TRAITS WAS RUN USING LOCALLY DOWNLOADED SUMSTATS ##################### 

### CHILDHOOD MALTREATMENT ###
CM_gwas <- read.delim("Retro_prospective_meta_childhoodmaltreatment.txt.gz", sep = " ")

# format exposure data
CM_gwas$Phenotype <- "Maltreatment"
CM_gwas$N <- 185414

# reduce the size of the dataset to only suggestive SNPs (makes it faster to clump later)
CM_exp_dat <- CM_gwas %>% filter(P < 0.000005) # 1698

CM_exp_dat <- format_data(CM_exp_dat,
                          type = "exposure",
                          snp_col = "SNP",
                          beta_col = "BETA",
                          se_col = "SE",
                          effect_allele_col = "A1",
                          other_allele_col = "A2",
                          pval_col = "P",
                          samplesize_col = "N",
                          min_pval = 1e-200,
                          #z_col = "Z",
                          #info_col = "INFO",
                          chr_col = "CHR",
                          pos_col = "BP") 


CM_exp_dat <- clump_data(CM_exp_dat) 

# extract exposure SNPs from outcome dataset
CM_out <- read_outcome_data(
  snps = CM_exp_dat$SNP,
  filename = "multimorbidity_Nikpay.Scott.Daner_noUKBB_non-het.txtz",
  sep = " ",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  id_col = "MM",
  chr_col = "CHR",
  pos_col = "BP"
)

dat_cm <- harmonise_data(
  exposure_dat = CM_exp_dat, 
  outcome_dat = CM_out
)
res_cm <- mr(dat_cm) 

# add MR Egger intercept and heterogeneity test 
res_cm_all <- my_MR_tests(res_cm, dat_cm)


### and same as above for all remaining traits ### 

# save output
write.csv(res_cm_all, file = "MR-with-local-sumstats.csv", quote = F)



