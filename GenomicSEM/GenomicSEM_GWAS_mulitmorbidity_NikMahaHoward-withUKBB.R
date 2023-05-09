
# A multivariate GWAS of psycho-cardiometabolic multimorbidity - version WITH UKBB
# Genomic SEM tutorial available at https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects 

# install and load packages
#install.packages("devtools")
library(devtools)
#install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
packageVersion("GenomicSEM")
# ‘0.0.5’

## STEP 1: Munge the summary statistics 
# CAD Nikpay: maf present, info present 
# T2D Mahajan: maf not present, info not present
# MD with UKBB: maf present, info present (treating FRQ_A_116209 as EAF) 

munge(c("cad.add.160614.website.txt", "Mahajan.NatGenet2018b.T2D.European_mapped_CHR_ALL.txt", "PGC_UKB_23andMe_depression_genome-wide.txt.gz"), "w_hm3.snplist",trait.names=c("CAD_Nik","T2D_M", "MD_UKBB"), c(184305, 898130, 807553),info.filter = 0.9, maf.filter = 0.01) 
# CAD: maf present, info present 
# Mahajan: maf yes, info no
# DEP full: maf yes (called freq1), info no

traits <- c("CAD_Nik.sumstats.gz", "T2D_M.sumstats.gz", "MD_UKBB.sumstats.gz")
sample.prev <- c(.33,.09,.44)
population.prev <- c(.07,.10,.15)  # depression prevalence taken from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3131101/ 
ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"
trait.names<-c("CAD_Nik", "T2D_M", "MD_UKBB")
LDSCoutput_UKBB <- ldsc(traits=traits, sample.prev=sample.prev, population.prev=population.prev, ld=ld, wld=wld, trait.names=trait.names)

##Computing the correlation matrix
round(cov2cor(LDSCoutput_UKBB$S),2)

### Common factor model
commonfactor(covstruc = LDSCoutput_UKBB, estimation="DWLS")

files=c("cad.add.160614.website.txt", "Mahajan.NatGenet2018b.T2D.European_mapped_CHR_ALL.txt", "PGC_UKB_23andMe_depression_genome-wide.txt.gz")
ref= "reference.1000G.maf.0.005.txt"
trait.names=c("CAD_Nik", "T2D_M", "MD_UKBB")
se.logit=c(T,T,T) 
# CAD Nikpay - has logistic regression and 2*pnorm(abs(.013006)/.017324, lower.tail = FALSE) gives identical p-value to sumstats (0.4528) 
# Mahajan - 2*pnorm(abs(-0.098)/0.19, lower.tail = FALSE) gives identical p-value to sumstats (0.61) 
# MD_UKBB - formula on tutorial gives very similar (although not identical results), although meta-analysis section on paper mentions "using the log of the odds ratios and the standard errors of the log of the odds ratio"
info.filter=0.6
maf.filter=0.01

m_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=NULL,linprob=NULL,N=c(184305, 898130, 807553),betas=NULL,info.filter=info.filter,maf.filter=maf.filter,keep.indel=FALSE, parallel=T)

save(m_sumstats, file = "m_sumstats_CADnikpay.T2Dmaha.MDwithUKBB_revision_FULL.RData") # m_sumstats has 6,820,149 SNPs
#load("m_sumstats_CADnikpay.T2Dmaha.MDwithUKBB_revision_FULL.RData")

model<-"F1=~CAD_Nik + T2D_M + MD_UKBB
CAD_Nik~~a*CAD_Nik
T2D_M~~b*T2D_M
a > .001
b > .001
F1 ~ SNP"

# running in parallel using 11 cores and no GC as LD intercepts for all sumstats < 1
result<-userGWAS(covstruc = LDSCoutput_UKBB, SNPs = m_sumstats, estimation = "DWLS", model = model, sub=c("F1~SNP"), cores = 42, toler = 1e-50, SNPSE = FALSE, parallel = T,GC="none",MPI=FALSE,smooth_check=FALSE)
# should be 6820149 rows 

result <- as.data.frame(result)
# save summary statistics with all SNPs 
write.table(result, "multimorbidity_NikMahUKBB_noGC_revision_FULL.txt", row.names=FALSE, sep="\t", quote=FALSE)

library(tidyverse)
warning <- result %>% filter(warning != 0) 
# no warnings 


### caclculate effective sample size ###

#restrict to MAF of 40% and 10%
result2<-subset(result, result$MAF <= .4 & result$MAF >= .1)
#calculate expected sample size (N_hat)
N_hat<-mean(1/((2*result2$MAF*(1-result2$MAF))*result2$SE^2))
N_hat
# 562506.7 
# sample size for Nikpay, Mahajan and Howard full UKBB is 562,507 individuals

sig <- result %>% filter(Pval_Estimate < .00000005) # 389 SNPs under this threshold 








