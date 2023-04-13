
# A multivariate GWAS of psycho-cardiometabolic multimorbidity - version WITHOUT UKBB
# Genomic SEM tutorial available at https://github.com/GenomicSEM/GenomicSEM/wiki/5.-User-Specified-Models-with-SNP-Effects 

# install and load packages
#install.packages("devtools")
library(devtools)
#install_github("GenomicSEM/GenomicSEM")
require(GenomicSEM)
packageVersion("GenomicSEM")

## STEP 1: Munge the summary statistics 

munge(c("cad.add.160614.website.txt", "METAANALYSIS_DIAGRAM_SE1_mapped.txt", "daner_pgc_mdd_meta_w2_rmUKBB_full.gz"), "w_hm3.snplist",trait.names=c("CAD","T2D", "MD"), c(184305, 159208, 518370),info.filter = 0.9, maf.filter = 0.01) 
# CAD: maf present, info present 
# T2D: maf no, info no
# MD: maf present, info present (treating FRQ_A_116209 as the EAF in daner_pgc_mdd_meta_w2_rmUKBB_full.gz)

#Step 2: Run multivariable LDSC 

traits <- c("CAD.sumstats.gz", "T2D.sumstats.gz", "MD.sumstats.gz")
sample.prev <- c(.33,.20,.27)
population.prev <- c(.07,.10,.15)  #depression prevalence taken from: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3131101/
ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"
trait.names<-c("CAD","T2D","MD")
LDSCoutput <- ldsc(traits=traits, sample.prev=sample.prev, population.prev=population.prev, ld=ld, wld=wld, trait.names=trait.names)

##Computing the correlation matrix
round(cov2cor(LDSCoutput$S),2)

### Common factor model
commonfactor(covstruc = LDSCoutput, estimation="DWLS")


files=c("cad.add.160614.website.txt", "METAANALYSIS_DIAGRAM_SE1_mapped.txt", "daner_pgc_mdd_meta_w2_rmUKBB_full.gz")
ref= "reference.1000G.maf.0.005.txt"
trait.names=c("CAD", "T2D", "MD")
se.logit=c(T,T,T)  # confirmed that daner SE is SE of log(OR), so se.logit = T 
# CAD Nikpay - has logistic regression and 2*pnorm(abs(.013006)/.017324, lower.tail = FALSE) gives identical p-value to sumstats (0.4528) 
# same for Scot 2*pnorm(abs(0.0038)/0.016, lower.tail = FALSE) = 0.8122689
info.filter=0.6
maf.filter=0.01

m_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=NULL,linprob=NULL,N=c(184305, 159208, 518370),betas=NULL,info.filter=info.filter,maf.filter=maf.filter,keep.indel=FALSE, parallel=F)
# 7323398 rows

save(m_sumstats, file = "m_sumstats_CADnikpay.Scott.Daner_revision.RData") # m_sumstats has 7323398 SNPs   

# specifying common factor model and constraining residual variance to be positive
model<-"F1=~CAD + T2D + MD
F1 ~ SNP
CAD~~a*CAD
T2D~~b*T2D
a > .001
b > .001"

# running in parallel using 11 cores and no GC as LD intercepts for all sumstats ~ 1
result<-userGWAS(covstruc = LDSCoutput, SNPs = m_sumstats, estimation = "DWLS", model = model, sub=c("F1~SNP"), cores = 11, toler = 1e-50, SNPSE = FALSE, parallel = T,GC="none",MPI=FALSE,smooth_check=FALSE)

result <- as.data.frame(result)
# explore result
dim(result)
head(result)

# check for warnings
library(tidyverse)

warning <- result %>% filter(warning != 0)
errors <- result %>% filter(error != 0)
# no errors or warnings

# save summary statistics  with all SNPs 
write.table(result, "multimorbidity_Nikpay.Scott.Daner_noGC_revision.txt", row.names=FALSE, sep="\t", quote=FALSE)

# calculate effective sample size
#restrict to MAF of 40% and 10%
result2<-subset(result, result$MAF <= .4 & result$MAF >= .1)
#calculate expected sample size (N_hat)
N_hat<-mean(1/((2*result2$MAF*(1-result2$MAF))*result2$SE^2))
#  156717.4 sample size 

sig <- result %>% filter(Pval_Estimate < .00000005) # 1,312 SNPs under this threshold 






