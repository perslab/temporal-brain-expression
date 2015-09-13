library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())
wd <- "/Users/pascaltimshel/Dropbox/0_Projects/p_EA/git/temporal-brain-expression/src"
setwd(wd)


#===================================== GET CONFIG =====================================#
#======================================================================================#
source("CONFIG.R", echo=TRUE) 
#======================================================================================#



########################################### LOAD data ###################################
############ ** ASSOCIATED genes ** ##########
### rnaseq
#load("RData/null_RData_broad_rnaseq_associated_paired-ttest_priority_ish.RData") #time_elapsed, list.par_analysis

############ ** PRIORITIZED genes ** ##########
### rnaseq
load(file.path(path.main.analysis,"RData/null_RData_broad_rnaseq_prioritized_paired-ttest_priority_ish.RData")) #time_elapsed, list.par_analysis



############################# EXTRACT BROAD DATA #################################
############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.median.summary")
list.null.fits <- lapply(list.par_analysis, "[[", "list.null.fits")
### Generating data.frames - USING ldply!
df.null.mapping <- ldply(list.par_analysis, "[[", "df.null.mapping") # the following worked when "scalar variables" were saved in the par.analyze_null_genes list: df.null.mapping <- ldply(list.par_analysis, function(x) {data.frame(x[["n_mapped_genes"]], x[["n_unmapped_genes"]])}) 
############## Subsequent extractions
####### Extracting fits
list.null.fit.natal.paired <- lapply(list.null.fits, "[[", "fit.natal.paired") # *NEW PAIRED* 
list.null.fit.natal <- lapply(list.null.fits, "[[", "fit.natal")
list.null.fit.prioritized.higher <- lapply(list.null.fits, "[[", "fit.prioritized.higher")
list.null.list.fit.stage <- lapply(list.null.fits, "[[", "list.fit.stage")
############## Extracting from list.null.list.fit.stage
df.fit.stage<-ldply(seq_along(list.null.list.fit.stage), function(i) {ldply(list.null.list.fit.stage[[i]], function(fit.stage) {data.frame(t_statistic=fit.stage$statistic, perm=i)})})

############## Combining summary data frames
df.null.mean.summary <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary <- ldply(list.null.median.summary) # COMBINING list of data frames


############################# STATISTCAL TESTs #################################


################### PRIORITIZED genes #############
### *** REMEMBER to update the list.null.fit.natal before running the ttest

# EA RELEASE 5 | ONE SIDED TEST
#t = -5.6442, df = 145, p-value = 4.224e-08


### Empirical p-value: NATAL TEST
obs.t_statistic <- -5.6442 # EA associated genes, RNAseq
df.null.t_statistic <- ldply(list.null.fit.natal.paired, function(x) {t=x$statistic}, .id="id") # same as sapply(list.null.natal_fits, "[[", c("statistic"))
df.null.t_statistic
empirical.pvalue <- sum(abs(df.null.t_statistic$t) > abs(obs.t_statistic))/length(df.null.t_statistic$t) # calc empirical p-val, one sided, alternative="greater"
empirical.pvalue
p <- ggplot(data=df.null.t_statistic, aes(x=t)) + geom_histogram(binwidth=0.2)
p <- p + geom_vline(xintercept=c(obs.t_statistic), linetype="dotted", size=1)
p <- p + labs(title=paste("Null GWAS - prioritized genes. pval=", round(empirical.pvalue,3), sep=""), y="Number of tests", x="t-statistic")
p

ggsave("EAv2_empirical_distribution_rnaseq_prioritized_genes_prenatal-vs-postnatal-9x4.pdf", w=9, h=4)

