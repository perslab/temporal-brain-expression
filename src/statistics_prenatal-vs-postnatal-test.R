############### SYNOPSIS ###################
# This script does statistical testing on the RNAseq data set
# The script performs several different tests of differences in expression levels btw prenatal and postnatal developmental stages.
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools)

rm(list=ls())

#===================================== GET CONFIG =====================================#
#======================================================================================#
source("CONFIG.R", echo=TRUE) 
#======================================================================================#



############################# LOAD EXPRESSION DATA #################################
file.out.processed <- file.path(path.main.analysis, "RData/data.temporal-brain-expression.rnaseq_expression_processed.RData")
load(file.out.processed) # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt

str(df.expression_matrix.clean.melt)

########### Setting prioritization #########
file.gene_prioritization <- "../gene_lists/gene_prioritization.txt"


gene_list <- basename(file_path_sans_ext(file.gene_prioritization))


########### READ prioritization file ###########
df.gene_prioritization <- read.csv(file.gene_prioritization,h=T)
## Check mapping - compare to "df.expression_matrix.clean" NOT MOLTEN, will be too slow
sum(df.expression_matrix.clean$ensembl_gene_id %in% df.gene_prioritization[,1]) # mapped genes
sum(!df.gene_prioritization[,1] %in% df.expression_matrix.clean$ensembl_gene_id) # non-mapped genes
cat(as.character(df.gene_prioritization[!df.gene_prioritization[,1] %in% df.expression_matrix.clean$ensembl_gene_id,]), sep="\n") # print non-mapped genes

### Setting priorizied factor for MOLTEN df
df.expression_matrix.clean.melt$gene_type <- as.factor(ifelse(df.expression_matrix.clean.melt$ensembl_gene_id %in% df.gene_prioritization[,1], "prioritized", "other"))
table(df.expression_matrix.clean.melt$gene_type)

str(df.expression_matrix.clean.melt)


################## PREPARING FOR TEST ###############
## subsetting on prioritized/associated genes
df.expression_matrix.clean.melt.priori <- subset(df.expression_matrix.clean.melt, gene_type=="prioritized")

## calculating mean for each gene in prenatal and postnatal
df.natal.gene.mean <- ddply(df.expression_matrix.clean.melt.priori, .(ensembl_gene_id, natal), summarise,
                            mean=mean(value, na.rm=T))


################## T-test: prenatal vs postnatal | *PAIRED* ###############

### t.test
levels(df.natal.gene.mean$natal)
fit1 <- t.test(mean~natal, data=df.natal.gene.mean, alternative="less", paired=TRUE) # OBS: EAv2 is testing "less" (SCZ was testing alternative="greater")
fit1

fit1$p.value # --> 4.223849e-08

### calculation odd's ratio
df.or <- ddply(df.natal.gene.mean, .(natal), summarise,
               group_mean=mean(mean, na.rm=T))
df.or
or_log2_corrected <- (2^df.or[df.or$natal=="prenatal","group_mean"]-1)/(2^df.or[df.or$natal=="postnatal","group_mean"]-1)
#formula: or = ((2^mu1)-1)/((2^mu2)-1)
or_log2_corrected # --> 1.356436



######################################################################################
###################################### PLOTS #########################################

### plotting distribution of mean values
q1 <- ggplot(df.natal.gene.mean, aes(x=mean, fill=natal)) + geom_histogram(aes(y=..density..)) +geom_density(alpha=0.8)
q1
### plotting box plot
q2 <- ggplot(df.natal.gene.mean, aes(x=natal, y=mean, fill=natal)) + geom_boxplot()
q2

### The test above is a gene-level test. However, if one does not believe that the DEPICT-prioritized genes
### can reasonably be taken as a random sample of all causal genes, then one should perform a donor-level
### test. That is, we ask the question: "Fix the current list of prioritized genes. We observe that the expression
### of these genes in the brain is higher during certain developmental stages. Is this observation likely
### to hold up if the BrainSpain Developmental Transcriptome ascertainment procedure secures more donors?"

######################################################################################
###################### DONOR-LEVEL STATISTICAL TESTS #################################

# Run all data-loading, dataframe-creation, and variable-definition steps in
# graphics_genes_temporal_trajectories.R
# up till the first intance of pooledSD (after the creation of df.summary.stage).

### T-TEST OF A CONTRAST
### We test whether the linear combination of temporal-stage means
###   (1/6)*stage01 + (1/6)*stage02 + (1/6)*stage03 + (1/6)*stage04 + (1/6)*stage05 + (1/6)*stage06
### - (1/6)*stage07 - (1/6)*stage08 - (1/6)*stage09 - (1/6)*stage10 - (1/6)*stage11 - (1/6)*stage12
### is significantly different from zero.

N <- table(df.summary.ind$stage)/3 # divide by 3 because each individual appears 3 times (once per gene list)
factor <- 0
for(i in 1:12) factor <- factor + (1/(36*N[i]))
SEg <- pooledSD*sqrt( factor )

pre <- mean(subset(df.summary.stage,stage %in% c("s2a","s2b","s3a","s3b","s4","s5") & gene_list=="gene_prioritization.txt")$mean1)
post <- mean(subset(df.summary.stage,stage %in% c("s6","s7","s8","s9","s10","s11") & gene_list=="gene_prioritization.txt")$mean1)
pt((pre - post)/SEg,df=(length(unique(df.summary.ind$donor_id)) - length(unique(df.summary.ind$stage))),lower.tail=FALSE)*2

### T-TEST OF PRE AND POST
### One can test the hypothesis that all stages have the same mean and that the pre vs post
### difference has arisen by chance.

t.pre.post <- t.test(
subset(df.summary.ind,stage %in% c("s2a","s2b","s3a","s3b","s4","s5") & gene_list=="gene_prioritization.txt")$mean,
subset(df.summary.ind,stage %in% c("s6","s7","s8","s9","s10","s11") & gene_list=="gene_prioritization.txt")$mean)

### PERMUTATION TEST OF PRE AND POST
### Same null as above, except tested non-parametrically.

permut <- function(i){
	df <- subset(df.summary.ind,gene_list=="gene_prioritization.txt")
	timestamps <- sample(df$stage,nrow(df),replace=FALSE) # permutation (reordering) of temporal-stage labels
	prenatal <- mean(subset(df,timestamps %in% c("s2a","s2b","s3a","s3b","s4","s5"))$mean)
	postnatal <- mean(subset(df,!(timestamps %in% c("s2a","s2b","s3a","s3b","s4","s5")))$mean)
	return(prenatal - postnatal) }

P <- 10000000	
permut_results <- unlist(lapply(1:P,permut))
sum(permut_results>(t.pre.post$estimate[1] - t.pre.post$estimate[2]))/P

### ROBUSTNESS CHECK OF DONOR-LEVEL TESTS
### See if the pre vs post difference holds up in brain regions represented by donors at all temporal stages.
### Sometimes a stage is represented by only one donor. Use the non-contrast tests to get around this.

df.depict <- subset(df.gene_list,gene_list=="gene_prioritization.txt")

### AMY (amygdaloid complex)

df.AMY <- subset(df.depict,structure_acronym=="AMY")

df.ind.AMY <- ddply(df.AMY, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))                   

df <- df.ind.AMY
permut <- function(i){
	timestamps <- sample(df$stage,nrow(df),replace=FALSE) # permutation (reordering) of temporal-stage labels
	prenatal <- mean(subset(df,timestamps %in% c("s2a","s2b","s3a","s3b","s4","s5"))$median)
	postnatal <- mean(subset(df,!(timestamps %in% c("s2a","s2b","s3a","s3b","s4","s5")))$median)
	return(prenatal - postnatal) }

t.pre.post <- t.test(
subset(df,stage %in% c("s2a","s2b","s3a","s3b","s4","s5"))$median,
subset(df,stage %in% c("s6","s7","s8","s9","s10","s11"))$median)
	
P <- 1000	
permut_results <- unlist(lapply(1:P,permut))
sum(permut_results>(t.pre.post$estimate[1] - t.pre.post$estimate[2]))/P	

# # HIP (hippocampus)

df.HIP <- subset(df.depict,structure_acronym=="HIP")
df.ind.HIP <- ddply(df.HIP, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))
df <- df.ind.HIP

t.pre.post <- t.test(
subset(df,stage %in% c("s2a","s2b","s3a","s3b","s4","s5"))$median,
subset(df,stage %in% c("s6","s7","s8","s9","s10","s11"))$median)
	
P <- 1000	
permut_results <- unlist(lapply(1:P,permut))
sum(permut_results>(t.pre.post$estimate[1] - t.pre.post$estimate[2]))/P	

# ITC (inferolateral temporal cortex)

df.ITC <- subset(df.depict,structure_acronym=="ITC")
df.ind.ITC <- ddply(df.ITC, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))  
df <- df.ind.ITC

t.pre.post <- t.test(
subset(df,stage %in% c("s2a","s2b","s3a","s3b","s4","s5"))$median,
subset(df,stage %in% c("s6","s7","s8","s9","s10","s11"))$median)
	
P <- 1000	
permut_results <- unlist(lapply(1:P,permut))
sum(permut_results>(t.pre.post$estimate[1] - t.pre.post$estimate[2]))/P	

# MFC (anterior cingulate cortex, AKA medial prefrontal cortex)

df.MFC <- subset(df.depict,structure_acronym=="MFC")
df.ind.MFC <- ddply(df.MFC, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))   
df <- df.ind.MFC

t.pre.post <- t.test(
subset(df,stage %in% c("s2a","s2b","s3a","s3b","s4","s5"))$median,
subset(df,stage %in% c("s6","s7","s8","s9","s10","s11"))$median)
	
P <- 1000	
permut_results <- unlist(lapply(1:P,permut))
sum(permut_results>(t.pre.post$estimate[1] - t.pre.post$estimate[2]))/P	

# OFC (orbital prefrontal cortex)

df.OFC <- subset(df.depict,structure_acronym=="OFC")
df.ind.OFC <- ddply(df.OFC, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))    
df <- df.ind.OFC

t.pre.post <- t.test(
subset(df,stage %in% c("s2a","s2b","s3a","s3b","s4","s5"))$median,
subset(df,stage %in% c("s6","s7","s8","s9","s10","s11"))$median)
	
P <- 1000	
permut_results <- unlist(lapply(1:P,permut))
sum(permut_results>(t.pre.post$estimate[1] - t.pre.post$estimate[2]))/P	

# VFC (ventrolateral prefrontal cortex)

df.VFC <- subset(df.depict,structure_acronym=="VFC")
df.ind.VFC <- ddply(df.VFC, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE)) 
df <- df.ind.VFC

t.pre.post <- t.test(
subset(df,stage %in% c("s2a","s2b","s3a","s3b","s4","s5"))$median,
subset(df,stage %in% c("s6","s7","s8","s9","s10","s11"))$median)
	
P <- 1000	
permut_results <- unlist(lapply(1:P,permut))
sum(permut_results>(t.pre.post$estimate[1] - t.pre.post$estimate[2]))/P	