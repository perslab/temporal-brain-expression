############### SYNOPSIS ###################
# This script does statistical testing on the RNAseq data set
# The script performs a paired t-test to test for differences in expression levels for prenatal and postnatal developmental stages.
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


######################################################################################
################################# STATISTICAL TESTs ##################################
######################################################################################

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

