############### SYNOPSIS ###################
### DESCRIPTION
# The script loads PROCESSED RNAseq data.
# It was created to be able to handle multiple gene list for plotting. The script uses plyr to load in multiple gene lists.
# The script was used to generate PUBLICATION READY GRAPHICS

### REMARKS
# THIS VERSION TAKES THE MEDIAN MEAN GENE EXPRESSION VALUES.
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



########################################### LOAD NULL data ###################################
################################## ** ASSOCIATED genes ** ##########################
### rnaseq
load(file.path(path.main.analysis, "RData/null_RData_broad_rnaseq_associated_priority.RData")) #time_elapsed, list.par_analysis

############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.median.summary")
############## Combining summary data frames
df.null.mean.summary.assoc <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary.assoc <- ldply(list.null.median.summary) # COMBINING list of data frames
# THIS BELOW APROACH WAS USED TO GENERATE PLOT FOR *HIRCHHORN LAB PRESENTATION*
# df.null.median.summary.sem.assoc <- ddply(df.null.median.summary.assoc, c("stage"), summarise,
#                                          mean1 = mean(mean, na.rm=TRUE),
#                                          sd1   = sd(mean, na.rm=TRUE))
####
# First we find the stage average across structures FOR EACH PERMUTATION. 
# In this way we process each permutation as we did for SCZ prioritized genes
df.null.median.summary.per_perm.assoc <- ddply(df.null.median.summary.assoc, .(permutation, stage), summarise,
                                               mean_stage_across_structure = mean(mean, na.rm=TRUE),
                                               sd_stage_across_structure   = sd(mean, na.rm=TRUE))

df.null.median.summary.sem.assoc <- ddply(df.null.median.summary.per_perm.assoc, .(stage), summarise,
                                          mean1 = mean(mean_stage_across_structure, na.rm=TRUE),
                                          sd1   = sd(mean_stage_across_structure, na.rm=TRUE))

################################ ** PRIORITIZED genes ** #########################
### rnaseq
load(file.path(path.main.analysis, "RData/null_RData_broad_rnaseq_prioritized_priority.RData")) #time_elapsed, list.par_analysis

############### Extracting from list
list.null.mean.summary <- lapply(list.par_analysis, "[[", "df.null.mean.summary")
list.null.median.summary <- lapply(list.par_analysis, "[[", "df.null.median.summary")
############## Combining summary data frames
df.null.mean.summary.prio <- ldply(list.null.mean.summary) # COMBINING list of data frames
df.null.median.summary.prio <- ldply(list.null.median.summary) # COMBINING list of data frames
df.null.median.summary.per_perm.prio <- ddply(df.null.median.summary.prio, .(permutation, stage), summarise,
                                              mean_stage_across_structure = mean(mean, na.rm=TRUE),
                                              sd_stage_across_structure   = sd(mean, na.rm=TRUE))

df.null.median.summary.sem.prio <- ddply(df.null.median.summary.per_perm.prio, .(stage), summarise,
                                         mean1 = mean(mean_stage_across_structure, na.rm=TRUE),
                                         sd1   = sd(mean_stage_across_structure, na.rm=TRUE))



############################# READING GENE LISTs #################################
path.datafiles <- "../gene_lists"

###### Read SPECIFIC FILES:
filenames2read <- c("gene_associated.txt", "gene_nearest.txt", "gene_prioritization.txt")

files <- as.list(paste(path.datafiles, filenames2read, sep="/"))
names(files) <- filenames2read
files
list_of_data <- llply(files, read.csv)#row.names = 1 --> NO!, stringsAsFactors = FALSE

names(list_of_data)

extract_genes_from_molten_df <- function(df_gene_list) {
  print("done")
  df <- subset(df.expression_matrix.clean.melt, ensembl_gene_id %in% df_gene_list[,1])
}
df.gene_list <- ldply(list_of_data, extract_genes_from_molten_df, .id="gene_list")
## Converting .id=gene_list to factor
df.gene_list$gene_list <- as.factor(df.gene_list$gene_list) 
str(df.gene_list)
levels(df.gene_list$gene_list)


###################################### PROCESSING GENE lists ################################
###### Mean per stage/structure ######## 
df.summary <- ddply(df.gene_list, c("stage", "structure_acronym", "gene_list"), summarise,
                    mean = median(value, na.rm=TRUE),
                    sd   = sd(value, na.rm=TRUE))
## plyr magic for renaming factor level
levels(df.summary$gene_list)
df.summary$gene_list <- revalue(df.summary$gene_list, c("gene_associated.txt"="Associated Genes", "gene_nearest.txt"="Nearest Genes", "gene_prioritization.txt"="Prioritized Genes"))
levels(df.summary$gene_list)

###### Mean per stage - FINAL ##########
df.summary.sem <- ddply(df.summary, c("stage","gene_list"), summarise,
                        mean1 = mean(mean, na.rm=TRUE),
                        sd1   = sd(mean, na.rm=TRUE))


###################################### Calculating overall mean ################################
### *** Runtime ~ 10 s ***
df.all.sem <- ddply(ddply(df.expression_matrix.clean.melt, .(stage, structure_acronym), summarise, mean=median(value, na.rm=TRUE)), .(stage), summarise, mean1=mean(mean, na.rm=TRUE),  sd1=sd(mean, na.rm=TRUE))




#############################################################################################
################################### SOURCE STEP-by-STEP plot ################################
#############################################################################################
source("function_stepbystep_null_plots.R")
#############################################################################################
#############################################################################################


#### show all plots
list_of_alphas <- lapply(initialyze_alpha(), function(x) {x <- 1})
p <- plotMe(list_of_alphas)
p


###################################### STEP 1 - assoc_mean ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1

p <- plotMe(list_of_alphas)
p
p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,0,0,0), size=c(1,1,3,1,3)))); p
saveMe(p, "plot1_null")



###################################### STEP 2 - assoc_mean + null ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["assoc_null_line"]] <- 1
list_of_alphas[["assoc_null_rib"]] <- 0.3

p <- plotMe(list_of_alphas)
p
p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,0.1,0,0), size=c(1,1,3,1,3)))); p
saveMe(p, "plot2_null")

###################################### STEP 3 - prio_mean + null ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_null_line"]] <- 1
list_of_alphas[["prio_null_rib"]] <- 0.3

p <- plotMe(list_of_alphas)
p
p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,0,0,1,0.1), size=c(1,1,3,1,3)))); p
saveMe(p, "plot3_null")

###################################### STEP 4 - all ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["assoc_null_line"]] <- 1
list_of_alphas[["assoc_null_rib"]] <- 0.3
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_null_line"]] <- 1
list_of_alphas[["prio_null_rib"]] <- 0.3

p <- plotMe(list_of_alphas)
p
p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(alpha=c(1,1,0.1,1,0.1), size=c(1,1,3,1,3)))); p
saveMe(p, "plot4_null")



