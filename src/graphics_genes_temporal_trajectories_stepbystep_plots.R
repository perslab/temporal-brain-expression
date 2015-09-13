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
source("function_stepbystep_plots.R")
#############################################################################################
#############################################################################################


#### show all plots
list_of_alphas <- lapply(initialyze_alpha(), function(x) {x <- 1})
p <- plotMe(list_of_alphas)
p


###################################### STEP 0 ################################
list_of_alphas <- initialyze_alpha()
p <- plotMe(list_of_alphas, plot_vline=FALSE)
p
saveMe(p, "plot0")

###################################### STEP 1 ################################
list_of_alphas <- initialyze_alpha()
p <- plotMe(list_of_alphas, plot_vline=TRUE)
p
saveMe(p, "plot1")

###################################### STEP 2 - all ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1

p <- plotMe(list_of_alphas)
p
p <- p + guides(colour = guide_legend(override.aes = list(alpha=c(1,0,0,0,0), size=c(1,1,1,1,0.1)))); p
p
saveMe(p, "plot2")

###################################### STEP 3 - associated ################################
list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1

p <- plotMe(list_of_alphas)
p
p <- p + guides(colour = guide_legend(override.aes = list(alpha=c(1,1,0,0,0), size=c(1,1,1,1,0.1)))); p
p
saveMe(p, "plot3")

###################################### STEP 4 - Nereast ################################

list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["nearest"]] <- 1
p <- plotMe(list_of_alphas)
p
p <- p + guides(colour = guide_legend(override.aes = list(alpha=c(1,1,1,0,0), size=c(1,1,1,1,0.1)))); p
p
saveMe(p, "plot4")


###################################### STEP 5 - Prioritized mean ################################

list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
list_of_alphas[["assoc_mean"]] <- 1
list_of_alphas[["nearest"]] <- 1
list_of_alphas[["prio_mean"]] <- 1

p <- plotMe(list_of_alphas)
p
p <- p + guides(colour = guide_legend(override.aes = list(alpha=c(1,1,1,1,0), size=c(1,1,1,1,0.1)))); p
p
saveMe(p, "plot5")


###################################### STEP 6 - Prioritized mean + structures ################################

list_of_alphas <- initialyze_alpha()
list_of_alphas[["all"]] <- 1
# list_of_alphas[["assoc_mean"]] <- 1
# list_of_alphas[["nearest"]] <- 1
list_of_alphas[["prio_mean"]] <- 1
list_of_alphas[["prio_struc"]] <- 1
p <- plotMe(list_of_alphas)
p
p <- p + guides(colour = guide_legend(override.aes = list(alpha=c(1,0,0,1,1), size=c(1,1,1,1,0.1)))); p
p
saveMe(p, "plot6")

