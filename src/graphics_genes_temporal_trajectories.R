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
################################# PLOT - 1 SD ERROR BARS confidence #########################
#############################################################################################

########### PLOT IT! ###########

p <- ggplot()
p <- p + geom_line(data=subset(df.summary, gene_list == "Prioritized Genes"), aes(x=stage, y=mean, group=structure_acronym, color="Prioritized genes (structures)")) #linetype="Brain regions"
p
### Adding mean Prioritized
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='#d7191c', width=0.2)
p
### Adding mean ALL (df.all.sem)
p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1),color='black', width=0.2)
p
### Adding Nearest Genes
#p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, y=mean1, group=1, color="Nearest genes"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='sky blue', width=0.2)
#p
### Adding Associated Genes
#p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, y=mean1, group=1, color="Associated genes"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, ymax=mean1+sd1, ymin=mean1-sd1), color='orange', width=0.2)
#p


p <- p + scale_color_manual(name="Gene list", values=c(
                                "Prioritized genes (structures)"="gray", 
                                "Prioritized genes"="#d7191c",
                                "All genes"="black",
                                "Nearest genes"="sky blue",
                                "Associated genes"="orange",
                                guide='legend'))
p

###### Adding vertical line - prenatal vs. postnatal
p <- p + geom_vline(xintercept=6.5, color="black", linetype="dashed")
p

######### Adding x-tickmarks for stage
stage_converter <- c("s1"="Embryonic",
                     "s2a"="Early prenatal",
                     "s2b"="Early prenatal",
                     "s3a"="Early mid-prenatal",
                     "s3b"="Early mid-prenatal",
                     "s4"="Late mid-prenatal",
                     "s5"="Late prenatal",
                     "s6"="Early infancy",
                     "s7"="Late infancy",
                     "s8"="Early childhood",
                     "s9"="Late childhood",
                     "s10"="Adolescence",
                     "s11"="Adulthood")
p <- p + scale_x_discrete(name="", labels = stage_converter) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=rel(1.15)))
p


### VARIABLE
label.y <- expression(paste("Mean brain expression ", bgroup("(",log[2]~group("[",RPKM+1,"]"),")") ))
p <- p + labs(y=label.y)
p


### SAVE PLOT
filename.plot <- file.path(path.main.analysis, "plot_gene_temporal_trajectories_rnaseq-9x6.pdf")
filename.plot
ggsave(file=filename.plot, w=9, h=6)







#############################################################################################
###################################### PLOT  2- 95% confidence ##############################
# USED FOR PUBLICATION
#############################################################################################

### Option 1: use 95 % confidence interval
# a) using a normal distribution | qnorm(0.975) --> ~1.96
# b) using a t-distribution | qt(0.975,df=26-1)=2.059539
### Option 2: Standard error of mean
# a) SEM = SE/sqrt(n) = SE/sqrt(26)
# SEE also http://en.wikipedia.org/wiki/Standard_error#Assumptions_and_usage
### CONCLUSION
# We use Option 1a)
sem.factor <- qnorm(0.975) # 95% confidence interval (plus/minus)

p <- ggplot()
p <- p + geom_line(data=subset(df.summary, gene_list == "Prioritized Genes"), aes(x=stage, y=mean, group=structure_acronym, color="Prioritized genes (structures)")) #linetype="Brain regions"
p
### Adding mean Prioritized
p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sem.factor*sd1, ymin=mean1-sem.factor*sd1),color='#d7191c', width=0.2)
p
### Adding mean ALL (df.all.sem)
p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sem.factor*sd1, ymin=mean1-sem.factor*sd1),color='black', width=0.2)
p
### Adding Nearest Genes
#p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, y=mean1, group=1, color="Nearest genes"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, ymax=mean1+sem.factor*sd1, ymin=mean1-sem.factor*sd1), color='sky blue', width=0.2)
#p
### Adding Associated Genes
#p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, y=mean1, group=1, color="Associated genes"), linetype='solid', size=1)
#p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, ymax=mean1+sem.factor*sd1, ymin=mean1-sem.factor*sd1), color='orange', width=0.2)
#p


p <- p + scale_color_manual(name="Gene list", values=c(
  "Prioritized genes (structures)"="gray", 
  "Prioritized genes"="#d7191c",
  "All genes"="black",
  "Nearest genes"="sky blue",
  "Associated genes"="orange",
  guide='legend'))
p

###### Adding vertical line - prenatal vs. postnatal
p <- p + geom_vline(xintercept=6.5, color="black", linetype="dashed")
p

######### Adding x-tickmarks for stage
stage_converter <- c("s1"="Embryonic",
                     "s2a"="Early prenatal",
                     "s2b"="Early prenatal",
                     "s3a"="Early mid-prenatal",
                     "s3b"="Early mid-prenatal",
                     "s4"="Late mid-prenatal",
                     "s5"="Late prenatal",
                     "s6"="Early infancy",
                     "s7"="Late infancy",
                     "s8"="Early childhood",
                     "s9"="Late childhood",
                     "s10"="Adolescence",
                     "s11"="Adulthood")
p <- p + scale_x_discrete(name="", labels = stage_converter) + theme(axis.text.x = element_text(angle = 35, hjust = 1, size=rel(1.15)))
p


### VARIABLE
label.y <- expression(paste("Mean brain expression ", bgroup("(",log[2]~group("[",RPKM+1,"]"),")") ))
p <- p + labs(y=label.y)
p


### SAVE PLOT
filename.plot <- file.path(path.main.analysis, "plot_gene_temporal_trajectories_rnaseq-95confint-9x6.pdf")
filename.plot
ggsave(file=filename.plot, w=9, h=6)




