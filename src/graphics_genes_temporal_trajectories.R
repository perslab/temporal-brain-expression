############### SYNOPSIS ###################
### DESCRIPTION
# The script loads PROCESSED RNAseq data.
# It was created to be able to handle multiple gene list for plotting. The script uses plyr to load in multiple gene lists.
# The script was used to generate PUBLICATION READY GRAPHICS

### REMARKS
# TO COMPUTE THE MAIN TRAJECTORY, THIS VERSION TAKES THE MEDIAN EXPRESSION OF ALL LISTED GENES, AS MEASURED
# IN A BRAIN REGION FROM A GIVEN DONOR. THE MEAN OF THESE MEDIANS OVER BRAIN RGIONS IS THEN COMPUTED FOR
# EACH DONOR. THIS DONOR-LEVEL QUANTITY IS TREATED AS THE UNIT OF OBSERVATION.
############################################


library(plyr)
library(ggplot2)
library(reshape2)
library(tools)
library(car)

rm(list=ls())

#===================================== GET CONFIG =====================================#
#======================================================================================#
source("CONFIG.R", echo=TRUE) 
#======================================================================================#



############################# LOAD EXPRESSION DATA #################################
file.out.processed <- file.path(path.main.analysis, "RData/data.temporal-brain-expression.rnaseq_expression_processed.RData")
load(file.out.processed) # RNAseq AFTER PROCESSING | df.expression_matrix.clean, df.expression_matrix.clean.melt

str(df.expression_matrix.clean.melt)
# Each row of this dataframe corresponds to a measurement of transcript abundance for a
# given gene, brain region, and individual (who died at a certain temporal stage).
# The column "value" gives the measurement itself. The minimum value is zero and the median is
# 0.04, which means that the majority of measurements are very small.

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
# This dataframe is the same as the original expression data, except only containing genes in one
# of the three lists (genes in loci, nearest genes, DEPICT-prioritized genes). There is also an
# additional column giving the gene list in which the gene appears.

###################################### PROCESSING GENE lists ################################
###### Mean per stage/structure ######## 
df.summary <- ddply(df.gene_list, c("stage", "structure_acronym", "gene_list"), summarise,
                                   median = median(value, na.rm=TRUE),
                                   sd   = sd(value, na.rm=TRUE))
# This function call groups the data into cells defined by temporal stage, brain area, and gene list.
# It then calculates the median and SD of each cell. 
# This dataframe is used to plot the regional trajectories.
# The nesting of donated brain regions within individuals is neglected. If the only plotted feature is
# the central tendency, then this neglect scarcely affects the results.

df.summary.reg <- ddply(df.gene_list, c("stage","donor_id", "structure_acronym", "gene_list"), summarise,
                                   median = median(value, na.rm=TRUE),
                                   sd     = sd(value, na.rm=TRUE))
# Each row of this dataframe contains the median expression level of all listed genes (e.g.,
# DEPICT-prioritized genes) in a given brain region taken from a given individual (who died at a certain temporal stage).    

df.summary.ind <- ddply(df.summary.reg, c("stage","donor_id","gene_list"), summarise,
                        mean = mean(median, na.rm=TRUE),
                        sd   = sd(median, na.rm=TRUE))
# Each row of this dataframe gives the mean over brain regions of a given individual.

df.summary.stage <- ddply(df.summary.ind,c("stage","gene_list"),summarise,
                          mean1 = mean(mean,na.rm=TRUE),
                          sd1   = sd(mean,na.rm=TRUE))
                          
# Some of the temporal stages are represented by only 2 individuals. To construct reasonably narrow
# confidence intervals, we need to calculate a pooled SD if this feasible.
bartlett.test(mean~stage,subset(df.summary.ind,gene_list=="gene_prioritization.txt"))
leveneTest(mean~stage,subset(df.summary.ind,gene_list=="gene_prioritization.txt"))
# If both of these tests yield P > 0.05, then it is probably reasonable to proceed with a pooled SD.

N <- as.data.frame(table(df.summary.ind$stage)/3) # divide by 3 because each individual appears 3 times (once per gene list)
colnames(N) <- c("stage","N")
df.summary.stage <- merge(df.summary.stage,N,by="stage")
df.for.pooledSD <- subset(df.summary.stage,gene_list=="gene_prioritization.txt") # pick out the gene list of interest
pooledSD <- sqrt( sum(df.for.pooledSD$sd1^2 * (df.for.pooledSD$N - 1)) / (sum(df.for.pooledSD$N - 1)) )

# It is of interest to plot the trajectory of all genes in the database.
df.all.ind <- ddply(ddply(df.expression_matrix.clean.melt, .(stage, donor_id, structure_acronym), summarise,
median=median(value, na.rm=TRUE)), .(stage, donor_id), summarise, mean=mean(median, na.rm=TRUE),  sd=sd(median, na.rm=TRUE))
df.all.stage <- ddply(df.all.ind,c("stage"),summarise,mean1=mean(mean,na.rm=TRUE),sd1=sd(mean,na.rm=TRUE))

bartlett.test(df.all.ind$mean~df.all.ind$stage)
leveneTest(df.all.ind$mean~df.all.ind$stage)
# Both P > 0.05.

N <- as.data.frame(table(df.all.ind$stage))
colnames(N) <- c("stage","N")
df.all.stage <- merge(df.all.stage,N,by="stage")
pooledSDall <- sqrt( sum(df.all.stage$sd1^2 * (df.all.stage$N - 1)) / (sum(df.all.stage$N - 1)) )

########### PLOT IT! ###########

degfree <- length(unique(df.summary.ind$donor_id)) - length(unique(df.summary.ind$stage))

p <- ggplot()
p <- p + geom_line(data=subset(df.summary, gene_list == "gene_prioritization.txt"), 
aes(x=stage, y=median, group=structure_acronym, color="Prioritized genes (regions)")) #linetype="Brain regions"
p
### Adding mean Prioritized
p <- p + geom_line(data=subset(df.summary.stage,gene_list=="gene_prioritization.txt"),
aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=subset(df.summary.stage,gene_list=="gene_prioritization.txt"),aes(x=stage,
ymax=(mean1+(qt(.975,degfree)*pooledSD/sqrt(N))), 
ymin=(mean1-(qt(.975,degfree)*pooledSD/sqrt(N))),
color="Prioritized genes"), width=0.2)
p
### Adding mean ALL
p <- p + geom_line(data=df.all.stage, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1)
p <- p + geom_errorbar(data=df.all.stage, aes(x=stage,
ymax=(mean1+(qt(.975,degfree)*pooledSDall/sqrt(N))), 
ymin=(mean1-(qt(.975,degfree)*pooledSDall/sqrt(N))),
color="All genes"), width=0.2)
p

p <- p + scale_color_manual(name="Gene list", values=c(
  "Prioritized genes (regions)"="gray", 
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

### ROBUSTNESS CHECK
# Not all temporal stages are represented by donors of particular regions.
# If the listed genes are more highly expressed in one brain region and donors of that
# region only contribute to prenatal stages, then we have confounding.
# The regional trajectories follow the global trend and thus make such
# confounding unlikely. But let's do a more thorough check by closely inspecting
# the trajectories of just those regions represented by donors at all stages.

df.depict <- subset(df.gene_list,gene_list=="gene_prioritization.txt")
tab <- table(df.depict$stage,df.depict$structure_acronym)/146 # divide by 146 because each entry is stage/region/gene
# brain regions represented by donors at all temporal stages: AMY, HIP, ITC, MFC, OFC, VFC

# AMY (amygdaloid complex)

df.AMY <- subset(df.depict,structure_acronym=="AMY")

df.ind.AMY <- ddply(df.AMY, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))                   

df.stage.AMY <- ddply(df.ind.AMY,c("stage"),summarise,mean = mean(median,na.rm=TRUE),sd = sd(median,na.rm=TRUE))

indices <- which(is.na(df.stage.AMY$sd))
bartlett.test(median~stage,subset(df.ind.AMY,!(stage %in% df.stage.AMY$stage[indices])))
leveneTest(median~stage,subset(df.ind.AMY,!(stage %in% df.stage.AMY$stage[indices])))

N <- as.data.frame(table(df.ind.AMY$stage))
colnames(N) <- c("stage","N")
df.stage.AMY <- merge(df.stage.AMY,N,by="stage",sort=FALSE)
df.for.pooledSD <- subset(df.stage.AMY,is.na(sd)==FALSE)
pooledSD <- sqrt( sum(df.for.pooledSD$sd^2 * (df.for.pooledSD$N - 1)) / (sum(df.for.pooledSD$N - 1)) )

degfree <- length(unique(df.AMY$donor_id)) - length(unique(df.for.pooledSD$stage)) - 
( length(unique(df.AMY$stage)) - length(unique(df.for.pooledSD$stage)) )
# Subtract the extra term because the number of stages with fewer than 2 donors is equal to the number of donors not used
# to calculate the pooled SD.

p <- ggplot()
p <- p + geom_line(data=df.stage.AMY, aes(x=stage, y=mean, group=1, color="Prioritized genes"), linetype='solid', size=1)
p
p <- p + geom_errorbar(data=df.stage.AMY, aes(x=stage,
ymax=(mean+(qt(.975,degfree)*pooledSD/sqrt(N))), 
ymin=(mean-(qt(.975,degfree)*pooledSD/sqrt(N))),
color="Prioritized genes"), width=0.2)
p

# HIP (hippocampus)

df.HIP <- subset(df.depict,structure_acronym=="HIP")

df.ind.HIP <- ddply(df.HIP, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))                   

df.stage.HIP <- ddply(df.ind.HIP,c("stage"),summarise,mean = mean(median,na.rm=TRUE),sd = sd(median,na.rm=TRUE))

indices <- which(is.na(df.stage.HIP$sd))
bartlett.test(median~stage,subset(df.ind.HIP,!(stage %in% df.stage.HIP$stage[indices])))
leveneTest(median~stage,subset(df.ind.HIP,!(stage %in% df.stage.HIP$stage[indices])))

N <- as.data.frame(table(df.ind.HIP$stage))
colnames(N) <- c("stage","N")
df.stage.HIP <- merge(df.stage.HIP,N,by="stage",sort=FALSE)
df.for.pooledSD <- subset(df.stage.HIP,is.na(sd)==FALSE)
pooledSD <- sqrt( sum(df.for.pooledSD$sd^2 * (df.for.pooledSD$N - 1)) / (sum(df.for.pooledSD$N - 1)) )

degfree <- length(unique(df.HIP$donor_id)) - length(unique(df.for.pooledSD$stage)) - 
( length(unique(df.HIP$stage)) - length(unique(df.for.pooledSD$stage)) )

p <- ggplot()
p <- p + geom_line(data=df.stage.HIP, aes(x=stage, y=mean, group=1, color="Prioritized genes"), linetype='solid', size=1)
p
p <- p + geom_errorbar(data=df.stage.HIP, aes(x=stage,
ymax=(mean+(qt(.975,degfree)*pooledSD/sqrt(N))), 
ymin=(mean-(qt(.975,degfree)*pooledSD/sqrt(N))),
color="Prioritized genes"), width=0.2)
p

# ITC (inferolateral temporal cortex)

df.ITC <- subset(df.depict,structure_acronym=="ITC")

df.ind.ITC <- ddply(df.ITC, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))                   

df.stage.ITC <- ddply(df.ind.ITC,c("stage"),summarise,mean = mean(median,na.rm=TRUE),sd = sd(median,na.rm=TRUE))

indices <- which(is.na(df.stage.ITC$sd))
bartlett.test(median~stage,subset(df.ind.ITC,!(stage %in% df.stage.ITC$stage[indices])))
leveneTest(median~stage,subset(df.ind.ITC,!(stage %in% df.stage.ITC$stage[indices])))

N <- as.data.frame(table(df.ind.ITC$stage))
colnames(N) <- c("stage","N")
df.stage.ITC <- merge(df.stage.ITC,N,by="stage",sort=FALSE)
df.for.pooledSD <- subset(df.stage.ITC,is.na(sd)==FALSE)
pooledSD <- sqrt( sum(df.for.pooledSD$sd^2 * (df.for.pooledSD$N - 1)) / (sum(df.for.pooledSD$N - 1)) )

degfree <- length(unique(df.ITC$donor_id)) - length(unique(df.for.pooledSD$stage)) - 
( length(unique(df.ITC$stage)) - length(unique(df.for.pooledSD$stage)) )

p <- ggplot()
p <- p + geom_line(data=df.stage.ITC, aes(x=stage, y=mean, group=1, color="Prioritized genes"), linetype='solid', size=1)
p
p <- p + geom_errorbar(data=df.stage.ITC, aes(x=stage,
ymax=(mean+(qt(.975,degfree)*pooledSD/sqrt(N))), 
ymin=(mean-(qt(.975,degfree)*pooledSD/sqrt(N))),
color="Prioritized genes"), width=0.2)
p

# MFC (anterior cingulate cortex, AKA medial prefrontal cortex)

df.MFC <- subset(df.depict,structure_acronym=="MFC")

df.ind.MFC <- ddply(df.MFC, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))                   

df.stage.MFC <- ddply(df.ind.MFC,c("stage"),summarise,mean = mean(median,na.rm=TRUE),sd = sd(median,na.rm=TRUE))

indices <- which(is.na(df.stage.MFC$sd))
bartlett.test(median~stage,subset(df.ind.MFC,!(stage %in% df.stage.MFC$stage[indices])))
leveneTest(median~stage,subset(df.ind.MFC,!(stage %in% df.stage.MFC$stage[indices])))

N <- as.data.frame(table(df.ind.MFC$stage))
colnames(N) <- c("stage","N")
df.stage.MFC <- merge(df.stage.MFC,N,by="stage",sort=FALSE)
df.for.pooledSD <- subset(df.stage.MFC,is.na(sd)==FALSE)
pooledSD <- sqrt( sum(df.for.pooledSD$sd^2 * (df.for.pooledSD$N - 1)) / (sum(df.for.pooledSD$N - 1)) )

degfree <- length(unique(df.MFC$donor_id)) - length(unique(df.for.pooledSD$stage)) - 
( length(unique(df.MFC$stage)) - length(unique(df.for.pooledSD$stage)) )

p <- ggplot()
p <- p + geom_line(data=df.stage.MFC, aes(x=stage, y=mean, group=1, color="Prioritized genes"), linetype='solid', size=1)
p
p <- p + geom_errorbar(data=df.stage.MFC, aes(x=stage,
ymax=(mean+(qt(.975,degfree)*pooledSD/sqrt(N))), 
ymin=(mean-(qt(.975,degfree)*pooledSD/sqrt(N))),
color="Prioritized genes"), width=0.2)
p

# OFC (orbital prefrontal cortex)

df.OFC <- subset(df.depict,structure_acronym=="OFC")

df.ind.OFC <- ddply(df.OFC, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))                   

df.stage.OFC <- ddply(df.ind.OFC,c("stage"),summarise,mean = mean(median,na.rm=TRUE),sd = sd(median,na.rm=TRUE))

indices <- which(is.na(df.stage.OFC$sd))
bartlett.test(median~stage,subset(df.ind.OFC,!(stage %in% df.stage.OFC$stage[indices])))
leveneTest(median~stage,subset(df.ind.OFC,!(stage %in% df.stage.OFC$stage[indices])))

N <- as.data.frame(table(df.ind.OFC$stage))
colnames(N) <- c("stage","N")
df.stage.OFC <- merge(df.stage.OFC,N,by="stage",sort=FALSE)
df.for.pooledSD <- subset(df.stage.OFC,is.na(sd)==FALSE)
pooledSD <- sqrt( sum(df.for.pooledSD$sd^2 * (df.for.pooledSD$N - 1)) / (sum(df.for.pooledSD$N - 1)) )

degfree <- length(unique(df.OFC$donor_id)) - length(unique(df.for.pooledSD$stage)) - 
( length(unique(df.OFC$stage)) - length(unique(df.for.pooledSD$stage)) )

p <- ggplot()
p <- p + geom_line(data=df.stage.OFC, aes(x=stage, y=mean, group=1, color="Prioritized genes"), linetype='solid', size=1)
p
p <- p + geom_errorbar(data=df.stage.OFC, aes(x=stage,
ymax=(mean+(qt(.975,degfree)*pooledSD/sqrt(N))), 
ymin=(mean-(qt(.975,degfree)*pooledSD/sqrt(N))),
color="Prioritized genes"), width=0.2)
p

# VFC (ventrolateral prefrontal cortex)

df.VFC <- subset(df.depict,structure_acronym=="VFC")

df.ind.VFC <- ddply(df.VFC, c("stage","donor_id"), summarise, median = median(value, na.rm=TRUE), sd = sd(value, na.rm=TRUE))                   

df.stage.VFC <- ddply(df.ind.VFC,c("stage"),summarise,mean = mean(median,na.rm=TRUE),sd = sd(median,na.rm=TRUE))

indices <- which(is.na(df.stage.VFC$sd))
bartlett.test(median~stage,subset(df.ind.VFC,!(stage %in% df.stage.VFC$stage[indices])))
leveneTest(median~stage,subset(df.ind.VFC,!(stage %in% df.stage.VFC$stage[indices])))

N <- as.data.frame(table(df.ind.VFC$stage))
colnames(N) <- c("stage","N")
df.stage.VFC <- merge(df.stage.VFC,N,by="stage",sort=FALSE)
df.for.pooledSD <- subset(df.stage.VFC,is.na(sd)==FALSE)
pooledSD <- sqrt( sum(df.for.pooledSD$sd^2 * (df.for.pooledSD$N - 1)) / (sum(df.for.pooledSD$N - 1)) )

degfree <- length(unique(df.VFC$donor_id)) - length(unique(df.for.pooledSD$stage)) - 
( length(unique(df.VFC$stage)) - length(unique(df.for.pooledSD$stage)) )

p <- ggplot()
p <- p + geom_line(data=df.stage.VFC, aes(x=stage, y=mean, group=1, color="Prioritized genes"), linetype='solid', size=1)
p
p <- p + geom_errorbar(data=df.stage.VFC, aes(x=stage,
ymax=(mean+(qt(.975,degfree)*pooledSD/sqrt(N))), 
ymin=(mean-(qt(.975,degfree)*pooledSD/sqrt(N))),
color="Prioritized genes"), width=0.2)
p






