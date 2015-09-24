
###################################### SAVE PLOT ################################
saveMe <- function(plot, file_prefix, width=9, height=4) {
  filename.plot <- file.path(path.main.analysis, sprintf("%s_trajectories_rnaseq-95confint-%sx%s.pdf",file_prefix,width,height))
  filename.plot
  ggsave(file=filename.plot, plot=plot, width=width, height=height)
  print(paste("Saved plot:", filename.plot))
}


###################################### STEP 0 - final plot ################################
plotMe <- function(list_of_alphas, plot_vline=TRUE) {
  library("grid") # needed for unit() function
  
  #############################################################################################
  ######################################## PLOT 95% confidence ################################
  #############################################################################################
  
  ##########################
  ### Standard error of mean
  # SEM = SE/sqrt(n), in our case we have n="number of brain structures"=length(unique(df.summary$structure_acronym))=26
  # REFERENCE: http://en.wikipedia.org/wiki/Standard_error#Assumptions_and_usage
  
  ### 95 % confidence interval for the mean
  # FORMULA: mean +- SEM*q(alpha)
  # a) q(alpha) --> using a normal distribution | qnorm(0.975) --> ~1.96
  # b) q(alpha) --> using a t-distribution | qt(0.975,df=26-1)=2.059539
  
  
  stopifnot(length(unique(df.summary$structure_acronym))==26)
  sem.factor <<- qnorm(0.975)/sqrt(26) # GLOBAL VARIBLE | needed to when "printing" the ggplot outside the function.
   # this avoids the "object created inside function not found by ggplot" problem [e.g Error in eval(expr, envir, enclos) : object 'sem.factor' not found]
  print(paste("sem.factor:", sem.factor))
  ##########################
  
  p <- ggplot()
  ### Prioritized (for each structure)
  p <- p + geom_line(data=subset(df.summary, gene_list == "Prioritized Genes"), aes(x=stage, y=mean, group=structure_acronym, color="Prioritized genes (structures)"), alpha=list_of_alphas[["prio_struc"]])
  ### Adding mean Prioritized
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1, alpha=list_of_alphas[["prio_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1*sem.factor, ymin=mean1-sd1*sem.factor),color='#d7191c', width=0.2, alpha=list_of_alphas[["prio_mean"]])
  ### Adding mean ALL (df.all.sem)
  p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1, alpha=list_of_alphas[["all"]])
  p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1*sem.factor, ymin=mean1-sd1*sem.factor),color='black', width=0.2, alpha=list_of_alphas[["all"]])
  ### Adding Associated Genes
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, y=mean1, group=1, color="Associated genes"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, ymax=mean1+sd1*sem.factor, ymin=mean1-sd1*sem.factor), color='orange', width=0.2, alpha=list_of_alphas[["assoc_mean"]])
  ### Adding Nearest Genes
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, y=mean1, group=1, color="Nearest genes"), linetype='solid', size=1, alpha=list_of_alphas[["nearest"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Nearest Genes"), aes(x=stage, ymax=mean1+sd1*sem.factor, ymin=mean1-sd1*sem.factor), color='sky blue', width=0.2, alpha=list_of_alphas[["nearest"]])
  
  
  #### SETTING LEGEND
  p <- p + scale_color_manual(name="Gene list", values=c(
    "Prioritized genes (structures)"="gray", 
    "Prioritized genes"="#d7191c",
    "All genes"="black",
    "Nearest genes"="sky blue",
    "Associated genes"="orange",
    guide='legend'))
  
  ###### Adding vertical line - prenatal vs. postnatal
  if (plot_vline) {
    p <- p + geom_vline(xintercept=6.5, color="black", linetype="dashed")
  }
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
  
  
  ### VARIABLE
  label.y <- expression(paste("Mean brain expression ", bgroup("(",log[2]~group("[",RPKM+1,"]"),")") ))
  p <- p + labs(y=label.y)
  
  
  ##### Setting margins
  p <- p + theme(plot.margin=unit(c(5,0,0,10), "mm"))
  # ^^ margin around entire plot (unit with the sizes of the top, right, bottom, and left margins)
  
  ### MAIN FIG
  #p <- p + guides(colour = guide_legend(keywidth = 2, keyheight = 1, override.aes = list(size=c(1,1,1,1,0.1,1,1,0.1))))
  
  return(p)
}

#### Initialyze list of alpha
initialyze_alpha <- function() {
  names_plots <- c("prio_struc","prio_mean","all","assoc_mean","nearest")
  n_plots <- length(names_plots)
  list_of_alphas <- as.list(rep(0,n_plots)) # vector("list", 10) ---> did not work!
  names(list_of_alphas) <- names_plots
  return(list_of_alphas)
}