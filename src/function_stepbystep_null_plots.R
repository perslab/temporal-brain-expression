
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
  
  sem.factor <- qnorm(0.975)/2 # 95% confidence interval (plus/minus)
  print(sem.factor)
  # *OBS* CORRECTION! 2015-09-12!!!
  
  p <- ggplot()
  ### Adding mean Prioritized
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, y=mean1, group=1, color="Prioritized genes"), linetype='solid', size=1, alpha=list_of_alphas[["prio_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Prioritized Genes"), aes(x=stage, ymax=mean1+sd1*sem.factor, ymin=mean1-sd1*sem.factor),color='#d7191c', width=0.2, alpha=list_of_alphas[["prio_mean"]])
  ### Adding mean ALL (df.all.sem)
  p <- p + geom_line(data=df.all.sem, aes(x=stage, y=mean1, group=1, color="All genes"), linetype='solid', size=1, alpha=list_of_alphas[["all"]])
  p <- p + geom_errorbar(data=df.all.sem, aes(x=stage, ymax=mean1+sd1*sem.factor, ymin=mean1-sd1*sem.factor),color='black', width=0.2, alpha=list_of_alphas[["all"]])
  ### Adding Associated Genes
  p <- p + geom_line(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, y=mean1, group=1, color="Associated genes"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_mean"]])
  p <- p + geom_errorbar(data=subset(df.summary.sem, gene_list == "Associated Genes"), aes(x=stage, ymax=mean1+sd1*sem.factor, ymin=mean1-sd1*sem.factor), color='orange', width=0.2, alpha=list_of_alphas[["assoc_mean"]])

  ### Adding NULL ASSOCIATED (median) - ribbon!
  p <- p + geom_line(data=df.null.median.summary.sem.assoc, aes(x=stage, y=mean1, group=1, color="Associated genes (Null)"), linetype='solid', size=1, alpha=list_of_alphas[["assoc_null_line"]])
  #p <- p + geom_ribbon(data=df.null.median.summary.sem.assoc, aes(x=stage, group=1, ymin=mean1-sd1*sem.factor, ymax=mean1+sd1*sem.factor), alpha=list_of_alphas[["assoc_null_rib"]], fill='darkolivegreen4')
  p <- p + geom_ribbon(data=df.null.median.summary.sem.assoc, aes(x=stage, group=1, ymin=mean1-sd1, ymax=mean1+sd1), alpha=list_of_alphas[["assoc_null_rib"]], fill='darkolivegreen4')
  p
  ### Adding NULL PRIORITIZED (median) - ribbon!
  p <- p + geom_line(data=df.null.median.summary.sem.prio, aes(x=stage, y=mean1, group=1, color="Prioritized genes (Null)"), linetype='solid', size=1, alpha=list_of_alphas[["prio_null_line"]])
  #p <- p + geom_ribbon(data=df.null.median.summary.sem.prio, aes(x=stage, group=1, ymin=mean1-sd1*sem.factor, ymax=mean1+sd1*sem.factor), alpha=list_of_alphas[["prio_null_rib"]], fill='lightskyblue3')
  p <- p + geom_ribbon(data=df.null.median.summary.sem.prio, aes(x=stage, group=1, ymin=mean1-sd1, ymax=mean1+sd1), alpha=list_of_alphas[["prio_null_rib"]], fill='lightskyblue3')
  p
  
  
  #### SETTING LEGEND
  p <- p + scale_color_manual(name="Gene list", values=c(
    "Prioritized genes (structures)"="gray", 
    "Prioritized genes"="#d7191c",
    "All genes"="black",
    "Nearest genes"="sky blue",
    "Associated genes"="orange",
    "Associated genes (Null)"="darkolivegreen4",
    "Prioritized genes (Null)"="lightskyblue3",
    guide='legend'))
  
  #"Associated genes"="darkgoldenrod1",
  #"Nearest genes"="darkgoldenrod3",
  
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
  names_plots <- c("prio_mean","all","assoc_mean")
  names_plots <- c(names_plots, "prio_null_line","prio_null_rib","assoc_null_line","assoc_null_rib")
  n_plots <- length(names_plots)
  list_of_alphas <- as.list(rep(0,n_plots)) # vector("list", 10) ---> did not work!
  names(list_of_alphas) <- names_plots
  return(list_of_alphas)
}
