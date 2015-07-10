library(tools)

#===================================== USER INPUT =====================================#
#======================================================================================#

### path.rnaseq_data 
# Path *MUST* contain the files "expression_matrix.csv", "columns_metadata.csv", "rows_metadata.csv"
path.rnaseq_data <- "<INSERT PATH TO BrainSpan DATA HERE>" 
#path.rnaseq_data <- path.expand("~/Dropbox/0_Projects/p_EAv2/git/EAv2/data/141031/rnaseq")

### path.main.analysis
# Path *MUST* be writable. If it does not exists, it will be created. This directory will be used to output .RData
path.main.analysis <- "<INSERT PATH TO WHERE YOU WANT THE OUTPUT FILES TO BE WRITTEN>"
#path.main.analysis <- path.expand("~/Dropbox/0_Projects/p_EAv2/git/temporal-brain-expression/analysis")


#======================================================================================#
