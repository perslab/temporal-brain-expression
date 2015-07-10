# temporal-brain-expression
Scripts to analyze spatio temporal brain expression of BrainSpan Developmental Transcriptome data.

### BrainSpan Data
Click [here](http://www.brainspan.org/api/v2/well_known_file_download/267666525) to download the RNAseq data (genes_matrix_csv.zip, 62.2 MB).

For more information, see the BrainSpan Developmental Transcriptome webpage: http://www.brainspan.org/static/download.html

### Documentation
1. Download the BrainSpan data (click [here](http://www.brainspan.org/api/v2/well_known_file_download/267666525) to download).
2. Enter the correct paths in `CONFIG.R`. This is needed to load the BrainSpan data and control where the output files of the script goes.
3. Run `read_rnaseq_data.R` to load and process the BrainSpan data. The script generates two `.RData` files with the BrainSpan data.
4. Perform analysis
  1. Run `graphics_genes_temporal_trajectories.R` to produce gene trajectory plots.
  2. Run `statistics_prenatal-vs-postnatal-test.R` to compare expression levels in prenatal vs postnatal developmental stages.

Individual scripts contains comments explaining the analysis steps.

### System Requirements
* R (version >= 3.1)
* R packages: plyr, dplyr, ggplot2, reshape2 and tools
* 12 GB RAM

The scripts have only been tested on OSX and Linux operation systems.
