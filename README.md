# temporal-brain-expression
Scripts to analyze spatio temporal brain expression based on the BrainSpan Developmental Transcriptome data.

### BrainSpan Data
Click [here](http://www.brainspan.org/api/v2/well_known_file_download/267666525) to download the RNA-seq data (genes_matrix_csv.zip, 62.2 MB).

For more information, please refer to the [BrainSpan Developmental Transcriptome webpage](http://www.brainspan.org/static/download.html).

### Documentation
1. Download the [BrainSpan RNA-seq data](http://www.brainspan.org/api/v2/well_known_file_download/267666525) (same file as referred to above).
2. Enter the correct paths in `CONFIG.R`. This is needed to load the BrainSpan data and control where the output files will be written.
3. Change to the `src/` directory.
4. Run `R CMD BATCH read_rnaseq_data.R` to load and process the BrainSpan data. The script generates two `.RData` files with the BrainSpan data (Run time is approx. 15 min.).
5. Perform analysis
  1. Run `R CMD BATCH graphics_genes_temporal_trajectories.R` to produce gene trajectory plots.
  2. Run `R CMD BATCH statistics_prenatal-vs-postnatal-test.R` to compare expression levels in prenatal vs postnatal developmental stages.

Individual scripts contains comments explaining the analysis steps.

### System Requirements
* R (version >= 3.1)
* R packages: plyr, dplyr, ggplot2, reshape2 and tools
* 12 GB RAM

The scripts have only been tested on OSX and Linux operation systems.
