#!/usr/bin/env Rscript

## Run MatrixEQTL: SNV vs. SNV

# get command line arguments
args <- commandArgs(trailingOnly = T)

# install required packages if missing and load 
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package('data.table'); load_package('MatrixEQTL'); load_package('pych'); load_package('gtools')


# Linear model to use (modelANOVA, modelLINEAR, or modelLINEAR_CROSS)
useModel <- modelLINEAR;

# set base directory of input files 
base_dir <- args[1] 

# create prefix for output files
prfx <- args[2]
  
snv_matrix <- list.files(base_dir, pattern = 'snv_matrix.tsv')
snv_loc <- list.files(base_dir, pattern = '-snv_locations.tsv')
gen_loc <- list.files(base_dir, pattern = '-gene_locations.tsv')
covars <- list.files(base_dir, pattern = 'covariates.tsv')

# SNPa matrix and location filenames, SNPb matrix and location filenames
SNP_file_name <- expression_file_name <- file.path(base_dir, snv_matrix) 
snps_location_file_name <- file.path(base_dir, snv_loc)
gene_location_file_name <- file.path(base_dir, gen_loc)

# Covariates file name
# Set to character() for no covariates
covariates_file_name <- file.path(base_dir, covars)

# Output file name
output_file_name_cis <- tempfile()
output_file_name_tra <- tempfile()

# Only associations significant at this level will be saved
pvOutputThreshold_cis <- 1e-3
pvOutputThreshold_tra <- 1e-3

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance <- numeric()

# Distance for local gene-SNP pairs
cisDist <- 1e6


## Load genotype data
snps <- SlicedData$new()
snps$fileDelimiter <- "\t"      # the TAB character
snps$fileOmitCharacters <- "NA" # denote missing values;
snps$fileSkipRows <- 1          # one row of column labels
snps$fileSkipColumns <- 1       # one column of row labels
snps$fileSliceSize <- 2000      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)

# load gene expression data
gene <- SlicedData$new()
gene$fileDelimiter <- "\t"      # the TAB character
gene$fileOmitCharacters <- "NA" # denote missing values;
gene$fileSkipRows <- 1          # one row of column labels
gene$fileSkipColumns <- 1       # one column of row labels
gene$fileSliceSize <- 2000      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)

# load covariates
cvrt <- SlicedData$new();
cvrt$fileDelimiter <- "\t"      # the TAB character
cvrt$fileOmitCharacters <- "NA" # denote missing values;
cvrt$fileSkipRows <- 2          # one row of column labels
cvrt$fileSkipColumns <- 1       # one column of row labels
if(length(covariates_file_name) > 0) {
  cvrt$LoadFile(covariates_file_name)
}

# load snp and gene locations
snpspos <- as.data.frame(fread(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE))
genepos <- as.data.frame(fread(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE))

# run the analysis
me <- Matrix_eQTL_main(
  snps = snps, 
  gene = gene, 
  cvrt = cvrt,
  output_file_name = output_file_name_tra,
  pvOutputThreshold  = pvOutputThreshold_tra,
  useModel = useModel, 
  errorCovariance = errorCovariance, 
  verbose = TRUE, 
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos, 
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = "qqplot", # choose qq-plot 
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE)

cat("Analysis done in: ", me$time.in.sec, ' seconds', '\n')

# save qq-plot 
plot_file <- file.path(base_dir, paste0(prfx, '-matrixEQTL_results.pdf'))
pdf(plot_file)
plot(me)
dev.off()

# remove snpA = snpB pairs, format tables for plotting 
cis_eqtls <-  as.data.table(me$cis$eqtls)[snps != gene][, gene := as.character(unlist(gene))][, snps := as.character(unlist(snps))]
trans_eqtls <- as.data.table(me$trans$eqtls)[, gene := as.character(unlist(gene))][, snps := as.character(unlist(snps))]

# write output files to base directory
cis_file <- file.path(base_dir, paste0(prfx, '-cisEQTL_table.tsv'))
fwrite(cis_eqtls, cis_file, quote = F, sep = '\t')

trans_file <- file.path(base_dir, paste0(prfx, '-transEQTL_table.tsv'))
fwrite(trans_eqtls, trans_file, quote = F, sep = '\t')

gc()

if(file.info(cis_file)$size != 0 & file.info(trans_file)$size != 0){
  print(paste0(paste0(prfx, '-cisEQTL_table.tsv, '), 
               paste0(prfx, '-transEQTL_table.tsv, AND '),
               paste0(prfx, '-matrixEQTL_results.pdf SAVED TO '),
               base_dir))
} 


