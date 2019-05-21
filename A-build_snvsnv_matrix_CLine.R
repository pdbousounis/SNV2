#!/usr/bin/env Rscript

## Build snv/snv table for MatrixEQTL
gctorture(on = F)

# get command line arguments
args <- commandArgs(trailingOnly = T)

# install required packages if missing and load 
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("data.table"); load_package("dplyr"); 


# location of the readcount files
csv_dir <- args[1]
snv_files <- list.files(csv_dir, pattern='.csv')

# create prefix for output files
prfx <- args[2]

# create an empty data table to hold the combined readcounts
dtab <- data.table()


# load in the readcounts
for (CSVfile in snv_files) {
  # create a temporary table for each individual readcount file
  temp <- fread(file.path(csv_dir, CSVfile), select = c('CHROM', 'POS', 'REF', 'ALT', 'AlignedReads', 'R'))
   
  # select the necessary columns and format the sample column
  temp[, SAMPLE := gsub(x = AlignedReads, pattern = '\\.sorted', replacement = '')][, AlignedReads := NULL]
  
  # combine all of the files together, remove temp table
  dtab <- rbindlist(list(dtab, temp))
  rm(temp)
}
gc()

# create SNP identifier
dtab[, SNP := paste0(CHROM, ':', POS, '_', REF, '>', ALT)]

# get number of samples
num_samples <- length(unique(dtab$SAMPLE))

# spread VAF values (R) by Sample and SNV
dcat <- dcast(dtab, SNP ~ SAMPLE, value.var = 'R')

# remove SNVs with > 50% (NA, VAF = 1, or VAF = 0)
dcat[ , num_NA := rowSums(is.na(dcat))]
dcat[ , perc_NA := num_NA / num_samples] 

dcat[ , num_VAF1 := rowSums(dcat > 0.9, na.rm = T)]
dcat[ , perc_VAF1 := num_VAF1 / num_samples]

dcat[ , num_VAF0 := rowSums(dcat < 0.1, na.rm = T)]
dcat[ , perc_VAF0 := num_VAF0 / num_samples]

dcat[, perc_TOT := rowSums(.SD), .SDcols = c("perc_NA", "perc_VAF1", "perc_VAF0")]
gc()

dcat2 <- dcat[perc_TOT < 0.5][ , c("num_NA", "perc_NA", 
                                   "num_VAF1", "perc_VAF1", 
                                   "num_VAF0", "perc_VAF0", 
                                   "perc_TOT") := NULL] 

d_loc <- unique(dtab[SNP %in% dcat2$SNP][, .(SNP, CHROM, POS)])
d_gene_loc <- copy(d_loc)[, POS2 := POS]

# create/specify output directory
dir_out <- paste0(csv_dir, paste0(prfx, '_SNV2_Output'))

if(!dir.exists(dir_out)) {
  dir.create(dir_out)
}

# write SNV matrix
mat_file <- file.path(dir_out, paste0(prfx, '-snv_matrix.tsv'))
fwrite(dcat2, mat_file, quote = F, row.names = F, sep = '\t')

# write location matrices
loc_file <- file.path(dir_out, paste0(prfx, '-snv_locations.tsv'))
fwrite(d_loc, loc_file, quote = F, row.names = F, sep = '\t')
gen_loc_file <- file.path(dir_out, paste0(prfx, '-gene_locations.tsv'))
fwrite(d_gene_loc, gen_loc_file, quote = F, row.names = F, sep = '\t')
gc()

if(file.info(mat_file)$size != 0 & file.info(loc_file)$size != 0 & file.info(gen_loc_file)$size != 0){
  print(paste0(paste0(prfx, '-snv_matrix.tsv, '),
        paste0(prfx, '-snv_locations.tsv, AND '),
        paste0(prfx, '-gene_locations.tsv SAVED TO '),
        dir_out))
  } 


