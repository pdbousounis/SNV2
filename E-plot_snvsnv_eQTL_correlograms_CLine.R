#!/usr/bin/env Rscript

# install required packages if missing and load 
load_package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, dep = TRUE)
    if(!require(x, character.only = TRUE)) stop(paste0("Package: ", x, " not found"))
  }
}

load_package("data.table"); load_package("corrplot"); load_package("gtools"); load_package("pheatmap"); 
load_package("RColorBrewer")

args <- commandArgs(trailingOnly = T)

  ### set base directory 
  base_dir <- args[1]
  ### select chromosome number
  chrom <- args[2]
  
  # locate input files
  snv_matrix <- list.files(base_dir, pattern = '-snv_matrix')
  cis_eqtls <- list.files(base_dir, pattern = 'cisEQTL_table.tsv')
  trans_eqtls <- list.files(base_dir, pattern = 'transEQTL_table.tsv')
  
  # create prefix for output files
  prfx <- args[3]
  
  # import/merge cis and trans tables
  cis <- fread(file.path(base_dir, cis_eqtls))
  trans <- fread(file.path(base_dir, trans_eqtls))
  cistrans <- rbindlist(list(cis, trans))
  
  #order and filter merged table by FDR
  toplot <- cistrans[order(FDR)][FDR < 0.05]
  
  # rename merged table columns, create chromosome columns for each SNP (A and B)
  names(toplot) <- c('SNPa', 'SNPb', 'beta', 't-stat', 'pvalue', 'FDR')
  toplot[, chromA := tstrsplit(SNPa, ':', keep = 1)]
  toplot[, chromB := tstrsplit(SNPb, ':', keep = 1)]
  toplot[, SNPa  := paste0('chr', SNPa)][, SNPb := paste0('chr', SNPb)]
  
  #backup
  toplot_b <- copy(toplot)
  
  # import and format SNV matrix
  NT <- fread(file.path(base_dir, snv_matrix), header = F)
  NT[, V1 := paste0('chr', V1)]
  NT[1,1] <- 'SNP'
  tNT <- data.table::transpose(NT) 
  names(tNT) <- unlist(tNT[1,])
  tNT <- tNT[-1,]
  
  toplot_c <- toplot[chromA == chrom & chromB == chrom][order(FDR)]
  
  # extract SNVs to plot (top 200 unique)
  SNPa <- head(unique(mixedsort(toplot_c$SNPa)), 500)
  SNPb <- head(unique(mixedsort(toplot_c$SNPb)), 500)
  
  # convert SNP tables to numeric matrices
  tNTa <- data.matrix(tNT[ , ..SNPa]) 
  tNTb <- data.matrix(tNT[ , ..SNPb])
  
  # build the correlation table for the selected chromosome
  system.time(cormatNT <- cor(tNTa, tNTb, use = "pairwise.complete.obs", method="spearman"))
  gc()
  
  # import gene annotations file
  sea_file <- list.files(file.path(base_dir, '..'), pattern = 'SeattleSeqAnnotation', full.names = T)
  
  if(length(sea_file) == 0){
    stop("SEATTLE SEQ ANNOTATION FILE NOT FOUND", call. = F)
  }
  
  sea <- fread(sea_file, fill = T)[, SNP := paste0('chr', chromosome, ':', position, '_', referenceBase, '>', sampleGenotype)]
  
  sea_chrom <- unique(sea, by = 'SNP')[, .(SNP, geneList)][data.table(SNP = rownames(cormatNT)), on = 'SNP'][mixedorder(SNP)]
  
  # remove duplicates and rows with geneID = 'none'
  sea_chrom[, geneID := geneList] #gsub('none', '', sea_chrom$geneList)]
  sea_chrom$geneID[duplicated(sea_chrom$geneID)] <- ''
  
  # replace correlation table row SNP IDs with gene IDs
  rownames(cormatNT) <- unlist(sea_chrom$geneID)

  # Generate data (modified the mydf slightly)
  col1 <- brewer.pal(12, "Set3")

  mymat <- matrix(cormatNT, ncol=ncol(cormatNT))
  colnames(mymat) <- sea_chrom$SNP
  rownames(mymat) <- make.names(sea_chrom$geneID, unique = T) #seq(1, nrow(mymat), 1)
  
  mydf <- data.frame(row.names = rownames(mymat), #seq(1, nrow(mymat), 1), 
                     Gene = sea_chrom$geneList)
  
 
  at <- which(sea_chrom$geneID != '')
  for(i in seq(1, length(at), 1)){at[i] <- at[i] - 1}
  at <- at[-1]
  
  plot_file <- file.path(base_dir, paste0(prfx, '_', 'chr', chrom, '_correlogram.png'))
  # add row annotations
  pheatmap(mymat, main = paste0('\n', prfx, '\n cistrans SNV-SNV Correlations'),
           cluster_cols = F, cluster_rows = F, angle_col = 45,
           cellheight = 5, cellwidth = 5,
           annotation_row = mydf, annotation_names_row = F, annotation_legend = F,
           legend = F, #border_color = NA,
           gaps_row = at, gaps_col = at,
           labels_row = sea_chrom$geneID,
           show_colnames = F,
           file = plot_file, 
           fontsize = 20, fontsize_row = 10, fontsize_col = 4.5,
           width = 35, height = 30)

gc()

if(file.info(plot_file)[["size"]] > 0){
  print(paste0(plot_file, ' SAVED TO ', base_dir))
}
  