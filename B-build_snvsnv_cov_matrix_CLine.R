## Build covariates table for MatrixEQTL input
gctorture(on = F)

# install/load missing required packages 
if(!require('pacman')) {install.packages('pacman'); library(pacman)}
p_load(data.table, dplyr, factoextra)

# get command line arguments
args <- commandArgs(trailingOnly = T)

# get directories of input and reference files, and SNV matrix
base_dir <- args[1] 
ref_dir <- args[2]
snv_mat <- list.files(base_dir, pattern = '-snv_matrix.tsv')

# read in SNV matrix
genotypes <- fread(file.path(base_dir, snv_mat))
genotypes_w <- genotypes[, -1]

# create prefix for output file
prfx <- args[3]

# read in sample read and phenotype reference tables
run_table <- fread(file.path(ref_dir, 'SraRunTable_full.tsv'), 
                   select = c('Run', 'submitted_subject_id', 'body_site'))
phenotypes <- fread(file.path(ref_dir, 'phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.tsv'),
                   select = c('SUBJID', 'SEX', 'AGE', 'ETHNCTY'))

# subset reference tables and merge
run_table_w <- run_table[Run %in% names(genotypes)]
phenotypes_w <- phenotypes[SUBJID %in% run_table$submitted_subject_id] 

comb_cov <- phenotypes[run_table_w, on = c('SUBJID' = 'submitted_subject_id')]

# perform PCA (10 PCs)
pca <- prcomp(data = genotypes_w, ~., na.action = 'na.omit', scale = T)
pc <- pca$x
rot <- pca$rotation
rot <- rot[,1:10]
rot <- rot %>% as.data.frame()
rot$Run <- row.names(rot)
rot <- as.data.table(rot)

comb_cov <- rot[comb_cov, on = 'Run']
comb_cov <- comb_cov[, .(Run, SEX, AGE, ETHNCTY, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)]

n <- comb_cov$Run
ID <- names(comb_cov)[-1]

# transpose all but the first column (name)
comb_cov_t <- data.table::transpose(comb_cov[, -1])
colnames(comb_cov_t) <- n
comb_cov_t <- data.table(ID, comb_cov_t)
comb_cov_t <- comb_cov_t %>% select(ID, everything())

# write the covariates file
cov_file <- file.path(base_dir, paste0(prfx, '-PC10_covariates.tsv'))
fwrite(comb_cov_t, file = cov_file, quote = F, row.names = F, sep = '\t')

# plot PCA results
pc_w <- pc[,1:10]
pca_plot <- file.path(base_dir, paste0(prfx, '-pca_plots.pdf'))
pdf(pca_plot)
scree <- fviz_eig(pca)
npcs <- length(scree[["data"]][["dim"]])
plot(pca, npcs=npcs)
fviz_pca_var(pca,
              col.var = "contrib", # Color by contributions to the PC
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE     # Avoid text overlapping
 )
dev.off()

# clean up environment if files are successfully saved
if(file.info(cov_file)$size != 0 & file.info(pca_plot)$size != 0){
  print('Output saved, clearing workspace..')
  rm(list = ls())
} else {
  print('Error: one or more of the output files is empty!')
}

gc()
