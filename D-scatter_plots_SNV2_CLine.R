  library(pacman) 
  p_load(data.table, dplyr, ggplot2, ggpubr, gtools, stringr, tidyr)
  detach("package:psych", unload=TRUE)
  
  # get command line arguments
  args <- commandArgs(trailingOnly = T)
  
  # directory of output files from scripts A-C
  base_dir <- args[1]
  
  # import the eQTL table
  pattern <- paste0('*', args[2], 'EQTL_table')
  in_file <- fread(list.files(base_dir, pattern =  pattern, full.names = T), header = T)
  in_file <- in_file[mixedorder(in_file$FDR), ][FDR < 0.05]
  names(in_file)[1:2] <- c('SNPa', 'SNPb')
  
  # designate output file prefix
  prfx <- args[3]
  
  to_plot <- as.data.frame(in_file) %>% dplyr::select(SNPa, SNPb) 
  to_plot <- to_plot[seq(1, nrow(to_plot), 2), ]
  
  vafA_data <- fread(list.files(base_dir, pattern = 'snv_matrix', full.names = T), header = T) %>% filter(SNP %in% to_plot$SNPa)
  vafA_data <- vafA_data %>% gather(sample, vaf, 2:ncol(vafA_data)) %>% drop_na()
  names(vafA_data) <- c('SNPa', 'sample', 'vafA')
  
  vafB_data <- fread(list.files(base_dir, pattern = 'snv_matrix', full.names = T), header = T) %>% filter(SNP %in% to_plot$SNPb)
  vafB_data <- vafB_data %>% gather(sample, vaf, 2:ncol(vafB_data)) %>% drop_na()
  names(vafB_data) <- c('SNPb', 'sample', 'vafB')
  
  # join the tables
  df <- left_join(to_plot, vafA_data)
  df <- left_join(df, vafB_data)
  df <- df %>% drop_na() %>% unite(id, SNPa, SNPb)
  
  uq200FDR <- head(sort(unique(in_file$FDR)), 200)
  FDR2df <- in_file[which(FDR %in% uq200FDR), ] %>% unite(id, SNPa, SNPb)
  
  dff <- df[which(df$id %in% FDR2df$id), ]
  
  if(grepl('trans', pattern)){
    col1 <- 'coral1'
    col2 <- 'coral4'
  } else {
    col1 <- 'forestgreen'
    col2 <- 'darkslategrey'
  }
  
  p <- ggscatter(dff, x = "vafA", y = "vafB",
                fill = col1, color = col1, shape = 21, size = 2, # Points color, shape and size
                add = "reg.line",  # Add regression line
                add.params = list(color = col2, fill = alpha(col2), 0.5), # Customize reg. line
                cor.coef = T, # Add correlation coefficient. see ?stat_cor
                cor.coef.size = 10, cor.coef.coord = c(0,0.85),
                cor.coeff.args = list(method = "spearman", label.sep = "\n"),
                ylim = c(0,1), xlim = c(0,1), font.ytickslab = c(30), font.xtickslab = c(30))
  
  fname <- ifelse(grepl('*cis*', pattern), 
                 paste0(prfx, '-cisEQTL_scatter.png'),
                 paste0(prfx, '-transEQTL_scatter.png'))
  
  
  plot_file <- file.path(base_dir, fname)
  #height <- ceiling(nrow(in_file)/5)*10
  #ifelse(nrow(in_file < 5), width <- nrow(in_file)*10, width <- 50)
  png(plot_file, height = 400, width = 50, units = 'in', res = 120)
  facet(p, facet.by = "id", ncol = 5, scales = "free_y",
        panel.labs.font = list(face = 'bold', size = 30))
  dev.off()
  
  # clean up environment
  if(file.info(plot_file)$size != 0){
    print('Output saved, clearing workspace..')
    rm(list = ls())
    gc()
  }
  
