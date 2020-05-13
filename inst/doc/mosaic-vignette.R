## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  # collapse = FALSE,
  # comment = "#>"
  results = "hold"
)

## ------------------------------------------------------------------------
library(mosaic.find)

head(expressions_rna[,1:5])
head(expressions_pro[,1:5])

## ----fig.align="center",fig.height = 5, fig.width = 7,warning=FALSE------

library(ggplot2)

# making a function to plot expression, since we're going to use it again
# samp: sample name that we want to visualize
# expressions_rna, expressions_pro: our original RNA and protein expression data frames
# results: data frame of mosaic results
plot_exp <- function(samp, expressions_rna, expressions_pro, results = data.frame()){
  # color information
  rna_blue_high <- "#268beb"
  pro_red_high <- "#d61e1e"
  
  rna_blue_low <- "#9ebfde"
  pro_red_low <- "#d49d9d"
  
  # making a function to create automatic palettes for rna and protein vizualizations
  rna_col_func <- colorRampPalette(c(rna_blue_low, rna_blue_high))
  pro_col_func <- colorRampPalette(c(pro_red_low, pro_red_high))
  
  tp <- seq(2,48,by=2) # our time points
  num_reps <- 3 # number of replicates
  
  # set to plot either the mosaic results or the original data
  if (nrow(results) == 0){ 
    # we're plotting the original data
    dat_rna <- expressions_rna[,-1]
    dat_pro <- expressions_pro[,-1]
    
    # logicals for the position of our desired sample in each dataset
    gene_log_rna <- expressions_rna$Sample_Name == samp
    gene_log_pro <- expressions_pro$Sample_Name == samp
  } else {
    # last number that is a parameter in the results dataset
    end_num <- 33

    # we're plotting the processed original data and the fitted data
    dat_rna <- results[,(end_num+1):((end_num+1)+((length(tp)*num_reps)-1))]
    fit_rna <- results[,(end_num+(length(tp)*num_reps)+1):((end_num+(length(tp)*num_reps)+1)+(length(tp)-1))]
    dat_pro <- results[,(end_num+(length(tp)*num_reps)+1+length(tp)):(end_num+(length(tp)*num_reps*2)+1+length(tp)-1)]
    fit_pro <- results[,(end_num+(length(tp)*num_reps*2)+1+length(tp)):ncol(results)]
    
    # logicals for the position of our desired sample
    gene_log_rna <- gene_log_pro <- results$Gene_Name == samp
  }
  
  # our visualization data frame
  gg.df <- data.frame(
    # maximum over time points
    "rna_max" = sapply(seq(1,ncol(dat_rna), by = num_reps), function(x) max(dat_rna[gene_log_rna,x:(x+num_reps-1)], na.rm = T)),
    "pro_max" = sapply(seq(1,ncol(dat_pro), by = num_reps), function(x) max(dat_pro[gene_log_pro,x:(x+num_reps-1)], na.rm = T)),
    # minimum over time points
    "rna_min" = sapply(seq(1,ncol(dat_rna), by = num_reps), function(x) min(dat_rna[gene_log_rna,x:(x+num_reps-1)], na.rm = T)),
    "pro_min" = sapply(seq(1,ncol(dat_pro), by = num_reps), function(x) min(dat_pro[gene_log_pro,x:(x+num_reps-1)], na.rm = T)),
    "timen" = tp # time points
  )
  
  # colors for lines and shading
  col_vect <- c(
    "Original RNA" = rna_blue_low,
    "Original Protein" = pro_red_low
  )
  
  # only plot individual replicates if we're only plotting the original data
  if (nrow(results) == 0){
    # add replicate colors
    # add two to provide more variability
    rna_rep_col <- rna_col_func(num_reps+2)
    pro_rep_col <- pro_col_func(num_reps+2)
    
    # add replicates to df and color scale
    for (i in 1:num_reps){
      gg.df[,paste0("rna_rep",i)] <- as.numeric(dat_rna[gene_log_rna,seq(i, ncol(dat_rna), by=num_reps)])
      gg.df[,paste0("pro_rep",i)] <- as.numeric(dat_pro[gene_log_pro,seq(i, ncol(dat_pro), by=num_reps)])
    
      col_vect[paste0("RNA, Replicate ",i)] <- rna_rep_col[i+1]
      col_vect[paste0("Protein, Replicate ",i)] <- pro_rep_col[i+1]
    }
  } else {
    # otherwise, we plot the model fits
    gg.df$rna_fit <- as.numeric(fit_rna[gene_log_rna,])
    gg.df$pro_fit <- as.numeric(fit_pro[gene_log_pro,])
    
    # add the colors
    col_vect["Fit RNA"] <- rna_blue_high
    col_vect["Fit Protein"] <- pro_red_high
  }
  
  # name breaks for the colors in the legend
  col_name_breaks <- c(names(col_vect)[grepl("RNA", names(col_vect), fixed = T)],
                       names(col_vect)[grepl("Protein", names(col_vect), fixed = T)])
  
  # visualize, with shading for each omics type
  p <- ggplot(gg.df, aes(x = timen))+
    # setting up scales and labels
    scale_fill_manual("", values = col_vect, breaks = col_name_breaks)+
    scale_color_manual("", values = col_vect, breaks = col_name_breaks)+
    ylab("Expression")+
    xlab("Time (Hours)")+
    theme_bw()+
    theme(
      plot.title = element_text(hjust = .5),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )+
    scale_x_continuous(expand = c(0, 0))+
    # add shading for original rna and protein data
    geom_ribbon(aes(ymax = rna_max, ymin = rna_min, fill = "Original RNA"), alpha = .5)+
    geom_ribbon(aes(ymax = pro_max, ymin = pro_min, fill = "Original Protein"), alpha = .5)+
    ggtitle(samp)+ # title
    guides( # colors
      color = guide_legend(nrow = 2, byrow = T),
      fill = guide_legend(nrow = 2)
    )+
    NULL
  
  
  # only plot individual replicates if we're only plotting the original data
  if (nrow(results) == 0){
    # add each of the replicates
    for (i in 1:num_reps){
      p <- p +
        geom_line(aes_string(y = paste0("rna_rep",i), color = shQuote(paste0("RNA, Replicate ",i))))+
        geom_line(aes_string(y = paste0("pro_rep",i), color = shQuote(paste0("Protein, Replicate ",i))))+
        NULL
    }
  } else { # otherwise, plot the fits
    p <- p +
      geom_line(aes(y = rna_fit, color = "Fit RNA"), size = 1)+
      geom_line(aes(y = pro_fit, color = "Fit Protein"), size = 1)+
      NULL
  }
  
  return(p)
}

suppressWarnings(plot_exp("Sample_2", expressions_rna, expressions_pro)) # to ignore warnings for missing values


## ------------------------------------------------------------------------
# run mosaic_find()
results <- mosaic_find(
  rna  = expressions_rna, # rna data
  pro = expressions_pro, # protein data
  begin = 2, # first time point (hours)
  end = 48, # last time point (hours)
  resol = 2, # resolution (hours)
  num_reps = 3, # number of replicates
  paired = F, # unpaired, or unrelated, replicates
  low = 20, # lower limit when looking for rhythms (hours)
  high = 28 # lower limit when looking for rhythms (hours)
)

# look at the results!
head(results[,1:16], 5)

## ----fig.align="center",fig.height = 5, fig.width = 7--------------------

# making a function, since we'll be using this again later
plot_dist <- function(results){
  # first, we only want to consider significant expressions
  sig_rna <- results$BH_Adj_P_Value_RNA < .05
  sig_pro <- results$BH_Adj_P_Value_Protein < .05
  
  # then get the significant models for each omics type
  mod_rna <- results$Best_Model_RNA[sig_rna]
  mod_pro <- results$Best_Model_Protein[sig_pro]
  
  # form into the data frame for plotting
  gg.df <- data.frame(
    "Model" = c(mod_rna, mod_pro),
    "Omics_Type" = c(rep("RNA", length(mod_rna)), rep("Protein",length(mod_pro)))
  )
  
  # colors
  col_proper <- c("Linear" = "#b32515", # red
           "Exponential" = "#e39f22", # orange
           "ECHO" = "#3461eb", # blue
           "ECHO Joint" = "#5c7de0", # light blue
           "ECHO Linear" = "#7b2bcc", # purple
           "ECHO Linear Joint" = "#945fc9" # light purple
  )
  
  # make a bar chart of distributions for each model type
  ggplot(gg.df, aes(Omics_Type, fill = Model))+
    geom_bar()+
    theme_bw()+
    ylab("Total Significant Expressions")+
    xlab("Omics Type")+
    theme_bw()+
    theme(
      plot.title = element_text(hjust = .5),
      legend.position = "bottom",
      legend.direction = "horizontal"
    )+
    ggtitle("MOSAIC Model Distributions")+ # title
    guides( # colors
      fill = guide_legend(nrow = 2)
    )+
    scale_y_continuous(expand = expand_scale(mult = c(0, .1)))+
    scale_fill_manual("Model", values = col_proper)+
    NULL
}

# plot our model distributions!
plot_dist(results)

## ----fig.align="center",fig.height = 5, fig.width = 7--------------------

# plot our fitted data!
suppressWarnings(plot_exp("Sample_2", expressions_rna, expressions_pro, results)) # to ignore warnings for missing values

# print our parameters for our sample
results[results$Gene_Name == "Sample_2",]


## ------------------------------------------------------------------------
# run mosaic_find() with preprocessing
results <- mosaic_find(
  rna  = expressions_rna, # rna data
  pro = expressions_pro, # protein data
  begin = 2, # first time point (hours)
  end = 48, # last time point (hours)
  resol = 2, # resolution (hours)
  num_reps = 3, # number of replicates
  paired = F, # unpaired, or unrelated, replicates
  run_all_per = T, # now we're looking for any period within our time course
  rem_unexpr = T, # remove any unexpressed genes (we don't have any, but some time courses might)
  rem_unexpr_amt = 70, # if we have less than 70% expression, we want to remove the gene
  rem_unexpr_amt_below = 0, # unexpressed genes are those with expression equal to 0
  is_normal = T, # normalize data
  is_smooth = T # smooth the data
)

# look at the results!
head(results[,1:16], 5)


## ----fig.align="center",fig.height = 5, fig.width = 7--------------------

# plot our model distributions!
plot_dist(results)

# plot our fitted data!
suppressWarnings(plot_exp("Sample_2", expressions_rna, expressions_pro, results)) # to ignore warnings for missing values

# print our parameters for our sample
results[results$Gene_Name == "Sample_2",]


