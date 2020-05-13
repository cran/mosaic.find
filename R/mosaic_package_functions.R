
#' @title MOSAIC (Multi-Omics Selection with Amplitude Independent Criteria)
#' @description  Function to calculate the results for RNA and protein data using the MOSAIC (Multi-Omics Selection with Amplitude Independent Criteria) method, which uses joint modeling and model selection to find both oscillatory and non-oscillatory trends in time series data.
#'
#' @param rna data frame of RNA expressions with the following specifications: first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank (NA). Labels must have same exact labels as corresponding protein data. RNA expressions with no corresponding protein data will be removed and not calculated under the MOSAIC method.
#' @param pro data frame of protein expressions with the following specifications: first column has gene labels/names, and all other columns have expression data. This expression data must be ordered by time point then by replicate, and must have evenly spaced time points. Any missing data must have cells left blank (NA). Labels must have same exact labels as corresponding RNA data. Protein expressions with no corresponding RNA data will be removed and not calculated under the MOSAIC method.
#' @param begin first time point for dataset (in hours)
#' @param end last time point for dataset (in hours)
#' @param resol resolution of time points (in hours)
#' @param num_reps number of replicates
#' @param paired if replicate data, whether the replicates are related (paired) or not (unpaired)
#' @param low lower limit when looking for rhythms, in hours. Will not be used if finding rhythms of any length within timecouse (run_all_per is TRUE).
#' @param high upper limit when looking for rhythms, in hours. Will not be used if finding rhythms of any length within timecouse (run_all_per is TRUE).
#' @param run_all_per boolean which indicates whether or not rhythms of any length within timecourse should be searched for.
#' @param rem_unexpr boolean indicating whether genes with less than rem_unexpr_amt percent expression should not be considered
#' @param rem_unexpr_amt percentage of expression for which genes should not be considered if rem_unexpr is TRUE
#' @param rem_unexpr_amt_below cutoff for expression
#' @param is_normal boolean that indicates whether data should be normalized or not
#' @param is_smooth boolean that indicates whether data should be smoothed or not
#' @param harm_cut postive number indicating the cutoff for a gene to be considered harmonic
#' @param over_cut postive number indicating the cutoff for a gene to be considered repressed/overexpressed
#' @return results, a data frame which contains:
#'   \item{Gene_Name}{gene name}
#'   \item{Best_Model_(RNA or Protein)}{The selected model type for the RNA or protein data, from a choice of oscillatory (ECHO, ECHO Joint ECHO Linear, ECHO Linear Joint) and non-oscillatory (Linear, Exponential) models.}
#'   \item{P_Value_(RNA or Protein)}{Significance of MOSAIC fit for RNA or protein data, unadjusted.}
#'   \item{BH_Adj_P_Value_(RNA or Protein)}{Significance of MOSAIC fit for RNA or protein data, adjusted using the Benjamini-Hochberg criterion. Corrects for multiple hypothesis testing.}
#'   \item{P_Value_Joint}{Significance of MOSAIC joint fit if best model is ECHO Joint or ECHO Linear Joint, unadjusted.}
#'   \item{BH_Adj_P_Value_Joint}{Significance of MOSAIC joint fit if best model is ECHO Joint or ECHO Linear Joint, adjusted using the Benjamini-Hochberg criterion. Corrects for multiple hypothesis testing.}
#'   \item{P_Value_Linear_Slope_(RNA or Protein)}{Significance of linear slope of MOSAIC fit for RNA or protein if the best model is Linear, unadjusted.}
#'   \item{BH_Adj_P_Value_Linear_Slope_(RNA or Protein)}{Significance of linear slope of MOSAIC fit for RNA or protein if the best model is Linear, adjusted using the Benjamini-Hochberg criterion. Corrects for multiple hypothesis testing.}
#'   \item{AC_Coefficient_(RNA or Protein)}{Amplitude Change Coefficient for RNA or protein. Parameter which states the amount of amplitude change over time in the system. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models.}
#'   \item{Oscillation_Type_(RNA or Protein)}{States the expression's category based on forcing coefficient (forced, damped, harmonic) for RNA or protein. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models.}
#'   \item{Initial_Amplitude_(RNA or Protein)}{Parameter describing initial amplitude of expression for RNA or protein data. Used in Exponential, ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models.}
#'   \item{Radian_Frequency_(RNA or Protein)}{Parameter describing frequency of oscillations, in radians, for RNA or protein data. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models.}
#'   \item{Period_(RNA or Protein)}{States the time for one complete oscillation, in hours, for RNA or protein. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models.}
#'   \item{Phase_Shift_(RNA or Protein)}{Parameter describing the amount the oscillator is shifted, in radians, for RNA or protein. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models.}
#'   \item{Hours_Shifted_(RNA or Protein)}{Desribes the amount the oscillator is shifted in hours, calculated from phase shift and fitted period, for RNA or protein. This is the time of the first peak of the oscillation, relative to 0 as determined by the time course entered by the user. Used in ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models.}
#'   \item{Growth_Rate_(RNA or Protein)}{Parameter describing the exponential change in amplitude for RNA or protein. Used in Exponential models.}
#'   \item{Slope_(RNA or Protein)}{Parameter describing the linear slope for RNA or protein. Used in Linear, ECHO Linear, and ECHO Linear Joint models.}
#'   \item{Equilibrium_Value_(RNA or Protein)}{Parameter describing the center, i.e. the y-intercept at time point 0, as determined by the user supplied time course, for RNA or protein. Used in Linear, Exponential, ECHO, ECHO Joint, ECHO Linear, and ECHO Linear Joint models.}
#'   \item{Processed_(RNA or Protein)_(original data column name)}{Your original data for RNA or protein, after any selected preprocessing, for specified time point (TP), using the same column names as original data.}
#'   \item{Fitted_(RNA or Protein)_TPX.R}{MOSAIC's fitted data for RNA or protein for time point (TP) X, and replicate R.}
#' @export
#' @importFrom stats p.adjust
#' @examples
#' # for more elaboration, please see the vignette
#' # "expressions_rna" is the example mosaic.find data frame for RNA
#' # "expressions_pro" is the example mosaic.find data frame for protein
#' # no preprocessing, looking for rhythms between 20 and 28 hours
#' \donttest{
#' mosaic_find(rna  = expressions_rna, pro = expressions_pro,
#' begin = 2, end = 48, resol = 2, num_reps = 3, paired = F,
#' low = 20, high = 28)
#' }
mosaic_find <- function(rna, pro, begin, end, resol, num_reps,  paired, low = 20, high = 28, run_all_per = F, rem_unexpr = F, rem_unexpr_amt = 70, rem_unexpr_amt_below=0, is_normal = F, is_smooth = F, harm_cut = .03, over_cut = .15){
  # loading and setting up parameters ----

  # creating times sequence used for the genes
  begin <- as.numeric(sapply(begin, function(x) eval(parse(text=x)))) # beginning
  end <- as.numeric(sapply(end, function(x) eval(parse(text=x)))) # end
  resol <- as.numeric(sapply(resol, function(x) eval(parse(text=x)))) # resolution of data
  timen <- seq(begin,end,resol) # the times for cicadian rhythms

  num_reps <- as.numeric(num_reps) # number of replicates

  tied <- as.logical(paired) # the type of replicate
  if (num_reps == 1){ # one replicate, default to true paired-ness
    tied <- TRUE
  }

  # subset to only run the fits in both
  both_n <- intersect(pro[,1], rna[,1])
  genes_rna <- rna[rna[,1] %in% both_n,]
  genes_pro <- pro[pro[,1] %in% both_n,]
  # order them the same
  genes_rna <- genes_rna[order(genes_rna[,1]),]
  genes_pro <- genes_pro[order(genes_pro[,1]),]

  if (nrow(genes_rna)==1){
    # creating a constant row and adding it to genes
    add_row <- data.frame(matrix(0L, 1, ncol(genes_rna)))
    add_row[1,1] <- "not considering"
    colnames(add_row) <- colnames(genes_rna)
    genes_rna <- rbind(genes_rna,add_row)

    colnames(add_row) <- colnames(genes_pro)
    genes_pro <- rbind(genes_pro,add_row)

    add_one <- TRUE # marker for appropriate displays
  } else{
    add_one <- FALSE # marker for appropriate displays
  }

  # figuring out whether a range is wanted, adjusting accordingly
  if (run_all_per){ # run all periods, a "free run"
    # low end
    if (resol >= 1){
      low <- 2*pi/resol
      low_end <- resol
    }
    else{ # if the begining is >=0, smallest period available is 1
      low <- 2*pi/1
      low_end <- 1
    }

    # high end
    high <- 2*pi/(resol*length(timen))
    high_end <- (resol*length(timen))
  } else{ # there is a high and low input
    low_input <- as.numeric(sapply(low, function(x) eval(parse(text=x))))
    low <- 2*pi/low_input
    low_end <- low_input

    high_input <- as.numeric(sapply(high, function(x) eval(parse(text=x))))
    high <- 2*pi/high_input
    high_end <- high_input
  }

  # preprocessing ----

  # removing unexpressed genes
  rem_unexpr <- rem_unexpr # indicator for removing unexpressed genes
  rem_unexpr_amt <- (rem_unexpr_amt)/100 # threshold for removing unexpressed genes, converted to a decimal
  rem_unexpr_amt_below <- abs(rem_unexpr_amt_below)
  # if yes, check for genes that are unexpressed before preprocessing
  if (rem_unexpr){
    rem_unexpr_vect_rna <- genes_unexpressed_all(genes_rna,rem_unexpr_amt,rem_unexpr_amt_below)
    rem_unexpr_vect_pro <- genes_unexpressed_all(genes_pro,rem_unexpr_amt,rem_unexpr_amt_below)
  } else{
    rem_unexpr_vect_rna <- rem_unexpr_vect_pro <- rep(FALSE,nrow(genes_rna))
  }

  rem_unexpr_combine <- (rem_unexpr_vect_rna | rem_unexpr_vect_pro)

  # normalize and store original data
  if (is_normal){
    norm_list <- normalize_all(genes_rna)
    genes_rna <- norm_list$dat

    norm_list <- normalize_all(genes_pro)
    genes_pro <- norm_list$dat
  }

  # make averages
  avg_genes_rna <- avg_all_rep(num_reps, genes_rna, timen)
  avg_genes_pro <- avg_all_rep(num_reps, genes_pro, timen)

  # smoothing, if necessary
  if (is_smooth){
    is_weighted <- T

    if (tied){ # if paired replicates
      genes_rna <- smoothing_all_tied(is_weighted, num_reps, genes_rna, avg_genes_rna, timen)
      genes_pro <- smoothing_all_tied(is_weighted, num_reps, genes_pro, avg_genes_pro, timen)
    } else{ # if unpaired replicates
      genes_rna <- smoothing_all_untied(is_weighted, num_reps, genes_rna, avg_genes_rna, timen)
      genes_pro <- smoothing_all_untied(is_weighted, num_reps, genes_pro, avg_genes_pro, timen)
    }

    # need to retake the averages because new smoothed data
    avg_genes_rna <- avg_all_rep(num_reps, genes_rna, timen)
    avg_genes_pro <- avg_all_rep(num_reps, genes_pro, timen)
  }

  # run mosaic ----

  # creating the final data frame

  all_param <- c("AC_Coefficient", "Oscillation_Type", "Period", "Hours_Shifted", "Initial_Amplitude", "Radian_Frequency", "Phase_Shift", "Growth_Rate", "Slope", "Equilibrium_Value")

  final_df <- data.frame(matrix(NA, nrow = nrow(genes_rna), ncol = 1+32+(ncol(genes_rna)-1)+(ncol(genes_pro)-1)+(length(timen))*2))
  colnames(final_df)[1] <- c("Gene_Name")
  colnames(final_df)[2:33] <-
    c("Best_Model_RNA", "Best_Model_Protein",
      "P_Value_RNA", "BH_Adj_P_Value_RNA",
      "P_Value_Protein", "BH_Adj_P_Value_Protein",
      "P_Value_Joint", "BH_Adj_P_Value_Joint",
      "P_Value_Linear_Slope_RNA", "BH_Adj_P_Value_Linear_Slope_RNA",
      "P_Value_Linear_Slope_Protein", "BH_Adj_P_Value_Linear_Slope_Protein",
      paste0(all_param,"_RNA"),
      paste0(all_param,"_Protein")
    )
  colnames(final_df)[34:ncol(final_df)] <- c(paste0("Processed_RNA_",colnames(genes_rna[-1])),
                                             paste0("Fitted_RNA_",paste0("TP",timen)),
                                             paste0("Processed_Protein_",colnames(genes_pro[-1])),
                                             paste0("Fitted_Protein_",paste0("TP",timen)))
  # add the data we already know
  # gene names
  final_df[,1] <- genes_rna[,1]
  # data
  final_df[34:(33+(length(timen)*num_reps))] <- genes_rna[,-1]
  final_df[(33+(length(timen)*num_reps)+length(timen)+1):(33+(length(timen)*num_reps)+length(timen) + (length(timen)*num_reps))] <- genes_pro[,-1]

  # no parallelism -- difficult in a package

  # run mosaic
  for (current_gene in 1:nrow(genes_rna)){
    fin <- multi_pipeline(current_gene, timen, resol, num_reps, final_df[current_gene,], genes_rna, genes_pro, avg_genes_rna, avg_genes_pro, rem_unexpr_combine, harm_cut, over_cut, low, high)
    final_df[current_gene,] <-fin$final_df
  }

  # final_df <- plyr::rbind.fill(lapply(fin, `[[`, 1))

  # change data frame values back to numeric
  final_df[,c(14,16:24,26:33)] <- sapply(final_df[,c(14,16:24,26:33)], as.numeric)

  # adjust p-values
  final_df$BH_Adj_P_Value_RNA <- p.adjust(final_df$P_Value_RNA, method = "BH")
  final_df$BH_Adj_P_Value_Protein <- p.adjust(final_df$P_Value_Protein, method = "BH")
  final_df$BH_Adj_P_Value_Joint <- p.adjust(final_df$P_Value_Joint, method = "BH")
  final_df$BH_Adj_P_Value_Linear_Slope_RNA <- p.adjust(final_df$P_Value_Linear_Slope_RNA , method = "BH")
  final_df$BH_Adj_P_Value_Linear_Slope_Protein <- p.adjust(final_df$P_Value_Linear_Slope_Protein , method = "BH")

  if (add_one){ # remove the added empty row
    final_df <- final_df[-nrow(final_df),]
  }

  return(final_df)
}
