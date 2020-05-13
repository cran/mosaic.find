#' Synthetic RNA expression data for 600 genes.
#'
#' A dataset containing the names and RNA expression values for 600 synthetically
#' generated samples. This example data has time points from 2 to 48 hours with
#' 2 hour resolution and 3 replicates. Random missing data is also included. Expression
#' labels are the same as the corresponding protein dataset, expressions_pro. Synthetic
#' data was created by randomly selecting parameters for the ECHO, ECHO Joint, ECHO
#' Linear, ECHO Linear Joint, linear, and exponential models, then adding random
#' noise to each expression, as described in H. De los Santos, et al. (2020).
#' There is a 2:1 ratio of non-oscillatory to oscillatory models. See linked paper in
#' vignette for more information.
#'
#' Note the data format: its first column first column has gene labels/names, and
#' all other columns have expression data. This expression data is ordered by
#' time point then by replicate, and has evenly spaced time points. Any missing
#' data has cells left blank (NA). Labels have same exact labels as corresponding protein
#' data.
#'
#' @format A data frame with 600 rows and 73 variables (column 1: sample labels, columns to 2 to 73: numerical values for RNA expression in the format TPX.Y (time point X, replicate Y)).
#'
"expressions_rna"

#' Synthetic protein expression data for 600 genes.
#'
#' A dataset containing the names and protein expression values for 600 synthetically
#' generated samples. This example data has time points from 2 to 48 hours with
#' 2 hour resolution and 3 replicates. Random missing data is also included. Expression
#' labels are the same as the corresponding protein dataset, expressions_rna. Synthetic
#' data was created by randomly selecting parameters for the ECHO, ECHO Joint, ECHO
#' Linear, ECHO Linear Joint, linear, and exponential models, then adding random
#' noise to each expression, as described in H. De los Santos, et al. (2020).
#' There is a 2:1 ratio of non-oscillatory to oscillatory models. See linked paper in
#' vignette for more information.
#'
#' Note the data format: its first column first column has gene labels/names, and
#' all other columns have expression data. This expression data is ordered by
#' time point then by replicate, and has evenly spaced time points. Any missing
#' data has cells left blank (NA). Labels have same exact labels as corresponding RNA
#' data.
#'
#' @format A data frame with 600 rows and 73 variables (column 1: sample labels, columns to 2 to 73: numerical values for protein expression in the format TPX.Y (time point X, replicate Y)).
#'
"expressions_pro"
