#' Constructs a binary sociomatrix from a valued sociomatrix given thresholds
#' @param sociomatrix a n times n matrix
#' @param type either A: import AND export, B: import OR export, C: only import, D: only export
#' @param threshold a percentage corresponding to the type such that the valued edge is designated as 1
#' 
#' @return a matrix with values 0 or 1 indicating if value was above the specified threshold


custom_trade_to_binary <- function(sociomatrix, type = "A", threshold = 0.01) {
  if (!(type %in% c("A", "B", "C", "D"))) stop("Wrong type specified.")
  if (!is.matrix(sociomatrix)) stop("Sociomatrix not specified as matrix.")
  if (!(threshold <= 1 & threshold >= 0)) stop("Threshold not between 0 and 1.")


  # row i -> col j, export of i, import of j
  # row_sums exports of i, col_sums imports of j
  row_sums <- apply(sociomatrix, MARGIN = 1, FUN = sum, na.rm = TRUE)
  col_sums <- apply(sociomatrix, MARGIN = 2, FUN = sum, na.rm = TRUE)


  # matrix vector comparison is column-wise
  exp_relevance <- sociomatrix >= row_sums * threshold
  imp_relevance <- t(t(sociomatrix) >= col_sums * threshold)


  if (type == "A") {
    sociomatrix_binary <- exp_relevance & imp_relevance
  } else if (type == "B") {
    sociomatrix_binary <- exp_relevance | imp_relevance
  } else if (type == "C") {
    sociomatrix_binary <- imp_relevance
  } else if (type == "D") {
    sociomatrix_binary <- exp_relevance
  }

  names(sociomatrix_binary) <- names(sociomatrix)

  return(sociomatrix_binary)
}