#' Calculate import relevance of trade relationship in a given year
#' 
#' 
#'
#' @param trd matrix of trade flow (export from row to column actor)
#' @param gdp vector of GDPs
#' @param threshold 
#'
#' @return
#' @export
#'
#' @examples
trade_import_gdp_relevance <- function(trd, gdp, threshold = 0.004) {
  imports <- t(trd)
  n <- dim(imports)[2]
  gdp_rep <- replicate(n, gdp)
  res <- imports / gdp_rep
  
  if (is.null(threshold)) {
    return(res)
  } else {
    res <- (res >= threshold)
    res <- 1*res
    return(res)
  }
}