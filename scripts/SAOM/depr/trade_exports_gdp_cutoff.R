#' trade_exports_gdp_cutoff
#' 
#' Binary Trade measure based on Exports and GDP of both partners
#' Formula I(trade)_ij = 2 * (exp_ij / (gdp_i + gdp_j) >= threshold)
#' 
#'
#' @param trd 2-dim exports matrix, country i to country j
#' @param gdp vector of gdp for countries
#' @param threshold
#'
#' @return 2-dim binary array
trade_exports_gdp_cutoff <- function(trd, gdp, threshold = 0.001) {
  # calculate trade measure
  n <- length(gdp)
  gdp_sum <- replicate(n, gdp)
  gdp_sum <- gdp_sum + t(gdp_sum)
  res <- 2 * trd/gdp_sum
  # apply threshold
  res <- (res >= threshold)*1
  return(res)
}
