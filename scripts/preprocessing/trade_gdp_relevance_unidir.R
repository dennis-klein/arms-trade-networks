#' trade_gdp_relevance_unidir
#' 
#' Unidirectional, binary measure of relevance of the trade relationship
#' for a country
#' 
#' Based on the measure of trade openness (trade flow / GDP)
#' 
#' Formula:
#' I(trade)_ij = ((exp_ij + imp_ij) / gdp_i) >= threshold
#' 
#'
#' @param trd 2-dim trade (exports) matrix, country i to country j
#' @param gdp vector of gdp for countries
#' @param threshold
#'
#' @return 2-dim binary array
trade_gdp_relevance_unidir <- function(trd, gdp, threshold = 0.01) {
  n <- length(gdp)
  trd_sum <- trd + t(trd)
  gdp_rep <- replicate(n, gdp)
  res <- trd_sum / gdp_rep
  res <- (res >= threshold)*1
  return(res)
}
