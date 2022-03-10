#' trade_exports_gdp_cutoff
#' trade_flow_gdp_cutoff
#' 
#' Binary Trade measure based on exports, imports and GDP.
#' Formula: Indicator_i = (flow_ij/gdp_i >= threshold).
#' Flow can be defined as total flow (import and export) or
#' only export (equivalent to only import because of symmetry).
#' To get the final indicator for the relation between i and j
#' OR or AND operators and can be applied.
#'
#' @param trd 2-dim trade (exports) matrix, country i to country j
#' @param gdp vector of gdp for the countries
#' @param threshold threshold for the flow-to-GDP ratio
#' @param flow "both" (default) or "exports", see above
#' @param operation "OR" (default) or "AND", see above
#'
#' @return 2-dim binary array
trade_flow_gdp_cutoff <- function(trd, gdp, threshold = 0.01, flow="both", operation = "OR") {
  if (flow == "both") {
    flow_val <- trd + t(trd)
  } else if (flow == "imports") {
    flow_val <- t(trd)
  } else if (flow == "exports") {
    flow_val <- trd
  } else {
    stop("Wrong flow type.")
  }
  
  n <- length(gdp)
  gdp_i <- replicate(n, gdp)
  gdp_j <- t(gdp_i)
  
  rel_i <- flow_val/gdp_i
  rel_j <- flow_val/gdp_j
  
  rel_i <- (rel_i >= threshold)
  rel_j <- (rel_j >= threshold)
  
  if (operation == "OR") {
    rel <- rel_i | rel_j
  } else if (operation == "AND") {
    rel <- rel_i & rel_j
  } else {
    stop("Wrong operation.")
  }
  rel <- rel*1
}
