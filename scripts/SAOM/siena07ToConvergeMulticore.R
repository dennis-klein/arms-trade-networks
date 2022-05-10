require(RSiena)
require(parallel)

siena07ToConvergenceMulticore <- function(alg, dat, eff, save_dir, ans_id,
                                 ans0=NULL, threshold=0.25, ...){
  n.clus <- detectCores() - 1
  
  numr <- 0
  ans <- siena07(alg, data=dat, effects=eff, prevAns=ans0,
                 useCluster = TRUE, nbrNodes = n.clus, initC = TRUE,
                 ...) # the first run
  repeat {
    saveRDS(ans, file=path(save_dir, paste0("tmp_fit_", ans_id, "_", numr), ext = "rds")) # to be safe
    numr <- numr+1 # count number of repeated runs
    tm <- ans$tconv.max # convergence indicator
    cat(numr, tm,"\n") # report how far we are
    if (is.na(tm)) {break} # calculation of convergence ratio not possible
    if (tm < threshold) {break} # success
    if (tm > 2) {break} # divergence without much hope of good return
    if (numr > 4) {break} # now it has lasted too long
    ans <- siena07(alg, data=dat, effects=eff, prevAns=ans,
                   useCluster = TRUE, nbrNodes = n.clus, initC = TRUE, ...)
  }
  saveRDS(ans, file=path(save_dir, paste0("final_fit_", ans_id), ext = "rds"))
  return(ans)
}