library(network)
library(ergm)
library(btergm)
library(multilayer.ergm)


rm(list = ls(all.names = TRUE))
source("utils/utils.R")
path = data_path.get()


# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/cdist.RData"))
load(file.path(path, "out/milit_exp.RData"))

conflict = readRDS(file.path(path, "out/conflict_intrastate.rds"))
nmc_cinc = readRDS(file.path(path, "out/nmc_cinc.rds")) * 100
polity = readRDS(file.path(path, "out/polity.rds"))
gdp = readRDS(file.path(path, "out/gdp.rds"))
conflict_intrastate = readRDS(file.path(path, "out/conflict_intrastate.rds"))
sipri_tiv = readRDS(file.path(path, "out/sipri_tiv.rds"))
trade = readRDS(file.path(path, "out/baci_aggregated.rds"))
load(file.path(path, "out/atop_alliance.RData"))


# determine time period
start = 1995
end = 2017
maxlag = 5
iterations = 500


# helper matrices
set.seed(1234)
A = matrix(c(1,0,0,0), nrow = 2)
B = matrix(c(0,0,0,1), nrow = 2)


for (vs in c(2)){
  
  # measure computing time 
  start_time = Sys.time()
  
  # init results list
  results = list()
  simulations = list()
  
  
  for(i in (start+maxlag):end){
  
    
    # adapt indices
    i1 = i - 1949
    i2 = i - 1994
    
    included = EX[, i1] == 1
    
    
    # init empty frames
    net = list() 
    free = list() 
    help = list() 
    pathdep_l1 = list() 
    pathdep_l2 = list() 
    polity_diff_l1 = list() 
    polity_diff_l2 = list() 
    alliance_l1 = list() 
    alliance_l2 = list() 
    lgdp_in_l1 = list()
    lgdp_out_l1 = list()
    lgdp_in_l2 = list()
    lgdp_out_l2 = list()
    conflict_in_l1 = list()
    conflict_in_l2 = list()
    cinc_out_l1 = list()
    cinc_out_l2 = list()
    
    
    # include 3 lagged networks
    # layer 1: arms, layer 2: trade + binarise
    n = sum(included)
    rstr = matrix(c(rep(c(rep(1, n), rep(0, n)), n), 
                    rep(c(rep(0, n), rep(1, n)), n)), 
                  nrow = 2*n, byrow = T)
    diag(rstr) = 0 
    
    free[[4]] = network(rstr, directed = T)
    free[[3]] = network(rstr, directed = T)
    free[[2]] = network(rstr, directed = T)
    free[[1]] = network(rstr, directed = T)
    
    
    # dependent multi-layer network
    net[[4]] = to.multiplex(1 * (sipri_tiv[[i1]][included, included] > 0), 
                            1 * (trade[[i2]][included, included] > 0), 
                            output = "network")
    
    net[[3]] = to.multiplex(1 * (sipri_tiv[[i1 - 1]][included, included] > 0), 
                            1 * (trade[[i2 - 1]][included, included] > 0), 
                            output = "network")
    
    net[[2]] = to.multiplex(1 * (sipri_tiv[[i1 - 2]][included, included] > 0), 
                            1 * (trade[[i2 - 2]][included, included] > 0), 
                            output = "network")
    
    net[[1]] = to.multiplex(1 * (sipri_tiv[[i1 - 3]][included, included] > 0), 
                            1 * (trade[[i2 - 3]][included, included] > 0), 
                            output = "network")
  
    
    # define path dependency
    tmp1 = 1 * ((sipri_tiv[[i1 - 1]][included, included] + sipri_tiv[[i1 - 2]][included, included]) > 0)
    tmp2 = 1 * ((trade[[i2 - 1]][included, included] + trade[[i2 - 2]][included, included]) > 0)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    pathdep_l1[[4]] = tmp1 ; pathdep_l2[[4]] = tmp2
    
    tmp1 = 1 * ((sipri_tiv[[i1 - 2]][included, included] + sipri_tiv[[i1 - 3]][included, included]) > 0)
    tmp2 = 1 * ((trade[[i2 - 2]][included, included] + trade[[i2 - 3]][included, included]) > 0)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    pathdep_l1[[3]] = tmp1 ; pathdep_l2[[3]] = tmp2
    
    tmp1 = 1 * ((sipri_tiv[[i1 - 3]][included, included] + sipri_tiv[[i1 - 4]][included, included]) > 0)
    tmp2 = 1 * ((trade[[i2 - 3]][included, included] + trade[[i2 - 4]][included, included]) > 0)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    pathdep_l1[[2]] = tmp1 ; pathdep_l2[[2]] = tmp2
    
    tmp1 = 1 * ((sipri_tiv[[i1 - 4]][included, included] + sipri_tiv[[i1 - 5]][included, included]) > 0)
    tmp2 = 1 * ((trade[[i2 - 4]][included, included] + trade[[i2 - 5]][included, included]) > 0)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    pathdep_l1[[1]] = tmp1 ; pathdep_l2[[1]] = tmp2
    
    
    # include exogenous covariates lagged by 1 time period for each net
    
    # atop alliance
    tmp1 = kronecker(A, atop_alliance[[i1 - 1]][included, included]) 
    tmp2 = kronecker(B, atop_alliance[[i1 - 1]][included, included])
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    alliance_l1[[4]] = tmp1; alliance_l2[[4]] = tmp2
    
    tmp1 = kronecker(A, atop_alliance[[i1 - 2]][included, included]) 
    tmp2 = kronecker(B, atop_alliance[[i1 - 2]][included, included])
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    alliance_l1[[3]] = tmp1; alliance_l2[[3]] = tmp2
    
    tmp1 = kronecker(A, atop_alliance[[i1 - 3]][included, included]) 
    tmp2 = kronecker(B, atop_alliance[[i1 - 3]][included, included])
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    alliance_l1[[2]] = tmp1; alliance_l2[[2]] = tmp2
    
    tmp1 = kronecker(A, atop_alliance[[i1 - 4]][included, included]) 
    tmp2 = kronecker(B, atop_alliance[[i1 - 4]][included, included])
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    alliance_l1[[1]] = tmp1; alliance_l2[[1]] = tmp2
    
    
    # absolute polity diff
    tmp1 = kronecker(A, abs(outer(polity[included, i1-1], polity[included, i1-1],'-')))
    tmp2 = kronecker(B, abs(outer(polity[included, i1-1], polity[included, i1-1],'-')))
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    polity_diff_l1[[4]] = tmp1
    polity_diff_l2[[4]] = tmp2
    
    tmp1 = kronecker(A, abs(outer(polity[included, i1-2], polity[included, i1-2],'-')))
    tmp2 = kronecker(B, abs(outer(polity[included, i1-2], polity[included, i1-2],'-')))
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    polity_diff_l1[[3]] = tmp1
    polity_diff_l2[[3]] = tmp2
    
    tmp1 = kronecker(A, abs(outer(polity[included, i1-3], polity[included, i1-3],'-')))
    tmp2 = kronecker(B, abs(outer(polity[included, i1-3], polity[included, i1-3],'-')))
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    polity_diff_l1[[2]] = tmp1
    polity_diff_l2[[2]] = tmp2
    
    tmp1 = kronecker(A, abs(outer(polity[included, i1-4], polity[included, i1-4],'-')))
    tmp2 = kronecker(B, abs(outer(polity[included, i1-4], polity[included, i1-4],'-')))
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    polity_diff_l1[[1]] = tmp1
    polity_diff_l2[[1]] = tmp2
    
    
    # distance
    log_cdist_l1 = kronecker(A, log(cdist[included, included] + 1)) 
    log_cdist_l2 = kronecker(B, log(cdist[included, included] + 1))
    colnames(log_cdist_l1) <- rownames(log_cdist_l1) <- 1:(2*n)
    colnames(log_cdist_l2) <- rownames(log_cdist_l2) <- 1:(2*n)
    
    
    # gdp in
    tmp1 = matrix(log(gdp[included, i1-1]), n, n, byrow = TRUE)
    tmp2 = matrix(log(gdp[included, i1-1]), n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    lgdp_in_l1[[4]] = tmp1 ; lgdp_in_l2[[4]] = tmp2
    
    tmp1 = matrix(log(gdp[included, i1-2]), n, n, byrow = TRUE)
    tmp2 = matrix(log(gdp[included, i1-2]), n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    lgdp_in_l1[[3]] = tmp1 ; lgdp_in_l2[[3]] = tmp2
    
    tmp1 = matrix(log(gdp[included, i1-3]), n, n, byrow = TRUE)
    tmp2 = matrix(log(gdp[included, i1-3]), n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    lgdp_in_l1[[2]] = tmp1 ; lgdp_in_l2[[2]] = tmp2
    
    tmp1 = matrix(log(gdp[included, i1-4]), n, n, byrow = TRUE)
    tmp2 = matrix(log(gdp[included, i1-4]), n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    lgdp_in_l1[[1]] = tmp1 ; lgdp_in_l2[[1]] = tmp2
    
    
    # gdp out
    tmp1 = matrix(log(gdp[included, i1-1]), n, n, byrow = FALSE)
    tmp2 = matrix(log(gdp[included, i1-1]), n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    lgdp_out_l1[[4]] = tmp1 ; lgdp_out_l2[[4]] = tmp2
    
    tmp1 = matrix(log(gdp[included, i1-2]), n, n, byrow = FALSE)
    tmp2 = matrix(log(gdp[included, i1-2]), n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    lgdp_out_l1[[3]] = tmp1 ; lgdp_out_l2[[3]] = tmp2
    
    tmp1 = matrix(log(gdp[included, i1-3]), n, n, byrow = FALSE)
    tmp2 = matrix(log(gdp[included, i1-3]), n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    lgdp_out_l1[[2]] = tmp1 ; lgdp_out_l2[[2]] = tmp2
    
    tmp1 = matrix(log(gdp[included, i1-4]), n, n, byrow = FALSE)
    tmp2 = matrix(log(gdp[included, i1-4]), n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    lgdp_out_l1[[1]] = tmp1 ; lgdp_out_l2[[1]] = tmp2
    
    
    # conflict in
    tmp1 = matrix(conflict_intrastate[included, i1-1], n, n, byrow = TRUE)
    tmp2 = matrix(conflict_intrastate[included, i1-1], n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    conflict_in_l1[[4]] = tmp1 ; conflict_in_l2[[4]] = tmp2
    
    tmp1 = matrix(conflict_intrastate[included, i1-2], n, n, byrow = TRUE)
    tmp2 = matrix(conflict_intrastate[included, i1-2], n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    conflict_in_l1[[3]] = tmp1 ; conflict_in_l2[[3]] = tmp2
    
    tmp1 = matrix(conflict_intrastate[included, i1-3], n, n, byrow = TRUE)
    tmp2 = matrix(conflict_intrastate[included, i1-3], n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    conflict_in_l1[[2]] = tmp1 ; conflict_in_l2[[2]] = tmp2
    
    tmp1 = matrix(conflict_intrastate[included, i1-4], n, n, byrow = TRUE)
    tmp2 = matrix(conflict_intrastate[included, i1-4], n, n, byrow = TRUE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    conflict_in_l1[[1]] = tmp1 ; conflict_in_l2[[1]] = tmp2
    
    
    # cinc out
    tmp1 = matrix(nmc_cinc[included, i1-1], n, n, byrow = FALSE)
    tmp2 = matrix(nmc_cinc[included, i1-1], n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    cinc_out_l1[[4]] = tmp1 ; cinc_out_l2[[4]] = tmp2
    
    tmp1 = matrix(nmc_cinc[included, i1-2], n, n, byrow = FALSE)
    tmp2 = matrix(nmc_cinc[included, i1-2], n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    cinc_out_l1[[3]] = tmp1 ; cinc_out_l2[[3]] = tmp2
    
    tmp1 = matrix(nmc_cinc[included, i1-3], n, n, byrow = FALSE)
    tmp2 = matrix(nmc_cinc[included, i1-3], n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    cinc_out_l1[[2]] = tmp1 ; cinc_out_l2[[2]] = tmp2
    
    tmp1 = matrix(nmc_cinc[included, i1-4], n, n, byrow = FALSE)
    tmp2 = matrix(nmc_cinc[included, i1-4], n, n, byrow = FALSE)
    tmp1 = kronecker(A, tmp1) 
    tmp2 = kronecker(B, tmp2)
    colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
    colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
    cinc_out_l1[[1]] = tmp1 ; cinc_out_l2[[1]] = tmp2
    
    
    # estimate model
    if (vs == 1) { model = net ~ 
      mutual(same = "layer.mem", diff = TRUE) +
      gwidegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
      gwodegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
      edges_layer(layer = 1) +
      gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
      gwdsp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
      edgecov(alliance_l1) +
      edgecov(polity_diff_l1) +
      edgecov(log_cdist_l1) +
      edgecov(cinc_out_l1) +
      edgecov(lgdp_in_l1) + 
      edgecov(lgdp_out_l1) +
      edges_layer(layer = 2) +
      gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
      gwdsp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
      edgecov(alliance_l2) +
      edgecov(polity_diff_l2) +
      edgecov(log_cdist_l2) + 
      edgecov(cinc_out_l2) +
      edgecov(lgdp_in_l2) +
      edgecov(lgdp_out_l2) 
    }
    
    if (vs == 2) { model = net ~ 
      mutual(same = "layer.mem", diff = TRUE) +
      gwidegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
      gwodegree(decay = 1.5, fixed = TRUE, attr = "layer.mem") +
      edges_layer(layer = 1) +
      gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
      gwdsp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
      edgecov(alliance_l1) +
      edgecov(polity_diff_l1) +
      edgecov(log_cdist_l1) +
      edgecov(cinc_out_l1) +
      edgecov(lgdp_in_l1) + 
      edgecov(lgdp_out_l1) +
      edges_layer(layer = 2) +
      gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
      gwdsp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
      edgecov(alliance_l2) +
      edgecov(polity_diff_l2) +
      edgecov(log_cdist_l2) + 
      edgecov(cinc_out_l2) +
      edgecov(lgdp_in_l2) +
      edgecov(lgdp_out_l2) +
      duplexdyad(c("e", "f", "g", "h"), layers = list(1, 2))
    }
    
    fit = btergm(model, parallel = "snow", ncpus = 4, verbose = FALSE, R = iterations)
    results[[i + 1 - start - maxlag]] = fit
    
    print(paste0(vs, ": Year ", i , " estimated."))
    
    #sim = gof(fit)
    
    # save
    #simulations[[i + 1 - start - maxlag]] = sim
  
  }
  
  
  end_time = Sys.time()
  duration = end_time - start_time
  
  if (vs == 1) {
    save(model, start, end, maxlag, duration, iterations, results, 
         file = file.path(path, "models/ERGM/model_estimated_F.RData"))
    }
  
  if (vs == 2) {
    save(model, start, end, maxlag, duration, iterations, results, 
         file = file.path(path, "models/ERGM/model_estimated_G.RData"))
  }
  
}

