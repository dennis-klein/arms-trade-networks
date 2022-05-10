library(network)
library(ergm)
library(btergm)
library(multilayer.ergm)
library(PRROC)


# init setup
rm(list = ls(all.names = TRUE))
set.seed(1234)
source("utils/utils.R")
source("utils/construct_header.R")
source("utils/custom_trade_to_binary.R")
path = data_path.get()


# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/cdist.RData"))
load(file.path(path, "out/atop_alliance.RData"))
load(file.path(path, "out/milit_exp.RData"))
nmc_cinc = readRDS(file.path(path, "out/nmc_cinc.rds"))
polity = readRDS(file.path(path, "out/polity.rds"))
gdp = readRDS(file.path(path, "out/gdp.rds"))
gdppc = readRDS(file.path(path, "out/gdppc.rds"))
sipri_tiv = readRDS(file.path(path, "out/sipri_tiv.rds"))
trade = readRDS(file.path(path, "out/baci_aggregated.rds"))


# number of simulations for prediction in unconditional
n_pred <- 100
version <- "A"
type <- "D"
start <- 1995
end <- 2018
maxlag <- 4


# selection of countries included
included <- rowSums(EX[, (start:end) - 1949]) == length(start:end)
n <- sum(included)


# some data transformations
# sipri military expenditure is in constant USD$m
# cinc is index which sums up to 0 -> better interpretability using pct points
milit_exp <- milit_exp * 1000000
nmc_cinc100 <- nmc_cinc * 100


# setup constraints
free <- matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2 * n, byrow = T)
diag(free) <- 0

free_arms <- free
free_arms[(n + 1):(2 * n), (n + 1):(2 * n)] <- 0


# matrices for kronecker products
A <- matrix(c(1, 0, 0, 0), nrow = 2)
B <- matrix(c(0, 0, 0, 1), nrow = 2)


# save predictions
scores <- data.frame(
  "model" = character(),
  "year" = numeric(),
  "metric" = character(),
  "type" = character(),
  "value" = numeric()
)


# out-of-samples predictions
for (year in (start+maxlag):(end-1)){
  
  load(file = paste0(path, "/models/ERGM/", version, "_", year , "_complete.RData"))
  
  
  # set correct indices corresponding data source
  i1 = year - 1949 # Cornelius Replication Data 1950:2018
  i2 = year - 1994 # CEPII BACI 1995:2018
  

  # load corresponding t+1 covariates
  future = to.multiplex(1* (sipri_tiv[[i1+1]][included, included] > 0), 
                        custom_trade_to_binary(trade[[i2+1]][included, included], type = type, threshold = 0.01), 
                        output = "network")
  current = to.multiplex(1* (sipri_tiv[[i1]][included, included] > 0), 
                        custom_trade_to_binary(trade[[i2+1]][included, included], type = type, threshold = 0.01), 
                        output = "network")
  
  # distance
  log_cdist_layer1 = kronecker(A, log(cdist[included, included] + 1)) 
  log_cdist_layer2 = kronecker(B, log(cdist[included, included] + 1))
  colnames(log_cdist_layer1) <- rownames(log_cdist_layer1) <- 1:(2*n)
  colnames(log_cdist_layer2) <- rownames(log_cdist_layer2) <- 1:(2*n)
  
  # gdp out
  tmp1 = matrix(log(gdp[included, i1]), n, n, byrow = FALSE)
  tmp2 = matrix(log(gdp[included, i1]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_out_layer1 = tmp1 ; log_gdp_out_layer2 = tmp2
  
  # gdp in
  tmp1 = matrix(log(gdp[included, i1]), n, n, byrow = TRUE)
  tmp2 = matrix(log(gdp[included, i1]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdp_in_layer1 = tmp1 ; log_gdp_in_layer2 = tmp2
  
  # atop alliance
  tmp1 = kronecker(A, atop_alliance[[i1]][included, included]) 
  tmp2 = kronecker(B, atop_alliance[[i1]][included, included])
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  alliance_layer1 = tmp1; alliance_layer2 = tmp2
  
  # absolute polity diff
  tmp1 = kronecker(A, abs(outer(polity[included, i1], polity[included, i1],'-')))
  tmp2 = kronecker(B, abs(outer(polity[included, i1], polity[included, i1],'-')))
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  polity_diff_layer1 = tmp1; polity_diff_layer2 = tmp2
  
  # military expenditure out
  tmp1 = matrix(log(milit_exp[included, i1]), n, n, byrow = FALSE)
  tmp2 = matrix(log(milit_exp[included, i1]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_out_layer1 = tmp1 ; log_milit_exp_out_layer2 = tmp2
  
  # military expenditure in
  tmp1 = matrix(log(milit_exp[included, i1]), n, n, byrow = TRUE)
  tmp2 = matrix(log(milit_exp[included, i1]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_in_layer1 = tmp1 ; log_milit_exp_in_layer2 = tmp2
  
  # define path dependency
  tmp1 = 1 * (sipri_tiv[[i1]][included, included] > 0)
  tmp2 = custom_trade_to_binary(trade[[i2]][included, included], type = type, threshold = 0.01)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  pathdep_layer1 = tmp1 ; pathdep_layer2 = tmp2
  
  
  # just use "free" observations
  #114*113 sanity check
  obs_total <- as.matrix.network(future)[which(free == 1, arr.ind = T)]
  obs_arms <- as.matrix.network(future)[which(free_arms == 1, arr.ind = T)]
  
  
  # predict t+1 without cross-layer dependence
  model_future = update(model, future ~ . - duplexdyad(c("e", "h"), layers = list(1, 2)))
  coef = coef(fit1)
  names(coef) = gsub("[[i]]", "", names(coef), fixed = T)
  
  
  pred_c_total <- predict(model_future, theta = coef, conditional = T, type = "response", 
                    output = "matrix", constraints = ~fixallbut(free))

  pred_c_total <- pred_c_total[which(free == 1, arr.ind = T)]
  
  
  #pred_c_arms <- predict(model_future, theta = coef, conditional = T, type = "response", 
                   # output = "matrix", constraints = ~fixallbut(free_arms))
  
  #pred_c_arms <- pred_c_arms[which(free == 1, arr.ind = T)]
  
  
  pred_total <- predict(model_future, theta = coef, conditional = F, type = "response", 
                  output = "matrix", nsim = n_pred, constraints = ~fixallbut(free))
  
  pred_total <- pred_total[which(free == 1, arr.ind = T)]
  
  
  pred_arms <- predict(model_future, theta = coef, conditional = F, type = "response", 
                  output = "matrix", nsim = n_pred, constraints = ~fixallbut(free_arms))
  
  pred_arms <- pred_arms[which(free_arms == 1, arr.ind = T)]
  
  

  fg <- pred_c_total[obs_total == 1]
  bg <- pred_c_total[obs_total == 0]

  pr_c_total <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  bs_c_total <- mean((pred_c_total - obs_total)**2)


  #fg <- pred_c_arms[obs_arms == 1]
  #bg <- pred_c_arms[obs_arms == 0]

  #pr_c_arms <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  #bs_c_arms <- mean((pred_c_arms - obs_arms)**2)


  fg <- pred_total[obs_total == 1]
  bg <- pred_total[obs_total == 0]

  pr_total <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  bs_total <- mean((pred_total - obs_total)**2)


  fg <- pred_arms[obs_arms == 1]
  bg <- pred_arms[obs_arms == 0]

  pr_arms <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  bs_arms <- mean((pred_arms - obs_arms)**2)
  
  
  scores = rbind(scores, c("restricted", year, "brier", "total (cond.)", bs_c_total))
  scores = rbind(scores, c("restricted", year, "brier", "total", bs_total))
  scores = rbind(scores, c("restricted", year, "brier", "arms", bs_arms))
  scores = rbind(scores, c("restricted", year, "pr", "total (cond.)", pr_c_total))
  scores = rbind(scores, c("restricted", year, "pr", "total", pr_total))
  scores = rbind(scores, c("restricted", year, "pr", "arms", pr_arms))
   
      
  # predict t+1 with cross-layer dependence
  model_future = update(model_future, ~ . + duplexdyad(c("e", "h"), layers = list(1, 2)))
  coef = coef(fit2)
  names(coef) = gsub("[[i]]", "", names(coef), fixed = T)
  
  
  pred_c_total <- predict(model_future, theta = coef, conditional = T, type = "response", 
                          output = "matrix", constraints = ~fixallbut(free))
  
  pred_c_total <- pred_c_total[which(free == 1, arr.ind = T)]
  
  
  #pred_c_arms <- predict(model_future, theta = coef, conditional = T, type = "response", 
                         #output = "matrix", constraints = ~fixallbut(free_arms))
  
  #pred_c_arms <- pred_c_arms[which(free == 1, arr.ind = T)]
  
  
  pred_total <- predict(model_future, theta = coef, conditional = F, type = "response", 
                        output = "matrix", nsim = n_pred, constraints = ~fixallbut(free))
  
  pred_total <- pred_total[which(free == 1, arr.ind = T)]
  
  
  pred_arms <- predict(model_future, theta = coef, conditional = F, type = "response", 
                       output = "matrix", nsim = n_pred, constraints = ~fixallbut(free_arms))
  
  pred_arms <- pred_arms[which(free_arms == 1, arr.ind = T)]
  
  
  fg <- pred_c_total[obs_total == 1]
  bg <- pred_c_total[obs_total == 0]
  
  pr_c_total <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  bs_c_total <- mean((pred_c_total - obs_total)**2)
  
  
  #fg <- pred_c_arms[obs_arms == 1]
  #bg <- pred_c_arms[obs_arms == 0]
  
  #pr_c_arms <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  #bs_c_arms <- mean((pred_c_arms - obs_arms)**2)
  
  
  fg <- pred_total[obs_total == 1]
  bg <- pred_total[obs_total == 0]
  
  pr_total <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  bs_total <- mean((pred_total - obs_total)**2)
  
  
  fg <- pred_arms[obs_arms == 1]
  bg <- pred_arms[obs_arms == 0]
  
  pr_arms <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = T)$auc.integral
  bs_arms <- mean((pred_arms - obs_arms)**2)
  
  
  scores = rbind(scores, c("full", year, "brier", "total (cond.)", bs_c_total))
  scores = rbind(scores, c("full", year, "brier", "total", bs_total))
  scores = rbind(scores, c("full", year, "brier", "arms", bs_arms))
  scores = rbind(scores, c("full", year, "pr", "total (cond.)", pr_c_total))
  scores = rbind(scores, c("full", year, "pr", "total", pr_total))
  scores = rbind(scores, c("full", year, "pr", "arms", pr_arms))
  
  print(paste0("Year ", year , " finished."))
}
scores$value <- as.numeric(scores$value)
scores$year <- strtoi(scores$year)
colnames(scores) <- c("model", "year", "metric", "type", "value")
save(scores, file = paste0(path, "/models/ERGM/", version, "_predictive_scores.RData"))

 


