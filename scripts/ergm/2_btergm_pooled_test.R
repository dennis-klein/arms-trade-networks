#------------------------------------------------------------------------------#
# Model estimation: only the year 2003
#------------------------------------------------------------------------------#
library(network)
library(ergm)
library(multilayer.ergm)
library(texreg)
library(kableExtra)
library(stargazer)

rm(list = ls(all.names = TRUE))


# setup
source("utils/utils.R")
source("utils/construct_header.R")
source("utils/custom_trade_to_binary.R")
path <- data_path.get()
set.seed(1234)


# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/cdist.RData"))
load(file.path(path, "out/atop_alliance.RData"))
load(file.path(path, "out/milit_exp.RData"))
nmc_cinc <- readRDS(file.path(path, "out/nmc_cinc.rds"))
polity <- readRDS(file.path(path, "out/polity.rds"))
gdp <- readRDS(file.path(path, "out/gdp.rds"))
gdppc <- readRDS(file.path(path, "out/gdppc.rds"))
sipri_tiv <- readRDS(file.path(path, "out/sipri_tiv.rds"))
trade <- readRDS(file.path(path, "out/baci_aggregated.rds"))



# trade data starting 1995
# A is if import AND export are above threshold
# "D" is export relevance, "C" is import relevance
# compare both results
start <- 1995
end <- 2018
maxlag <- 1


# setup standard errors
n_boot <- 100
n_sim <- 100


# selection of countries included
included <- rowSums(EX[, (start:end) - 1949]) == length(start:end)
n <- sum(included)


# some data transformations
# sipri military expenditure is in constant USD$m
# cinc is index which sums up to 0 -> better interpretability using pct points
milit_exp <- milit_exp * 1000000
nmc_cinc100 <- nmc_cinc * 100


# matrices for kronecker products
A <- matrix(c(1, 0, 0, 0), nrow = 2)
B <- matrix(c(0, 0, 0, 1), nrow = 2)


# setup constraints
free <- matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), nrow = 2 * n, byrow = T)
diag(free) <- 0


# init covariates
netC <- list()
netD <- list()
log_gdppc_out_layer1 <- list()
log_gdppc_out_layer2 <- list()
log_gdppc_in_layer1 <- list()
log_gdppc_in_layer2 <- list()
alliance_layer1 <- list()
alliance_layer2 <- list()
polity_diff_layer1 <- list()
polity_diff_layer2 <- list()
log_milit_exp_out_layer1 <- list()
log_milit_exp_out_layer2 <- list()
log_milit_exp_in_layer1 <- list()
log_milit_exp_in_layer2 <- list()
pathdep_layer1 <- list()
pathdep_layer2C <- list()
pathdep_layer2D <- list()


# distance
log_cdist_layer1 <- kronecker(A, log(cdist[included, included] + 1))
log_cdist_layer2 <- kronecker(B, log(cdist[included, included] + 1))
colnames(log_cdist_layer1) <- rownames(log_cdist_layer1) <- 1:(2 * n)
colnames(log_cdist_layer2) <- rownames(log_cdist_layer2) <- 1:(2 * n)


# fill matrices
for (year in (start + maxlag):end){
  
  # set correct indices corresponding data source
  i1 = year - 1949 # Cornelius Replication Data 1950:2018
  i2 = year - 1994 # CEPII BACI 1995:2018
  j = year - (start + 1 - 1)
  
  
  # dependent multi-layer network
  netC[[j]] = to.multiplex(
    1*(sipri_tiv[[i1]][included, included] > 0), 
    custom_trade_to_binary(trade[[i2]][included, included], type = "C", threshold = 0.01), 
    output = "network"
  )
  
  netD[[j]] = to.multiplex(
    1*(sipri_tiv[[i1]][included, included] > 0), 
    custom_trade_to_binary(trade[[i2]][included, included], type = "D", threshold = 0.01), 
    output = "network"
  )
 
  # include exogenous covariates lagged by 1 time period for each net
  # gdppc out
  tmp1 = matrix(log(gdppc[included, i1-1]), n, n, byrow = FALSE)
  tmp2 = matrix(log(gdppc[included, i1-1]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdppc_out_layer1[[j]] = tmp1 ; log_gdppc_out_layer2[[j]] = tmp2
  
  # gdppc in
  tmp1 = matrix(log(gdppc[included, i1-1]), n, n, byrow = TRUE)
  tmp2 = matrix(log(gdppc[included, i1-1]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_gdppc_in_layer1[[j]] = tmp1 ; log_gdppc_in_layer2[[j]] = tmp2
  
  # atop alliance
  tmp1 = kronecker(A, atop_alliance[[i1 - 1]][included, included]) 
  tmp2 = kronecker(B, atop_alliance[[i1 - 1]][included, included])
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  alliance_layer1[[j]] = tmp1; alliance_layer2[[j]] = tmp2
  
  # absolute polity diff
  tmp1 = kronecker(A, abs(outer(polity[included, i1-1], polity[included, i1-1],'-')))
  tmp2 = kronecker(B, abs(outer(polity[included, i1-1], polity[included, i1-1],'-')))
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  polity_diff_layer1[[j]] = tmp1; polity_diff_layer2[[j]] = tmp2
  
  # military expenditure out
  tmp1 = matrix(log(milit_exp[included, i1-1]), n, n, byrow = FALSE)
  tmp2 = matrix(log(milit_exp[included, i1-1]), n, n, byrow = FALSE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_out_layer1[[j]] = tmp1 ; log_milit_exp_out_layer2[[j]] = tmp2
  
  # military expenditure in
  tmp1 = matrix(log(milit_exp[included, i1-1]), n, n, byrow = TRUE)
  tmp2 = matrix(log(milit_exp[included, i1-1]), n, n, byrow = TRUE)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  log_milit_exp_in_layer1[[j]] = tmp1 ; log_milit_exp_in_layer2[[j]] = tmp2

  # define path dependency
  tmp1 = 1 * (sipri_tiv[[i1-1]][included, included] > 0)
  tmp2 = custom_trade_to_binary(trade[[i2-1]][included, included], type = "C", threshold = 0.01)
  tmp3 = custom_trade_to_binary(trade[[i2-1]][included, included], type = "D", threshold = 0.01)
  tmp1 = kronecker(A, tmp1) 
  tmp2 = kronecker(B, tmp2)
  tmp3 = kronecker(B, tmp3)
  colnames(tmp1) <- rownames(tmp1) <- 1:(2*n)
  colnames(tmp2) <- rownames(tmp2) <- 1:(2*n)
  colnames(tmp3) <- rownames(tmp3) <- 1:(2*n)
  pathdep_layer1[[j]] = tmp1 ; pathdep_layer2C[[j]] = tmp2; pathdep_layer2D[[j]] = tmp3
  
}


# specify models
modelC <- netC ~
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  mutual(same = "layer.mem", diff = TRUE) +
  gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
  gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
  duplexdyad(c("e", "f", "g", "h", "i"), layers = list(1, 2)) +
  edgecov(log_cdist_layer1) +
  edgecov(log_cdist_layer2) +
  edgecov(log_gdppc_in_layer1) +
  edgecov(log_gdppc_in_layer2) +
  edgecov(log_gdppc_out_layer1) +
  edgecov(log_gdppc_out_layer2) +
  edgecov(alliance_layer1) +
  edgecov(alliance_layer2) +
  edgecov(polity_diff_layer1) +
  edgecov(polity_diff_layer2) +
  edgecov(log_milit_exp_in_layer1) +
  edgecov(log_milit_exp_in_layer2) +
  edgecov(log_milit_exp_out_layer1) +
  edgecov(log_milit_exp_out_layer2) +
  edgecov(pathdep_layer1) +
  edgecov(pathdep_layer2C)

modelD <- netD ~
  edges_layer(layer = 1) +
  edges_layer(layer = 2) +
  mutual(same = "layer.mem", diff = TRUE) +
  gwidegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
  gwodegree(decay = 1, fixed = TRUE, attr = "layer.mem") +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 1) +
  gwesp_layer(decay = 1.5, fixed = TRUE, layer = 2) +
  duplexdyad(c("e", "f", "g", "h", "i"), layers = list(1, 2)) +
  edgecov(log_cdist_layer1) +
  edgecov(log_cdist_layer2) +
  edgecov(log_gdppc_in_layer1) +
  edgecov(log_gdppc_in_layer2) +
  edgecov(log_gdppc_out_layer1) +
  edgecov(log_gdppc_out_layer2) +
  edgecov(alliance_layer1) +
  edgecov(alliance_layer2) +
  edgecov(polity_diff_layer1) +
  edgecov(polity_diff_layer2) +
  edgecov(log_milit_exp_in_layer1) +
  edgecov(log_milit_exp_in_layer2) +
  edgecov(log_milit_exp_out_layer1) +
  edgecov(log_milit_exp_out_layer2) +
  edgecov(pathdep_layer1) +
  edgecov(pathdep_layer2D)


# estimate models
fit_C = btergm(modelC, parallel = "snow", ncpus = 2, verbose = T, R = n_boot)
fit_D = btergm(modelD, parallel = "snow", ncpus = 2, verbose = T, R = n_boot)

cat("Estimation of models finished. \n")
cat("Simulating networks. \n")


if (FALSE) {
  # simulate confidence intervals
  sim_C = simulate(fit_C, index = 1, nsim = n_sim, output = "stats", constraints = ~fixallbut(free), verbose = F)
  sim_D = simulate(fit_D, index = 1, nsim = n_sim, output = "stats", constraints = ~fixallbut(free), verbose = F)
  
  
  for (i in 2:23){
    i = as.numeric(i)
    tmp = simulate(fit_C, index = i, nsim = n_sim, output = "stats", constraints = ~fixallbut(free), verbose = F)
    sim_C = sim_C + tmp
    
    tmp = simulate(fit_D, index = i, nsim = n_sim, output = "stats", constraints = ~fixallbut(free), verbose = F)
    sim_D = sim_D + tmp
  }
  
  
  # compute covariance
  cov_C = cov(sim_C)
  cov_D = cov(sim_D)
}


#------------------------------------------------------------------------------#
# Goodness of Fit Assessment
#------------------------------------------------------------------------------#

# functions for gof
esp_layer1 <- function(mat, ...) {
  d <- summary(as.matrix(mat)[1:114, 1:114] ~ esp(0:(nrow(mat) - 2)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Edge-wise shared partners Layer 1"
  return(d)
}

esp_layer2 <- function(mat, ...) {
  d <- summary(as.matrix(mat)[115:228, 115:228] ~ esp(0:(nrow(mat) - 2)))
  names(d) <- 0:(length(d) - 1)
  attributes(d)$label <- "Edge-wise shared partners Layer 1"
  return(d)
}

geodesic_layer1 <- function(mat, ...) {
  
  fillup <- function(x, another.length) {  # fill up x if shorter
    difference <- length(x) - another.length
    inf.value <- x[length(x)]
    if (difference < 0) {  # x is shorter
      x <- x[1:(length(x) - 1)]
      x <- c(x, rep(0, abs(difference)), inf.value)
    } else if (difference > 0) {
      x <- x[1:(length(x) - difference)]
      x <- c(x, inf.value)
    }
    return(x)
  }
  
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat <- mat[1:114, 1:114]
  g <- fillup(ergm::ergm.geodistdist(network(as.matrix(mat), directed = F)), nrow(mat) - 1)
  attributes(g)$label <- "Geodesic distances Layer 1"
  return(g)
}

geodesic_layer2 <- function(mat, ...) {
  
  fillup <- function(x, another.length) {  # fill up x if shorter
    difference <- length(x) - another.length
    inf.value <- x[length(x)]
    if (difference < 0) {  # x is shorter
      x <- x[1:(length(x) - 1)]
      x <- c(x, rep(0, abs(difference)), inf.value)
    } else if (difference > 0) {
      x <- x[1:(length(x) - difference)]
      x <- c(x, inf.value)
    }
    return(x)
  }
  
  mat <- as.matrix(mat)
  mat[is.na(mat)] <- 0
  mat <- mat[115:228, 115:228]
  g <- fillup(ergm::ergm.geodistdist(network(as.matrix(mat), directed = F)), nrow(mat) - 1)
  attributes(g)$label <- "Geodesic distances Layer 2"
  return(g)
}


stats <- c(esp_layer1, esp_layer2, geodesic_layer1, geodesic_layer2)
gof_C <- gof(fit_C, statistics = stats, nsim = 100, parallel = "snow", ncpus = 2, verbose = T)
gof_D <- gof(fit_D, statistics = stats, nsim = 100, parallel = "snow", ncpus = 2, verbose = T)



#------------------------------------------------------------------------------#
# Save results
#------------------------------------------------------------------------------#

# save
save.image(file = paste0(path, "/models/ergm/2_tergm_pooled.RData"))

# load
#load(file = paste0(path, "/models/ergm/2_tergm_pooled.RData"))




#------------------------------------------------------------------------------#
# Plot tables for report
#------------------------------------------------------------------------------#

custom.coef.names <- c(
  "Weapons: Edges",
  "Conventional: Edges",
  "Weapons: Mutuality",
  "Conventional: Mutuality",
  "Weapons: GWIDEG.1",
  "Conventional: GWIDEG.1",
  "Weapons: GWODEG.1",
  "Conventional: GWODEG.1",
  "Weapons: GWESP.1.5",
  "Conventional: GWESP.1.5",
  "E",
  "F",
  "G",
  "H",
  "I",
  "Weapons: log Distance",
  "Conventional: log Distance",
  "Weapons: log GDP per capita in",
  "Conventional: log GDP per capita in",
  "Weapons: log GDP per capita out",
  "Conventional: log GDP per capita out",
  "Weapons: Alliance",
  "Conventional: Alliance",
  "Weapons: Polity Diff. (abs)",
  "Conventional: Polity Diff. (abs)",
  "Weapons: log Military Exp. in",
  "Conventional: log Military Exp. in",
  "Weapons: log Military Exp. out",
  "Conventional: log Military Exp. out",
  "Weapons: Path Dependency",
  "Conventional: Path Dependency"
)


screenreg(
  list(fit_res, fit_full), 
  caption = "Multilayer temporal ERG model for the period 1996-2018. Maximum pseudolikelihood estimates and 95 pct. confidence intervals constructed using 1000 bootstrap iterations for each year.",
  custom.model.names = c("Restricted","Full"),
  custom.coef.names = custom.coef.names, 
  label = "table:model1", 
  groups = list("Network Effects" = 1:6, "Covariates" = 7:18, "Cross-Layer Network Effects" = 19:20),
  dcolumn = TRUE, single.row = T, fontsize = "small", include.nobs = F,
  float.pos = "!htp"
)


if(FALSE){
  texreg(
    list(fit_res, fit_full), 
    caption = "Multilayer temporal ERG model for the period 1996-2018. Maximum pseudolikelihood estimates and 95 pct. confidence intervals constructed using 1000 bootstrap iterations for each year.",
    custom.model.names = c("Restricted","Full"),
    custom.coef.names = custom.coef.names, 
    label = "table:ergm1", 
    groups = list("Network Effects" = 1:6, "Covariates" = 7:22, "Cross-Layer Network Effects" = 23:24),
    dcolumn = TRUE, single.row = T, fontsize = "small", include.nobs = F,
    float.pos = "!htp"
  )
}

# output with kableExtra
df <- data.frame(coef = custom.coef.names)
df$res_est <- c(fit_res@coef, NA, NA)
df$res_lcib <- c(confint(fit_res)[, 2], NA, NA)
df$res_ucib <- c(confint(fit_res)[, 3], NA, NA)
df$res_lci <- df$res_est - c(qnorm(0.975) * sqrt(diag(cov_res)), NA, NA)
df$res_uci <- df$res_est + c(qnorm(0.975) * sqrt(diag(cov_res)), NA, NA)

df$full_est <- c(fit_full@coef)
df$full_lcib <- c(confint(fit_full)[, 2])
df$full_ucib <- c(confint(fit_full)[, 3])
df$full_lci <- df$full_est - qnorm(0.975) * sqrt(diag(cov_full))
df$full_uci <- df$full_est + qnorm(0.975) * sqrt(diag(cov_full))

df[, 2:ncol(df)] <- round(df[, 2:ncol(df)], 3)
df[is.na(df)] <- "-"


# add indicator for statistical significance
ind <- !(df$res_lcib <= 0 & 0 <= df$res_ucib)
ind[19:20] <- FALSE
df$res_est[ind] <- paste0(df$res_est[ind], "*")
df$res_est[!ind] <- paste0(df$res_est[!ind], " ")

ind <- !(df$full_lcib <= 0 & 0 <= df$full_ucib)
df$full_est[ind] <- paste0(df$full_est[ind], "*")
df$full_est[!ind] <- paste0(df$full_est[!ind], " ")
names <- c("", rep(c("estimate", "lower", "upper"), 2))


# output latex table
sink(file = "figures/ergm_model_2.txt")
kable(df[, c(1, 2, 3, 4, 7, 8, 9)], 
      format = "latex", booktabs = T, label = "ergm_model_1",
      col.names = names, align = "lrrrrrr", 
      caption = "Multilayer temporal ERG model (directed) for the period 1996-2018.") %>%
  kable_styling(latex_options = c("scale_down")) %>%
  add_header_above(c(" ", "Import" = 3, "Export" = 3)) %>%
  pack_rows("Network Effects", 1, 6) %>%
  pack_rows("Covariates", 7, 18) %>%
  pack_rows("Cross-Layer Network Effects", 19, 20) %>%
  footnote(general = "* Null hypothesis value outside the confidence interval. Maximum pseudolikelihood estimates and 95 pct. confidence intervals constructed using 1000 bootstrap iterations for each year. Cross-layer statistics as defined in \ref{fig:crosslayer}.", 
           threeparttable = T)
sink()





#------------------------------------------------------------------------------#
# Plot gof stats for report
#------------------------------------------------------------------------------#

pdf(file = "figures/ergm_gof_2.pdf", width = 10, height = 0.9 * sqrt(2) * 10)
par(mfrow = c(3,2))
plot(gof_full)
dev.off()

