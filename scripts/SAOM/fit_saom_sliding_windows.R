# Fitting a SAOM sliding windows model

require(RSiena)
require(fs)
source("scripts/SAOM/siena07ToConverge.R")

fit_saom_sliding_windows <- function(win_size, model_id, data_obj_path, save_dir,
                                     test_mode = FALSE) {
  if (dir.exists(path(save_dir, model_id))) stop("Model already exsists")
  
  # setup
  load(data_obj_path)
  model_dir <- path(save_dir, model_id)
  dir.create(model_dir)
  models <- list()
  
  # dimensions
  actors <- dimnames(arm[[1]])[[1]] # actors
  act <- 1:length(actors)
  obs <- length(arm) # observations
  windows <- lapply(1:(obs-win_size+1), function(x) x:(x+win_size-1))
  if (obs < win_size) stop("Too few observations for window")
  
  # test_mode
  if (test_mode) {
    message("TEST MODE")
    windows <- sample(windows, 2)
    act_sub <- c(sample(act, 10), c(4, 8, 21, 31, 36, 39, 47, 89, 91, 108, 109, 110))
    actors <- actors[act_sub]
  }
  
  # save some model information
  model_meta <- list(
    actors = actors,
    windows = windows,
    years = names(arm)
  )
  saveRDS(model_meta, path(model_dir, "model_meta", ext = "rds"))
  
  i <- 1
  for (win in windows) {
    m <- sprintf("Fitting window %d of %d... window: %s",
                 i, length(windows), paste(win, collapse = " "))
    message(m)
    
    # set up RSiena data for window
    dyad_win <- head(win, -1)
    RSiena_vars <- list(
      # Dependent variables
      arm = sienaNet(lapply(arm[win], function(x) x[act, act])),
      trd = sienaNet(lapply(trd[win], function(x) x[act, act])),
      
      # Constant covariates
      # -
      
      # Varying covariates
      gdp_log = varCovar(gdp_log[act, win]),
      nmc_std = varCovar(nmc_std[act, win]),
      
      # Constant dyadic covariates
      cdist_std = coDyadCovar(cdist_std[act, act]),
      
      # Varying dyadic covariates
      pol_diff_std = varDyadCovar(lapply(pol_diff_std[dyad_win], function(x) x[act, act])),
      allied = varDyadCovar(lapply(allied[dyad_win], function(x) x[act, act])),
      trd_lag1 = varDyadCovar(lapply(trd_lag1[dyad_win], function(x) x[act, act])),
      trd_lag2 = varDyadCovar(lapply(trd_lag2[dyad_win], function(x) x[act, act])),
      trd_lag3 = varDyadCovar(lapply(trd_lag3[dyad_win], function(x) x[act, act])),
      arm_lag1 = varDyadCovar(lapply(arm_lag1[dyad_win], function(x) x[act, act])),
      arm_lag2 = varDyadCovar(lapply(arm_lag2[dyad_win], function(x) x[act, act])),
      arm_lag3 = varDyadCovar(lapply(arm_lag3[dyad_win], function(x) x[act, act]))
    )
    
    dat <- do.call(sienaDataCreate, RSiena_vars)
    
    # define effects
    eff <- getEffects(dat)
    ## within networks effects
    # standard effects (outdegree, reciprocity automatically added)
    eff <- includeEffects(eff, transTies, name = "arm", verbose = FALSE)
    eff <- includeEffects(eff, transTies, name = "trd", verbose = FALSE)
    # covariate effects
    eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "gdp_log", verbose = FALSE)
    eff <- includeEffects(eff, egoX, altX, name = "arm", interaction1 = "nmc_std", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "arm", interaction1 = "pol_diff_std", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "arm", interaction1 = "cdist_std", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "arm", interaction1 = "allied", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "arm", interaction1 = "arm_lag1", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "arm", interaction1 = "arm_lag2", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "arm", interaction1 = "arm_lag3", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "arm", interaction1 = "trd_lag1", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "arm", interaction1 = "trd_lag2", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "arm", interaction1 = "trd_lag3", verbose = FALSE)
    eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "gdp_log", verbose = FALSE)
    eff <- includeEffects(eff, egoX, altX, name = "trd", interaction1 = "nmc_std", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "trd", interaction1 = "pol_diff_std", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "trd", interaction1 = "cdist_std", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "trd", interaction1 = "allied", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "trd", interaction1 = "arm_lag1", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "trd", interaction1 = "arm_lag2", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "trd", interaction1 = "arm_lag3", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "trd", interaction1 = "trd_lag1", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "trd", interaction1 = "trd_lag2", verbose = FALSE)
    eff <- includeEffects(eff, X, name = "trd", interaction1 = "trd_lag3", verbose = FALSE)
    ## between networks effects
    eff <- includeEffects(eff, crprod, name = "arm", interaction1 = "trd", verbose = FALSE)
    eff <- includeEffects(eff, crprod, name = "trd", interaction1 = "arm", verbose = FALSE)
    eff <- includeEffects(eff, crprodRecip, name = "arm", interaction1 = "trd", verbose = FALSE)
    eff <- includeEffects(eff, crprodRecip, name = "trd", interaction1 = "arm", verbose = FALSE)
    
    # test mode effects
    if (test_mode) {
      eff <- getEffects(dat)
    }
    
    # run model
    alg <- sienaAlgorithmCreate(projname = NULL)
    ans_id <- sprintf("%03d_to_%03d", win[1], win[length(win)])
    ans_dir <- path(model_dir, paste0("tmp_fit_", ans_id))
    dir.create(ans_dir)
    if (!exists("ans")) ans <- NULL
    ans <- siena07ToConvergence(alg = alg, dat = dat, eff = eff,
                                ans0 = ans,
                                save_dir = ans_dir, ans_id = ans_id)
    # save SAOM
    models[[ans_id]] <- ans
    saveRDS(ans, file = path(model_dir, paste0("fit_", ans_id), ext = "rds"))
    
    # increment counter
    i <- i + 1
  }
  saveRDS(models, file = path(model_dir, "model", ext = "rds"))
}
