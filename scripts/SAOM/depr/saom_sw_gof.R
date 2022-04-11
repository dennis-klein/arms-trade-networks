library(RSiena)
library(tidyverse)
library(fs)

source("utils/utils.R")

dpath <- data_path.get()
model_id <- "model_220407"
model_path <- path(dpath, "models", "SAOM", model_id)
files <- list.files(path = model_path, pattern = "final_fit")
windows <- lapply(files, function(x) str_match(x, pattern = "_(\\d+)\\.")[2])
models <- lapply(files, function(x) readRDS(file = path(model_path, x)))

tmp <- models[[1]]

sim_models <- sienaAlgorithmCreate()