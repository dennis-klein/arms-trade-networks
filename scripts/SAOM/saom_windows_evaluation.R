library(RSiena)
library(ggplot2)

source("utils/utils.R")


dpath <- data_path.get()
model_folder <- "test_tmp_1647432770"
model <- readRDS(path(dpath, "models/SAOM", model_folder, "model.rds"))
model_meta <- readRDS(path(dpath, "models/SAOM", model_folder, "model_meta.rds"))
plots <- list()

tmp <- models[[1]]

tmp$theta
tmp$effects$effectName
tmp$se

dfs <- lapply(model, function(x) data.frame(coef = x$effects$effectName,
                                     value = x$theta,
                                     se = tmp$se))
window_years <- lapply(model_meta$windows, function(x) paste0(model_meta$years[x][1],
                                                             ":",
                                                             model_meta$years[x][length(model_meta$years[x])]))
for (i in 1:length(dfs)) {
  dfs[[i]]$years <- window_years[[i]]
}

df <- do.call(rbind, dfs)

ggplot(df, aes(x=years, y=value)) +
  geom_point() +
  geom_errorbar(aes(ymin = value-1.96*se, ymax = value+1.96*se), width = 0.1) +
  geom_hline(aes(yintercept=0), color="blue") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=90)) +
  facet_wrap(~ coef, scales = "free_y")
