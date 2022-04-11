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

# build data frame with estimates



res <- data.frame()
for (i in 1:length(models)) {
  mod <- models[[i]]
  df <- mod$effects
  df$est <- mod$theta
  df$se <- mod$se
  df$win <- i
  df$ci_low <- df$est - 1.96*df$se
  df$ci_up <- df$est + 1.96*df$se
  res <- rbind(res, df)
}


p1 <- res %>% 
  filter(functionType == "objective") %>% 
  # filter(name == "arm") %>% 
  ggplot(aes(x = win + 1998 - 1, y = est)) +
  geom_hline(aes(yintercept=0), color="#616161") +
  geom_point() +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_up), width = 0.1) +
  geom_line(stat="smooth", method = "loess", formula = y ~ x,
            size = 0.7, color = "blue", alpha = 0.7) +
  facet_wrap(effectName ~ ., ncol = 2, scales = "free") +
  labs(title = "Estimated Coefficients in Sliding Windows Model",
       caption = "4 year windows, data from 1998 to 2018") +
  xlab("Window starting year") +
  ylab("Estimate") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=0))
ggsave("figures/plot_saom_sw.pdf", plot = p1, width = 7, height = 25)

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
  facet_wrap(~ coef, scales = "free_y") +
