library(RSiena)
library(tidyverse)
library(fs)
library(patchwork)
library(lemon)

source("utils/utils.R")

dpath <- data_path.get()
model_id <- "model_220407"
model_path <- path(dpath, "models", "SAOM", model_id)
model_path <- path(dpath, "models/main-saom-models/saom_sw_220520_import")
files <- list.files(path = model_path, pattern = "final_fit")
windows <- lapply(files, function(x) str_match(x, pattern = "_(\\d+)\\.")[2])
models <- lapply(files, function(x) readRDS(file = path(model_path, x)))

models[[1]]$tconv.max

res <- data.frame()
for (i in 1:length(models)) {
  mod <- models[[i]]
  df <- mod$effects
  df$est <- mod$theta
  df$se <- mod$se
  df$win_number <- i
  df$ci_low <- df$est - 1.96*df$se
  df$ci_up <- df$est + 1.96*df$se
  res <- rbind(res, df)
}
res <- tibble(res)
win_size <- models[[1]]$observations
t_start <- 1998
res$win_start <- res$win_number - 1 + t_start
res$win_end <- res$win_start + win_size - 1
res$win <- sprintf("%i-%i", res$win_start, res$win_end)

tmp <- list(
  c("Rate", "Rate"),
  c("density", "Outdegree"),
  c("recip", "Reciprocity"),
  c("transRecTrip", "Transit. Recip. Triplets"),
  c("cycle3", "3-cycles"),
  c("gwespFF", "GWESP"),
  c("inPopSqrt", "Indegree Popularity (sqrt)"),
  c("inActSqrt", "Indegree Activity (sqrt)"),
  c("outActSqrt", "Outdegree Activity (sqrt)"),
  c("Xcdist", "Capital Distance"),
  c("Xpol_diff", "Abs. Polity Difference"),
  c("Xallied", "Alliance"),
  c("altXgdp_log", "GDP Alter (log)"),
  c("egoXgdp_log", "GDP Ego (log)"),
  c("altXmil_log", "Military Exp. Alter (log)"),
  c("egoXmil_log", "Military Exp. Ego (log)"),
  c("crprodtrd", "Conventional Trade"),
  c("crprodReciptrd", "Recip. Conventional Trade"),
  c("crprodarm", "Arms Trade"),
  c("crprodReciparm", "Recip. Arms Trade")
)
effect_key <- data.frame(do.call(rbind, tmp))
res$effect_identifier <- str_c(res$shortName, res$interaction1)
colnames(effect_key) <- c("id", "display")
res$effect_display <- plyr::mapvalues(
  res$effect_identifier,
  from = effect_key$id,
  to = effect_key$display
) %>% 
  factor(levels = effect_key$display)

struc_between <- c("crprod", "crprodRecip")
struc_within <- c("density", "recip", "transRecTrip", "cycle3", "gwespFF", "inPopSqrt", "inActSqrt", "outActSqrt")
cov_eff <- c("X", "egoX", "altX")

res$effect_type <- case_when(
  res$shortName %in% struc_within ~ "Within",
  res$shortName %in% struc_between ~ "Between",
  res$shortName %in% cov_eff ~ "Covariates"
) %>% 
  factor(levels = c("Covariates", "Within", "Between"))


plot_over_time <- function(net_filter, effect_selection, ncol, plot_title) {
  res %>% 
    filter(name == net_filter) %>% 
    filter(effect_identifier %in% effect_selection) %>% 
    ggplot(aes(x = win_start, y = est)) +
    geom_hline(aes(yintercept=0), color="#616161") +
    geom_point() +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_up), width = 0) +
    geom_line(aes(x = win_start, y = est), stat="smooth", method = "loess",
              formula = y ~ x, size = 0.7, color = "blue", alpha = 0.7) +
    facet_wrap(effect_display ~ ., ncol = ncol, scales = "free_y") +
    labs(title = plot_title) +
    xlab("") +
    ylab("") +
    scale_y_symmetric(mid = 0) +
    theme_minimal()
}

# Covariates Effects
p_arm_act <- plot_over_time("arm", c("egoXgdp_log", "altXgdp_log",
                                     "egoXmil_log", "altXmil_log"),
                            ncol = 4, plot_title = "Arms, Actor Covariates")
p_trd_act <- plot_over_time("trd", c("egoXgdp_log", "altXgdp_log",
                                     "egoXmil_log", "altXmil_log"),
                            ncol = 4, plot_title = "Trade, Actor Covariates")
p_arm_dya <- plot_over_time("arm", c("Xcdist", "Xpol_diff", "Xallied"), ncol = 3,
                            plot_title = "Arms, Dyadic Covariates")
p_trd_dya <- plot_over_time("trd", c("Xcdist", "Xpol_diff", "Xallied"), ncol = 3,
                            plot_title = "Trade, Dyadic Covariates")

p1 <- p_arm_act / p_trd_act / p_arm_dya / p_trd_dya +
  plot_annotation(tag_levels = "A")
# plot size dinA4 - 1 inch each side
ggsave("figures/sliding_windows_cov.pdf", plot = p1, width = 8.27, height = 10.69)

p1_slide <- p_trd_act / p_trd_dya +
  plot_annotation(tag_levels = "A")
ggsave("figures/sliding_windows_cov_slide.pdf", plot = p1_slide, width = 5*1.5, height = 3.8*1.5)



# Network structure Effects
p_arm_within <- plot_over_time("arm", c("density", "recip", "transRecTrip",
                                        "cycle3", "gwespFF", "inPopSqrt",
                                        "inActSqrt", "outActSqrt"),
                               ncol = 4, plot_title = "Arms, Within Effects")
p_trd_within <- plot_over_time("trd", c("density", "recip", "transRecTrip",
                                        "cycle3", "gwespFF", "inPopSqrt",
                                        "inActSqrt", "outActSqrt"),
                               ncol = 4, plot_title = "Trade, Within Effects")
p_arm_between <- plot_over_time("arm", c("crprodtrd", "crprodReciptrd"),
                               ncol = 4, plot_title = "Arms, Between Effects")
p_trd_between <- plot_over_time("trd", c("crprodarm", "crprodReciparm"),
                                ncol = 4, plot_title = "Trade, Between Effects")
p2 <- p_arm_within / p_trd_within / ( p_arm_between | p_trd_between ) +
  plot_annotation(tag_levels = "A") +
  plot_layout(heights = c(2,2,1))
ggsave("figures/sliding_windows_netstruc.pdf", plot = p2, width = 8.27, height = 10.69)


p2_slide <- (p_trd_within / p_trd_between) +
  plot_annotation(tag_levels = "A")
  plot_layout(widths = c(1.2,1))
ggsave("figures/sliding_windows_netstruc_slide.pdf", plot = p2_slide, width = 5*1.5, height = 3.8*1.5)

