library(tidyverse)
library(fs)
library(gt)
library(gtools)
library(RSiena)
library(parallel)
library(patchwork)
library(ggplot2)
library(gridExtra)

source("utils/utils.R")


### Results table

dpath <- data_path.get()
saom_multi <- readRDS(path(dpath, "models/main-saom-models/saom_multilevel_220510_import/final_fit_saom_multilevel_220510_import.rds"))
saom_simple <- readRDS(path(dpath, "models/main-saom-models/saom_simple_220511_import/final_fit_saom_simple_220511_import.rds"))

conv_multi <- saom_multi$tconv.max[1,1]
all_t_multi <- max(saom_multi$tconv)

conv_simple <- saom_simple$tconv.max[1,1]
all_t_simple <- max(saom_simple$tconv)


eff_multi <- saom_multi$effects
eff_multi$estimate <- saom_multi$theta
eff_multi$se <- saom_multi$se
eff_multi$model <- "multi"

eff_simple <- saom_simple$effects
eff_simple$estimate <- saom_simple$theta
eff_simple$se <- saom_simple$se
eff_simple$model <- "simple"

eff <- rbind(eff_multi, eff_simple)

eff$test_stat <- abs(eff$estimate / eff$se)
eff$p_val <- pnorm(eff$test_stat, lower.tail = FALSE)*2
eff$p_stars <- str_pad(stars.pval(eff$p_val), width = 3, side = "right")
eff$entry <- sprintf("% 2.2f (%2.2f)%s", eff$estimate, eff$se, eff$p_stars)

struc_between <- c("crprod", "crprodRecip")
struc_within <- c("density", "recip", "transRecTrip", "cycle3", "gwespFF", "inPopSqrt", "inActSqrt", "outActSqrt")
cov_eff <- c("X", "egoX", "altX")

eff$effect_type <- case_when(
  eff$shortName %in% struc_within ~ "Within",
  eff$shortName %in% struc_between ~ "Between",
  eff$shortName %in% cov_eff ~ "Covariates"
) %>% 
  factor(levels = c("Covariates", "Within", "Between"))


eff$effect_identifier <- str_c(eff$shortName, eff$interaction1)
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
colnames(effect_key) <- c("id", "display")
eff$effect_display <- plyr::mapvalues(
  eff$effect_identifier,
  from = effect_key$id,
  to = effect_key$display
) %>% 
  factor(levels = effect_key$display)

# eff$net_display <- plyr::mapvalues(
#   eff$name,
#   from = c("arm", "trd"),
#   to = c("Arms", "Trade")
# )

# eff %>%
#   tibble %>% 
#   filter(type == "eval") %>%
#   arrange(effect_type) %>% 
#   select(name, effect_type, effect_display, entry) %>% 
#   pivot_wider(names_from = net_display, values_from = entry)

t1 <- eff %>%
  tibble %>% 
  filter(type == "eval") %>%
  arrange(effect_type) %>% 
  select(name, effect_type, effect_display, entry, model) %>% 
  pivot_wider(names_from = model, values_from = entry, names_prefix = "entry_") %>% 
  pivot_wider(names_from = name, values_from = c(entry_multi, entry_simple)) %>% 
  gt(
    rowname_col = "effect_display",
    groupname_col = "effect_type"
  ) %>% 
  tab_spanner(
    label = "Restricted Model",
    columns = c(entry_simple_arm, entry_simple_trd)
  ) %>% 
  tab_spanner(
    label = "With Between Effects",
    columns = c(entry_multi_arm, entry_multi_trd)
  ) %>% 
  cols_label(
    entry_multi_arm = "Arms", entry_multi_trd = "Trade",
    entry_simple_arm = "Arms", entry_simple_trd = "Trade"
  ) %>% 
  cols_align(
    align = "center",
    columns = c(entry_multi_arm, entry_multi_trd, entry_simple_arm, entry_simple_trd)
  ) %>% 
  cols_move_to_start(
    columns = c(entry_simple_arm, entry_simple_trd, entry_multi_arm, entry_multi_trd)
  ) %>% 
  tab_stubhead(label = "Effect") %>% 
  fmt_missing(columns = everything(), rows = everything(), missing_text = "") %>% 
  tab_header(
    title = "Estimated Coeffcients for the Stochastic Actor-Oriented Models",
    subtitle = "Main model and restricted model without between networks effects"
  ) %>% 
  tab_footnote(
    footnote = sprintf("Overall max. conv. ratio: %.2f; all t-conv. ratios below %.2f.",
                       conv_simple, all_t_simple),
    locations = cells_column_spanners(
      spanners = c("Restricted Model")
    )
  ) %>% 
  tab_footnote(
    footnote = sprintf("Overall max. conv. ratio: %.2f; all t-conv. ratios below %.2f.",
                       conv_multi, all_t_multi),
    locations = cells_column_spanners(
      spanners = c("With Between Effects")
    )
  )
t1
t1 %>%
  as_latex() %>% 
  as.character() %>% 
  cat(file = "figures/table_saom_base_multi_import.tex")

# TODO get latex table header center

### GOFs

dat_simple <- readRDS(file = path(dpath, "models/main-saom-models/saom_simple_220511_import/data_object_saom_simple_220511_import.rds"))
eff_simple <- readRDS(file = path(dpath, "models/main-saom-models/saom_simple_220511_import/effects_object_saom_simple_220511_import.rds"))
eff_simple$initialValue[eff_simple$include] <- saom_simple$theta

dat_multi <- readRDS(file = path(dpath, "models/main-saom-models/saom_multilevel_220510_import/data_object_saom_multilevel_220510_import.rds"))
eff_multi <- readRDS(file = path(dpath, "models/main-saom-models/saom_multilevel_220510_import/effects_object_saom_multilevel_220510_import.rds"))
eff_multi$initialValue[eff_multi$include] <- saom_multi$theta

# simulations, multicore
# n.clus <- detectCores() - 1
# 
# alg_sim <- sienaAlgorithmCreate(projname = "alg_sim", cond = FALSE,
#                                 useStdInits = FALSE, nsub = 0, simOnly = TRUE, n3 = 100)
# 
# sims_simple <- siena07(alg_sim, data = dat_simple, effects = eff_simple, batch = TRUE,
#                        returnDeps = T,
#                        useCluster = TRUE, nbrNodes = n.clus, initC = TRUE)
# sims_multi <- siena07(alg_sim, data = dat_multi, effects = eff_multi, batch = TRUE,
#                        returnDeps = T,
#                       useCluster = TRUE, nbrNodes = n.clus)
# 
# saveRDS(sims_simple, file = path(dpath, "models/main-saom-models/saom_simple_220511_import/sims_saom_simple_220511_import.rds"))
# saveRDS(sims_multi, file = path(dpath, "models/main-saom-models/saom_multilevel_220510_import/sims_saom_multilevel_220510_import.rds"))

sims_simple <- readRDS(file = path(dpath, "models/main-saom-models/saom_simple_220511_import/sims_saom_simple_220511_import.rds"))
sims_multi <- readRDS(file = path(dpath, "models/main-saom-models/saom_multilevel_220510_import/sims_saom_multilevel_220510_import.rds"))


stats <- c(OutdegreeDistribution, IndegreeDistribution, TriadCensus)
dep_vars <- c("arm", "trd")
gofs <- expand_grid(stats, dep_vars)


# gof_out_arm <- sienaGOF(sims_multi, OutdegreeDistribution, verbose = T, join = T, varName = "arm")
# gof_out_trd <- sienaGOF(sims_multi, OutdegreeDistribution, verbose = T, join = T, varName = "trd")
# gof_in_arm <- sienaGOF(sims_multi, IndegreeDistribution, verbose = T, join = T, varName = "arm")
# gof_in_trd <- sienaGOF(sims_multi, IndegreeDistribution, verbose = T, join = T, varName = "trd")
# gof_tri_arm <- sienaGOF(sims_multi, TriadCensus, verbose = T, join = T, varName = "arm")
# gof_tri_trd <- sienaGOF(sims_multi, TriadCensus, verbose = T, join = T, varName = "trd")
# 
# save(gof_out_arm, gof_out_trd,
#      gof_in_arm, gof_in_trd,
#      gof_tri_arm, gof_tri_trd,
#      file = path(dpath, "models/main-saom-models/saom_multilevel_220510_import/gofs.RData"))

load(file = path(dpath, "models/main-saom-models/saom_multilevel_220510_import/gofs.RData"))

p1 <- plot(gof_out_arm)
p2 <- plot(gof_out_trd)
p3 <- plot(gof_in_arm)
p4 <- plot(gof_in_trd)
p5 <- plot(gof_tri_arm)
p6 <- plot(gof_tri_trd)

# ggsave("figures/sliding_windows_netstruc.pdf", plot = p2, width = 8.27, height = 10.69)

p_slide <- grid.arrange(p1, p3, p5, p2, p4, p6, nrow = 2)
ggsave("figures/gofs_multi_slide.pdf", plot = p_slide, width = 5*2.5, height = 3.8*2.5)

p_report <- grid.arrange(p1, p2, p3, p4, p5, p6, nrow = 3)
ggsave("figures/gofs_multi_report.pdf", plot = p_report, width = 1.5*8.27, height = 1.5*10.69)
