library(tidyverse)
library(fs)
library(gt)
library(gtools)

source("utils/utils.R")

dpath <- data_path.get()
saom_multi <- readRDS(path(dpath, "models/SAOM/saom_multilevel_220412/final_fit_saom_multilevel_220412.rds"))
saom_simple <- readRDS(path(dpath, "models/SAOM/saom_simple_220412/final_fit_saom_simple_220412.rds"))

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
  fmt_missing(columns = everything(), rows = everything(), missing_text = "")
t1 %>%
  as_latex() %>% 
  as.character() %>% 
  cat(file = "figures/table_saom_base_multi.tex")
t1 
