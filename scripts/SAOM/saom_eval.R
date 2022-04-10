library(tidyverse)
library(RSiena)
library(gt)
library(fs)
source("utils/utils.R")

dpath <- data_path.get()
saom <- readRDS(path(dpath, "models/SAOM/model_220331/tmp_fit_model_220331_1", ext = "rds"))
tt <- sienaTimeTest(saom)

eff_no_rate <- saom$effects %>% 
  filter(!str_detect(effectName, "rate"))

eff <- saom$effects
eff$estimate <- saom$theta
eff$se <- saom$se

struc_between <- c("crprod", "crprodRecip")
struc_within <- c("density", "recip", "transRecTrip", "cycle3", "gwespFF", "inPopSqrt", "inActSqrt", "outActSqrt")
cov_eff <- c("X", "egoX", "altX")

eff$struc_within <- eff$shortName %in% struc_within
eff$struc_between <- eff$shortName %in% struc_between
eff$cov_eff <- eff$shortName %in% cov_eff
eff$test_stat <- eff$estimate / eff$se
eff$sig5 <- abs(eff$test_stat) > qnorm(0.975)
eff$sig1 <- abs(eff$test_stat) > qnorm(0.995)
eff$sig01 <- abs(eff$test_stat) > qnorm(0.9995)

eff %>%
  tibble %>% 
  filter(!str_detect(effectName, "rate")) %>% 
  select(name, shortName, estimate, se, interaction1) %>% 
  pivot_wider(names_from = name, values_from = c("estimate", "se")) %>% 
  View()

eff %>%
  tibble %>% 
  filter(!str_detect(effectName, "rate")) %>% 
  mutate(shortNameInter1 = paste(shortName, interaction1, sep = "_")) %>% 
  select(name, estimate, se, shortNameInter1) %>% 
  pivot_wider(names_from = name, values_from = c("estimate", "se")) %>% 
  gt(rowname_col = "shortNameInter1") %>% 
  tab_row_group(
    label = "Covariates",
    rows = c("egoX_gdp_log", "altX_gdp_log", "egoX_mil_log", "altX_mil_log",
             "X_cdist", "X_pol_diff", "X_allied")
  ) %>% 
  tab_row_group(
    label = "Within",
    rows = paste0(struc_within, "_")
  ) %>% 
  tab_row_group(
    label = "Between",
    rows = c("crprod_trd", "crprodRecip_trd", "crprod_arm", "crprodRecip_arm")
  ) %>% 
  tab_stubhead(label = "Effect") %>% 
  tab_spanner(
    label = "Arms Trade",
    columns = c(estimate_arm, se_arm)
  ) %>% 
  tab_spanner(
    label = "Conventional Trade",
    columns = c(estimate_trd, se_trd)
  ) %>% 
  cols_label(
    estimate_arm = "Estimate",
    se_arm = "SE"
  ) %>% 
  fmt_number(
    columns = c(estimate_trd, se_trd, estimate_arm, se_arm)
  )


eff %>% 
  filter(name == "arm") %>% 
  filter(struc_within) %>% 
  ggplot(aes(x=effectName, y=estimate)) +
  geom_point() +
  coord_flip()



eff %>% 
  filter(!str_detect(effectName, "rate")) %>% 
  mutate(est_std = estimate / se,
         arm_eff = str_detect(effectName ,"arm:")) %>% 
  ggplot(aes(x=effectName, y=est_std)) +
  geom_point() +
  coord_flip() +
  facet_grid(arm_eff ~ .) +
  theme_bw() +
  facet_grid(arm_eff ~ .)

xtable(saom, type = "txt")

# struc_within <- c(
#   "outdegree (density)",
#   "reciprocity",
#   "transitive recipr. triplets",
#   "3-cycles",
#   "GWESP I -> K -> J (69)",
#   "indegree - popularity (sqrt)",
#   "indegree - activity (sqrt)",
#   "outdegree - activity (sqrt)"
# )
# cov_eff <- c(
#   "cdist",
#   "pol_diff",
#   "allied",
#   "gdp_log alter",
#   "gdp_log ego",
#   "mil_log alter",
#   "mil_log ego",
# )
# struc_between <- c(
#   "trd", "reciprocity with trd ", "arm", "reciprocity with arm"
# )