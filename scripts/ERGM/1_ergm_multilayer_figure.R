#------------------------------------------------------------------------------#
# Replication code for the ERGM Analysis
# MCMC Assessment and Figures
#------------------------------------------------------------------------------#

library(network)
library(ergm)
library(multilayer.ergm)
library(texreg)
library(kableExtra)
library(stargazer)
library(coda)
library(plyr)
library(dplyr)
library(tidyr)
library(scales)


rm(list = ls(all.names = TRUE))


# setup
source("utils/utils.R")
source("utils/custom_trade_to_binary.R")
path <- data_path.get()
set.seed(1234)


# load
load(file = paste0(path, "/models/ERGM/ergm_results_2003_mcmle.RData"))



#------------------------------------------------------------------------------#
# MCMC Assessment
#------------------------------------------------------------------------------#

# output mcmc statistics
pdf(file = "figures/ergm_mcmc_1C_2003.pdf", width = 10, height = 0.9*sqrt(2)*10)
par(mfrow = c(6,2))
plot(fit_C$sample, smooth = TRUE, auto.layout = FALSE)
dev.off()

pdf(file = "figures/ergm_mcmc_1Cnn_2003.pdf", width = 10, height = 0.9*sqrt(2)*10)
par(mfrow = c(6,2))
plot(fit_C_nn$sample, smooth = TRUE, auto.layout = FALSE)
dev.off()

pdf(file = "figures/ergm_mcmc_1D_2003.pdf", width = 10, height = 0.9*sqrt(2)*10)
par(mfrow = c(6,2))
plot(fit_D$sample, smooth = TRUE, auto.layout = FALSE)
dev.off()

pdf(file = "figures/ergm_mcmc_1Dnn_2003.pdf", width = 10, height = 0.9*sqrt(2)*10)
par(mfrow = c(6,2))
plot(fit_D_nn$sample, smooth = TRUE, auto.layout = FALSE)
dev.off()


# Gweke Diagnostic?
mcmc.diagnostics(fit_C, which = c("burnin"))
mcmc.diagnostics(fit_C_nn, which = c("burnin"))
mcmc.diagnostics(fit_D, which = c("burnin"))
mcmc.diagnostics(fit_D_nn, which = c("burnin"))


# Effective Sample Size
neff_C <- data.frame("neff" = coda::effectiveSize(as.mcmc.list(fit_C$sample)), "Name" = names(coda::effectiveSize(as.mcmc.list(fit_C$sample))))
neff_C_nn <- data.frame("neff" = coda::effectiveSize(as.mcmc.list(fit_C_nn$sample)), "Name" = names(coda::effectiveSize(as.mcmc.list(fit_C_nn$sample))))
neff_D <- data.frame("neff" = coda::effectiveSize(as.mcmc.list(fit_D$sample)), "Name" = names(coda::effectiveSize(as.mcmc.list(fit_D$sample))))
neff_D_nn <- data.frame("neff" = coda::effectiveSize(as.mcmc.list(fit_D_nn$sample)), "Name" = names(coda::effectiveSize(as.mcmc.list(fit_D_nn$sample))))

neff_C %>% 
  rename("Import Dep. with" = neff) %>%
  full_join(neff_C_nn, by = "Name") %>%
  rename("Import Dep. without" = neff) %>%
  full_join(neff_D, by = "Name") %>%
  rename("Export Dep. with" = neff) %>%
  full_join(neff_D_nn, by = "Name") %>%
  rename("Export Dep. without" = neff)%>%
  relocate(Name) %>%
  pivot_longer(2:5, names_to = "Model", values_to = "eff") %>%
  group_by(Model) %>%
  summarise(Mean = mean(eff, na.rm = TRUE), SD = sd(eff, na.rm = TRUE)) -> neff




#------------------------------------------------------------------------------#
# Goodness of Fit Assessment
#------------------------------------------------------------------------------#

# recover observed cross layer statistics
obs_stats_C <- as.data.frame(as.list(summary(netC ~ duplexdyad(c("e", "f",  "h"), layers = list(1, 2)))))
obs_stats_D <- as.data.frame(as.list(summary(netD ~ duplexdyad(c("e", "f",  "h"), layers = list(1, 2)))))


# plot raster of distributions with observed and simulated cross-layer statistics
pdf(file = "figures/ergm_gof_duplex_2003.pdf", width = 10, height = 0.3*sqrt(2)*10)
par(mfrow = c(2,3))

plot(density(sim_cross_C[, "duplexdyad.e"]) , main = "Import Dependency Model: E", lty = 3)
lines(density(sim_cross_C_nn[, "duplexdyad.e"]), lty = 2)
abline(v = obs_stats_C$duplexdyad.e, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

plot(density(sim_cross_C[, "duplexdyad.f"]) , main = "Import Dependency Model: F", lty = 3)
lines(density(sim_cross_C_nn[, "duplexdyad.f"]), lty = 2)
abline(v = obs_stats_C$duplexdyad.f, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

plot(density(sim_cross_C[, "duplexdyad.h"]) , main = "Import Dependency Model: H", lty = 3)
lines(density(sim_cross_C_nn[, "duplexdyad.h"]), lty = 2)
abline(v = obs_stats_C$duplexdyad.h, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))


plot(density(sim_cross_D[, "duplexdyad.e"]) , main = "Export Dependency Model: E", lty = 3)
lines(density(sim_cross_D_nn[, "duplexdyad.e"]), lty = 2)
abline(v = obs_stats_D$duplexdyad.e, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

plot(density(sim_cross_D[, "duplexdyad.f"]) , main = "Export Dependency Model: F", lty = 3)
lines(density(sim_cross_D_nn[, "duplexdyad.f"]), lty = 2)
abline(v = obs_stats_D$duplexdyad.f, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

plot(density(sim_cross_D[, "duplexdyad.h"]) , main = "Export Dependency Model: H", lty = 3)
lines(density(sim_cross_D_nn[, "duplexdyad.h"]), lty = 2)
abline(v = obs_stats_D$duplexdyad.h, lty = 1)
legend("topright", legend = c("with", "without", "observed"), lty = c(3,2,1))

dev.off()


#------------------------------------------------------------------------------#
# Plot comparison of degree esp desp  
#------------------------------------------------------------------------------#

# C vs C_nn comparison
net_arms <- as.network(as.matrix.network(netC)[1:114, 1:114])
net_trade_C <- as.network(as.matrix.network(netC)[115:228, 115:228])
net_trade_D <- as.network(as.matrix.network(netD)[115:228, 115:228])

obs_stats_arms <- as.data.frame(as.list(summary(net_arms ~ edges + mutual + idegree(0:50) + odegree(0:50) + esp(0:50) + + dsp(0:50))))
obs_stats_C <- as.data.frame(as.list(summary(net_trade_C ~ edges + mutual + idegree(0:50) + odegree(0:50) + esp(0:50) + + dsp(0:50))))
obs_stats_D <- as.data.frame(as.list(summary(net_trade_D ~ edges + mutual + idegree(0:50) + odegree(0:50) + esp(0:50) + + dsp(0:50))))

sim_cross_C <- as.data.frame(sim_cross_C)
sim_cross_C_nn <- as.data.frame(sim_cross_C_nn)
sim_cross_D <- as.data.frame(sim_cross_D)
sim_cross_D_nn <- as.data.frame(sim_cross_D_nn)

sim_C <- as.data.frame(sim_C)
sim_C_nn <- as.data.frame(sim_C_nn)
sim_D <- as.data.frame(sim_D)
sim_D_nn <- as.data.frame(sim_D_nn)

geo_dist_arms <- ergm.geodistdist(net_arms)
names(geo_dist_arms) <- paste0("geodist.", names(geo_dist_arms))
geo_dist_C <- ergm.geodistdist(net_trade_C)
names(geo_dist_C) <- paste0("geodist.", names(geo_dist_C))
geo_dist_D <- ergm.geodistdist(net_trade_D)
names(geo_dist_D) <- paste0("geodist.", names(geo_dist_D))

range(sim_cross_C$edges_layer.1)
range(sim_cross_C$edges_layer.1)
range(sim_cross_C$edges_layer.1)
range(sim_cross_C$edges_layer.1)


## prepare data

# edges + reciprocity
bind_rows(list(sim_cross_C, sim_cross_C_nn)) %>%
  mutate(model = as.factor(c(rep("Import - with", 1000), rep("Import - without", 1000)))) -> pl_edges

# geodist distribution
bind_rows(list(sim_C, sim_C_nn)) %>%
  mutate(model = as.factor(c(rep("Import - with", 2000), rep("Import - without", 2000)))) -> pl_sim



#------------------------------------------------------------------------------#
# Figure Part 1
#------------------------------------------------------------------------------#

pdf(file = "figures/ergm_gof_import_models_pt1_2003.pdf",  width = 10, height = 0.9*sqrt(2)*10)
layout(matrix(c(1,2,3,4,
                5,5,5,6,
                7,7,7,8, 
                9,9,9,10,
                11,11,11,12), nrow = 5, ncol = 4, byrow = TRUE))

# Plot Edges
par(mar = c(5.1, 7.1, 4.1, 2.1))
boxplot(edges_layer.1 ~ model,  data = pl_edges, xlab = "Edges - MCW", ylab = "", horizontal=TRUE, las = 1, main = 1)
points(rep(obs_stats_arms$edges, 2), 1:2, pch = 19, col = "red", lwd = 4)

par(mar = c(5.1, 7.1, 4.1, 2.1))
boxplot(edges_layer.2 ~ model,  
        data = pl_edges, 
        xlab = "Edges - Conventional Trade", ylab = "", 
        horizontal=TRUE, 
        las = 1, main = 2)
points(c(obs_stats_C$edges, obs_stats_C$edges), 1:2, pch = 19, col = "red", lwd = 4)


# Plot Reciprocity
par(mar = c(5.1, 7.1, 4.1, 2.1))
boxplot(mutual.same.layer.mem.1 ~ model,  data = pl_edges, xlab = "Reciprocity - MCW", ylab = "", horizontal=TRUE, las = 1, main = 3)
points(rep(obs_stats_arms$mutual, 4), 4:1, pch = 19, col = "red", lwd = 4)

par(mar = c(5.1, 7.1, 4.1, 2.1))
boxplot(mutual.same.layer.mem.2 ~ model,  
        data = pl_edges, 
        xlab = "Reciprocity - Conventional Trade", ylab = "", 
        horizontal=TRUE, 
        las = 1, main = 4)
points(c(obs_stats_C$mutual, obs_stats_C$mutual), 1:2, pch = 19, col = "red", lwd = 4)


# Plot Geodist for ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("geodist.1", "geodist.2", "geodist.3", "geodist.4", "geodist.5", "geodist.6")
boxplot(tmp[, take], 
        xlab = "Geodesic Distance - MCW - with cross-layer effects", ylab = "",
        las = 1, main = 5)
points(1:length(take), geo_dist_arms[take], pch = 19, col = "red", lwd = 4)

boxplot(tmp[, c("geodist.Inf")], 
        xlab = "", ylab = "",
        las = 1, main = 6)
axis(side = 1, at = 1, label = "geodist.Inf")
points(1, geo_dist_arms[c("geodist.Inf")], pch = 19, col = "red", lwd = 4)


tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "Geodesic Distance - MCW - without cross-layer effects", ylab = "",
        las = 1, main = 7)
points(1:length(take), geo_dist_arms[take], pch = 19, col = "red", lwd = 4)

boxplot(tmp[, c("geodist.Inf")], 
        xlab = "", ylab = "",
        las = 1, main = 8)
axis(side = 1, at = 1, label = "geodist.Inf")
points(1, geo_dist_arms[c("geodist.Inf")], pch = 19, col = "red", lwd = 4)


# Plot Geodist for TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "Geodesic Distance - Conventional Trade - with cross-layer effects", ylab = "",
        las = 1, main = 9)
points(1:length(take), geo_dist_C[take], pch = 19, col = "red", lwd = 4)

boxplot(tmp[, c("geodist.Inf")], 
        xlab = "", ylab = "",
        las = 1, main = 10)
axis(side = 1, at = 1, label = "geodist.Inf")
points(1, geo_dist_C[c("geodist.Inf")], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "Geodesic Distance - Conventional Trade - without cross-layer effects", ylab = "",
        las = 1, main = 11)
points(1:length(take), geo_dist_C[take], pch = 19, col = "red", lwd = 4)

boxplot(tmp[, c("geodist.Inf")], 
        xlab = "", ylab = "",
        las = 1, main = 12)
axis(side = 1, at = 1, label = "geodist.Inf")
points(1, geo_dist_C[c("geodist.Inf")], pch = 19, col = "red", lwd = 4)


dev.off()



#------------------------------------------------------------------------------#
# Figure Part 2
#------------------------------------------------------------------------------#

## Plot Edgewise Shared Partners

pdf(file = "figures/ergm_gof_import_models_pt2_2003.pdf",  width = 10, height = 0.9*sqrt(2)*10)
layout(matrix(c(1,1,2,2,
                3,3,4,4,
                5,5,6,6,
                7,7,8,8), nrow = 4, ncol = 4, byrow = TRUE))
# Indegree ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("idegree0", "idegree1", "idegree2", "idegree3", "idegree4", "idegree5", "idegree6", "idegree7", "idegree8", "idegree9", "idegree10")
boxplot(tmp[, take], 
        xlab = "Indegree - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 1)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "Indegree - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 2)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# Indegree TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "Indegree - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 3)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "Indegree - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 4)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)



# Outdegree ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("odegree0", "odegree1", "odegree2", "odegree3", "odegree4", "odegree5", "odegree6", "odegree7", "odegree8", "odegree9", "odegree10")
boxplot(tmp[, take], 
        xlab = "Outdegree - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 5)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "Outdegree - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 6)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# Outdegree TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "Outdegree - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 7)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "Outdegree - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 8)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

dev.off()



#------------------------------------------------------------------------------#
# Figure Slides Part 1
#------------------------------------------------------------------------------#

pdf(file = "figures/ergm_gof_import_models_sl1_2003.pdf", width = 10, height = 6)
layout(matrix(c(1,1,2,2,
                3,3,4,4
), nrow = 2, ncol = 4, byrow = TRUE))

# Outdegree ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("odegree0", "odegree1", "odegree2", "odegree3", "odegree4", "odegree5", "odegree6", "odegree7", "odegree8", "odegree9", "odegree10")
boxplot(tmp[, take], 
        xlab = "Outdegree - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 1)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "Outdegree - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 2)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# Outdegree TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "Outdegree - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 3)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "Outdegree - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 4)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

dev.off()


#------------------------------------------------------------------------------#
# Figure Part 3
#------------------------------------------------------------------------------#

# Plot ESP
pdf(file = "figures/ergm_gof_import_models_pt3_2003.pdf",  width = 10, height = 0.9*sqrt(2)*10)
layout(matrix(c(1,1,2,2,
                3,3,4,4,
                5,5,6,6,
                7,7,8,8), nrow = 4, ncol = 4, byrow = TRUE))

# ESP ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("esp0", "esp1", "esp2", "esp3" ,"esp4" ,"esp5", "esp6" ,"esp7", "esp8", "esp9" ,"esp10" ,"esp11" ,"esp12" ,"esp13", "esp14", "esp15", "esp16", "esp17" ,"esp18" ,"esp19", "esp20")
boxplot(tmp[, take], 
        xlab = "ESP - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 1)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "ESP - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 2)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# ESP TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "ESP - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 3)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "ESP - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 4)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)


# DSP ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("dsp0", "dsp1", "dsp2", "dsp3" ,"dsp4" ,"dsp5", "dsp6" ,"dsp7", "dsp8", "dsp9" ,"dsp10" ,"dsp11" ,"dsp12" ,"dsp13", "dsp14", "dsp15", "dsp16", "dsp17" ,"dsp18" ,"dsp19", "dsp20")
boxplot(tmp[, take], 
        xlab = "DSP - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 5)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "DSP - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 6)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# DSP TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "DSP - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 7)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "DSP - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 8)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

dev.off()



#------------------------------------------------------------------------------#
# Figure Slides Part 2
#------------------------------------------------------------------------------#

pdf(file = "figures/ergm_gof_import_models_sl2_2003.pdf",  width = 10, height = 6)
layout(matrix(c(1,1,2,2,
                3,3,4,4), nrow = 2, ncol = 4, byrow = TRUE))

# ESP ARMS
tmp <- filter(pl_sim, model == "Import - with" & network == 1)
take <- c("esp0", "esp1", "esp2", "esp3" ,"esp4" ,"esp5", "esp6" ,"esp7", "esp8", "esp9" ,"esp10" ,"esp11" ,"esp12" ,"esp13", "esp14", "esp15", "esp16", "esp17" ,"esp18" ,"esp19", "esp20")
boxplot(tmp[, take], 
        xlab = "ESP - MCW - with cross-layer effects", 
        ylab = "", las = 1, main = 1)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 1)
boxplot(tmp[, take], 
        xlab = "ESP - MCW - without cross-layer effects", 
        ylab = "", las = 1, main = 2)
points(1:length(take), obs_stats_arms[take], pch = 19, col = "red", lwd = 4)


# ESP TRADE
tmp <- filter(pl_sim, model == "Import - with" & network == 2)
boxplot(tmp[, take], 
        xlab = "ESP - Conventional Trade - with cross-layer effects", 
        ylab = "", las = 1, main = 3)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

tmp <- filter(pl_sim, model == "Import - without" & network == 2)
boxplot(tmp[, take], 
        xlab = "ESP - Conventional Trade - without cross-layer effects", 
        ylab = "", las = 1, main = 4)
points(1:length(take), obs_stats_C[take], pch = 19, col = "red", lwd = 4)

dev.off()


#------------------------------------------------------------------------------#
# Plot Results
#------------------------------------------------------------------------------#

custom.coef.names <- c(
  "Edges",
  "Edges",
  "Reciprocity",
  "Reciprocity",
  "Gw Indegree (d = 1.5)",
  "Gw Indegree (d = 1.5)",
  "Gw Outdegree (d = 1.5)",
  "Gw Outdegree (d = 1.5)",
  "GWESP Outgoing Two-path (d = 0.69)",
  "GWESP Outgoing Two-path (d = 0.69)",
  "E",
  "F",
  "H",
  "Distance (log)",
  "Distance (log)",
  "GDP in (log)",
  "GDP in (log)",
  "GDP out (log)",
  "GDP out (log)",
  "Alliance",
  "Alliance",
  "Polity Diff. (abs)",
  "Polity Diff. (abs)",
  "Military Expenditure in (log)",
  "Military Expenditure in (log)",
  "Military Expenditure out (log)",
  "Military Expenditure out (log)",
  "Path Dependency",
  "Path Dependency", 
  "Path Dependency"
)

custom.coef.names2 <- c(
  "Edges",
  "Edges",
  "Reciprocity",
  "Reciprocity",
  "Gw Indegree (d = 1.5)",
  "Gw Indegree (d = 1.5)",
  "Gw Outdegree (d = 1.5)",
  "Gw Outdegree (d = 1.5)",
  "GWESP Outgoing Two-path (d = 0.69)",
  "GWESP Outgoing Two-path (d = 0.69)",
  "Distance (log)",
  "Distance (log)",
  "GDP in (log)",
  "GDP in (log)",
  "GDP out (log)",
  "GDP out (log)",
  "Alliance",
  "Alliance",
  "Polity Diff. (abs)",
  "Polity Diff. (abs)",
  "Military Expenditure in (log)",
  "Military Expenditure in (log)",
  "Military Expenditure out (log)",
  "Military Expenditure out (log)",
  "Path Dependency",
  "Path Dependency", 
  "Path Dependency"
)


sink("figures/ergm_effectivesize_2003.txt")
kbl(neff, booktabs = T,  format = "latex", 
    digits = 0,  escape = F, linesep = "",
    caption = "Effective Sample Size of the MCMC Sample Statistics", label = "ergm_effectivesize_2003") %>%
  kable_styling(font_size = 11, full_width = FALSE) %>%
  kable_styling(latex_options = "HOLD_position")
sink()


# output 
sink(file = "figures/ergm_estimates_1CD_2003.txt")
texreg(list(fit_C, fit_D),
       single.row = T,
       float.pos = "H",
       use.packages = FALSE,
       custom.coef.names = custom.coef.names, 
       custom.model.names = c("Import Dep.", "Export Dep."),
       booktabs = T, dcolumn = T, include.nobs = F,
       reorder.coef = c(1, 3, 5, 7, 9, 14, 16, 18, 20, 22, 24, 26, 28,
                        2, 4, 6, 8, 10, 15, 17, 19, 21, 23, 25, 27, 29, 
                        11, 12, 13), 
       groups = list("Layer 1: Arms Trade" = 1:13,
                     "Layer 2: Conventional Trade" = 14:26, 
                     "Cross Layer Network Effects" = 27:29),
       caption = "MERGM results for two-layer network of weapons and import (left) or export (right) trade dependency in the year 2003.", 
       label = "tab:ergm_estimates_model1C", 
       custom.note = "Estimates based on Monte Carlo MLE. Standard Errors in parenthesis.\\newline%stars. ")
sink()


# output 
sink(file = "figures/ergm_estimates_1CnnDnn_2003.txt")
texreg(list(fit_C_nn, fit_D_nn),
       single.row = T, 
       use.packages = FALSE,
       float.pos = "H",
       custom.coef.names = custom.coef.names2, 
       custom.model.names = c("Import Dep.", "Export Dep."),
       booktabs = T, dcolumn = T, include.nobs = F,
       reorder.coef = c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25,
                        2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26),
       groups = list("Layer 1: Arms Trade" = 1:13,
                     "Layer 2: Conventional Trade" = 14:26),
       caption = "MERGM results for two-layer network of weapons and import (left) or export (right) trade dependency in the year 2003 - estimated without cross-layer effects.", 
       label = "tab:ergm_estimates_model1CnnDnn", 
       custom.note = "Estimates based on Monte Carlo MLE. Standard Errors in parenthesis.\\newline%stars. ")
sink()





