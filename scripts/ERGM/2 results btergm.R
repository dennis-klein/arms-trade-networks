library(network)
library(xtable)
library(btergm)
library(multilayer.ergm)
library(ggplot2)
library(Cairo)
library(patchwork)
library(lemon)


rm(list = ls(all.names = TRUE))
source("utils/utils.R")
path = data_path.get()


# load model
version = "A"
load(file = file.path(path, paste0("models/ERGM/", version, "_estimates.RData")))


# gather estimates and confidence intervals for the full model 
tmp = list()

for (i in 1:length(results_full)){ 
  year = start + maxlag + i - 1
  tmp[[i]] <- data.frame("estimate" = results_full[[i]][, 1],
                         "lci_boot" = results_full[[i]][, 2],
                         "uci_boot" = results_full[[i]][, 3],
                         "se" = sqrt(diag(cov_full[[i]])),
                         "lci" = results_full[[i]][, 1] - qnorm(0.975) * sqrt(diag(cov_full[[i]])),
                         "uci" = results_full[[i]][, 1] + qnorm(0.975) * sqrt(diag(cov_full[[i]])),
                         "year" = year, 
                         "name" = gsub("[[i]]", "", rownames((results_full[[i]])), fixed = T))
}

confidence = do.call(rbind, tmp)


# some issues with se derived from simulations:
summary(confidence)


# plot estimates with base R: two different graphics for ex and end each
all = unique(confidence$name)
group1 = unique(confidence$name)[1:10]
group2 = unique(confidence$name)[11:26]
group3 = unique(confidence$name)[27:28]


CairoPDF(paste0("figures/2 temporal mergm covariates all ", version, ".pdf"),
  paper = "a4", width = 7, height = 10, pointsize = 9
)
par(mfrow = c(4, 2))
for (object in all) {
  tmp <- confidence[confidence$name == object, ]
  plot(tmp$year, tmp$estimate,
    type = "b", pch = 20,
    ylab = "Estimate", xlab = "Year",
    main = object,
    ylim = c(min(tmp$lci_boot), max(tmp$uci_boot)),
    cex.axis = 0.9
  )
  arrows(
    x0 = tmp$year, y0 = tmp$lci_boot,
    x1 = tmp$year, y1 = tmp$uci_boot,
    code = 3, angle = 90, length = 0.02
  )
  if ((min(tmp$lci_boot) <= 0) | (0 <= max(tmp$uci_boot))) {
    abline(h = 0, col = "grey")
  }
}
dev.off()


plots <- list()
i <- 0
for (object in group3) {
  i <- i + 1
  tmp <- confidence[confidence$name == object, ]
  plots[[i]] <- ggplot(tmp, aes(x = year, y = estimate)) +
    geom_pointline() +
    geom_errorbar(aes(ymin = lci_boot, ymax = uci_boot), width = 0.3) +
    labs(x = "Year", y = "Estimate", title = paste(object)) +
    geom_hline(yintercept = 0, color = "grey") +
    theme_classic()
}

ggsave(
  wrap_plots(plots) +
    plot_annotation(
      title = "Cross Layer Effects",
      caption = "A Reciprocity and B Reinforcement",
      tag_levels = "A"
    ),
  device = cairo_pdf, width = 12, height = 4, pointsize = 9, units = "in",
  filename = paste0("figures/2 temporal mergm covariates group3 ", version, ".pdf")
)



plots <- list()
i <- 0
for (object in group2) {
  i <- i + 1
  tmp <- confidence[confidence$name == object, ]
  plots[[i]] <- ggplot(tmp, aes(x = year, y = estimate)) +
    geom_pointline() +
    geom_errorbar(aes(ymin = lci_boot, ymax = uci_boot), width = 0.3) +
    labs(x = "Year", y = "Estimate", title = paste(object)) +
    geom_hline(yintercept = 0, color = "grey") +
    theme_classic()
}

ggsave(
  wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = "Fixed Effects",
      caption = "",
      tag_levels = "A"
    ),
  device = cairo_pdf, width = 12, height = 18, pointsize = 9, units = "in",
  filename = paste0("figures/2 temporal mergm covariates group2 ", version, ".pdf"),
)


plots <- list()
i <- 0
for (object in group1) {
  i <- i + 1
  tmp <- confidence[confidence$name == object, ]
  plots[[i]] <- ggplot(tmp, aes(x = year, y = estimate)) +
    geom_pointline() +
    geom_errorbar(aes(ymin = lci_boot, ymax = uci_boot), width = 0.3) +
    labs(x = "Year", y = "Estimate", title = paste(object)) +
    geom_hline(yintercept = 0, color = "grey") +
    theme_classic()
}

ggsave(
  wrap_plots(plots, ncol = 2) +
    plot_annotation(
      title = "Fixed Effects",
      caption = "",
      tag_levels = "A"
    ),
  device = cairo_pdf, width = 12, height = 18, pointsize = 9, units = "in",
  filename = paste0("figures/2 temporal mergm covariates group1 ", version, ".pdf")
)



# Compute Brier Score for two models and plot them
# we only care about the block diagonal networks and exclude diagonals
free = matrix(c(rep(c(rep(1, n), rep(0, n)), n), rep(c(rep(0, n), rep(1, n)), n)), 
              nrow = 2*n, ncol = 2*n, byrow = T)
diag(free) = 0 
score = matrix(NA, nrow = 4, ncol = length((start+maxlag+1):end), 
               dimnames = list(c("arms restricted", "arms full", 
                                 "trade restricted", "trade full"), 
                               c((start+maxlag+1):(end))))

for(i in seq_along(predictions_full)){
  
  future <- as.matrix.network(target[[i]])
  future[!free] <- 0
  
  #pred <- predictions_restricted[[i]]
  #pred[!free] <- 0
  #score[1, i] <- mean((as.sociomatrix(future)-pred)**2)
  
  pred <- predictions_full[[i]]
  pred[!free] <- 0
  score[2, i] <- mean((future[1:n, 1:n]-pred[1:n, 1:n])**2)
  score[4, i] <- mean((future[(n+1):(2*n), (n+1):(2*n)]- pred[(n+1):(2*n), (n+1):(2*n)])**2)
}


plot(score)



