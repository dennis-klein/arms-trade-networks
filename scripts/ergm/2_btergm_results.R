library(network)
library(xtable)
library(btergm)
library(multilayer.ergm)
library(ggplot2)
library(Cairo)
library(patchwork)
library(lemon)
library(scales)


rm(list = ls(all.names = TRUE))
source("utils/utils.R")
path = data_path.get()


# load model
version = "A"
load(file = file.path(path, paste0("models/ERGM/", version, "_estimates.RData")))
load(file = paste0(path, "/models/ERGM/", version, "_predictive_scores.RData"))


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



# first plot with all network statistics
plots <- list()
i <- 0
for (object in group1) {
  i <- i + 1
  tmp <- confidence[confidence$name == object, ]
  plots[[i]] <- ggplot(tmp, aes(x = year, y = estimate, ymin = lci_boot, ymax = uci_boot)) +
    geom_pointrange() +
    coord_flip() +
    labs(x = "Year", y = "Estimate", title = paste(object)) +
    theme(text = element_text(size = 14))
  
  if (!(object %in% c("edges_layer.2", "gwideg1.layer.mem1", "gwideg1.layer.mem2", "gwodeg1.layer.mem1"))) {
    plots[[i]] <- plots[[i]] +
      geom_hline(yintercept = 0, lty = 2)
  }
}
ggsave(
  wrap_plots(plots, nrow = 2, byrow = F) +
    plot_annotation(
      title = "",
      caption = "",
      tag_levels = "A"
    ),
  device = cairo_pdf, width = 18, height = 12, units = "in",
  filename = paste0("figures/2 temporal mergm covariates network ", version, ".pdf")
)


# second plot with part 1 and part 2
plots <- list()
i <- 0
for (object in group2) {
  i <- i + 1
  tmp <- confidence[confidence$name == object, ]
  plots[[i]] <- ggplot(tmp, aes(x = year, y = estimate, ymin = lci_boot, ymax = uci_boot)) +
    geom_pointrange() +
    coord_flip() +
    labs(x = "Year", y = "Estimate", title = paste(object)) +
    theme(text = element_text(size = 14))
  
  if (!(object %in% c("edgecov.pathdep_layer1", "edgecov.pathdep_layer2"))){
    plots[[i]] <- plots[[i]] + 
      geom_hline(yintercept = 0, lty = 2)
  }
}

ggsave(
  wrap_plots(plots[1:10], nrow = 2, byrow = F) +
    plot_annotation(
      title = "",
      caption = "",
      tag_levels = "A"
    ),
  device = cairo_pdf, width = 18, height = 12, units = "in",
  filename = paste0("figures/2 temporal mergm covariates part1 ", version, ".pdf"),
)


# plots of missing others and endogenous effects
for (object in group3) {
  i <- i + 1
  tmp <- confidence[confidence$name == object, ]
  plots[[i]] <- ggplot(tmp, aes(x = year, y = estimate, ymin = lci_boot, ymax = uci_boot)) +
    geom_pointrange() +
    coord_flip() +
    labs(x = "Year", y = "Estimate", title = paste(object)) +
    geom_hline(yintercept = 0, lty = 2) +
    theme(text = element_text(size = 14))
}

layout <- "
ACE#G
BDF#H
"
  
ggsave(
  wrap_plots(plots[11:18], byrow = F, design = layout) +
    plot_annotation(
      title = "",
      caption = "",
      tag_levels = "A"
    ),
  device = cairo_pdf, width = 18, height = 12, pointsize = 9, units = "in",
  filename = paste0("figures/2 temporal mergm covariates part2 ", version, ".pdf")
)



tmp <- subset(scores, metric == "brier" & type == "arms")

# plot scores 
plot_scores <- list()
plot_scores[[1]] <- ggplot(
  subset(scores, metric == "brier" & type == "arms"),
  aes(x = year, y = value, colour = model, group = interaction(model, type))
) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, color = "grey") +
  theme_classic() +
  scale_colour_viridis_d("Model") +
  labs(subtitle = "Brier Score: only Arms", x = "Year", y = "") +
  scale_y_continuous(expand = c(0, 0), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_continuous(expand = c(0, 0))


plot_scores[[2]] <- ggplot(
  subset(scores, metric == "brier" & type == "total"),
  aes(x = year, y = value, colour = model, group = interaction(model, type))
) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, color = "grey") +
  theme_classic() +
  scale_colour_viridis_d("Model") +
  labs(subtitle = "Brier Score: total", x = "Year", y = "") +
  scale_y_continuous(expand = c(0, 0), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_continuous(expand = c(0, 0))


plot_scores[[3]] <- ggplot(
  subset(scores, metric == "pr" & type == "arms"),
  aes(x = year, y = value, colour = model, group = interaction(model, type))
) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, color = "grey") +
  theme_classic() +
  scale_colour_viridis_d("Model") +
  labs(subtitle = "Precision Recall AUC: only Arms", x = "Year", y = "") +
  scale_y_continuous(expand = c(0, 0), labels = scales::number_format(accuracy = 0.01)) +
  scale_x_continuous(expand = c(0, 0))


plot_scores[[4]] <- ggplot(
  subset(scores, metric == "pr" & type == "total"),
  aes(x = year, y = value, colour = model, group = interaction(model, type))
) +
  geom_point() +
  geom_line() +
  theme_classic() +
  scale_colour_viridis_d("Model") +
  labs(subtitle = "Precision Recall AUC: total", x = "Year", y = "") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01))


ggsave(
  wrap_plots(plot_scores, ncol = 2) +
    plot_annotation(
      title = "Comparison in Predictive Capacity: Brier Score and Precision Recall AUC",
      caption = "",
      tag_levels = "A"
    ) +plot_layout(guides = "collect") & theme(legend.position = "bottom"),
  device = cairo_pdf, width = 12, height = 12, pointsize = 14, units = "in",
  filename = paste0("figures/2 predictive comparison ", version, ".pdf")
)


scores_sum <- aggregate(value ~ model + metric + type, data = scores, sum)
scores_sum <- arrange(scores_sum, type, metric)

xtable(scores_sum[1:4,])
