#------------------------------------------------------------------------------#
# A Network Plot for the Multilevel Network in 2003
#------------------------------------------------------------------------------#

# Corresponding Tutorials: 
# http://mr.schochastics.net/netVizR.html (for igraph / ggraph)
# http://blog.schochastics.net/post/visualizing-multilevel-networks-with-graphlayouts/

library(igraph)
library(ggraph)
library(graphlayouts)

rm(list = ls(all.names = TRUE))


# setup
source("utils/utils.R")
source("utils/construct_header.R")
source("utils/custom_trade_to_binary.R")
path <- data_path.get()
set.seed(1234)


# load data 
load(file.path(path, "out/EX.RData"))
load(file.path(path, "out/colony.RData"))
load(file.path(path, "out/country_list.RData"))
sipri_tiv <- readRDS(file.path(path, "out/sipri_tiv.rds"))
trade <- readRDS(file.path(path, "out/baci_aggregated.rds"))

thrshld1 <- 0
thrshld2 <- 100
year <- 2003
start <- 1995
end <- 2018

# selection of countries included: present over the complete time horizon (cf. soam)
included <- rowSums(EX[, (start:end) - 1949]) == length(start:end)
n <- sum(included)

mat1 <- (sipri_tiv[[year - 1949]][included, included] > thrshld1) * 1
mat2 <- (trade[[year - 1994]][included, included] > thrshld1) * 1

diag <- diag(1, n, n)

colnames(diag) <- country_list$V1[included]
rownames(diag) <- country_list$V1[included]

mat <- rbind(cbind(mat2, diag), cbind(diag, mat1))

net <- graph_from_adjacency_matrix(mat, mode = "directed", diag = FALSE)
net <- set_vertex_attr(net, "lvl", index = V(net), c(rep(1, n), rep(2, n)))
vertex_attr_names(net)


dev.new()
pdf(file = "figures/figure_descriptives_multilevel_2003.pdf", paper = "a4r", width = sqrt(2) * 10, height = 10 )
xy <- layout_as_multilevel(net, type = "fix1", FUN1 = layout_with_stress, alpha = 25, beta = 45)
ggraph(net, "manual", x = xy[, 1], y = xy[, 2]) +
  geom_edge_link0(
    aes(filter = (node1.lvl == 1 & node2.lvl == 1)),
    edge_colour = "dodgerblue3",
    alpha = 0.5,
    edge_width = 0.3
  ) +
  geom_edge_link0(
    aes(filter = (node1.lvl != node2.lvl)),
    alpha = 0.3,
    edge_width = 0.1,
    edge_colour = "black"
  ) +
  geom_edge_link0(
    aes(filter = (node1.lvl == 2 &
                    node2.lvl == 2)),
    edge_colour = "firebrick3",
    edge_width = 0.3,
    alpha = 0.5
  ) +
  geom_node_point(aes(shape = as.factor(lvl)), fill = "grey25", size = 3) +
  scale_shape_manual(values = c(21, 22)) +
  theme_graph() +
  coord_cartesian(clip = "off", expand = TRUE) +
  theme(legend.position = "none")

dev.off()




