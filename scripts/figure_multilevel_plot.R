# A Network Plot for the Multilevel Network in 2003

library(igraph)
library(ggraph)
library(graphlayouts)


set.seed(1234)
rm(list = ls(all.names = TRUE))


# load data 
load("data/out/EX.RData")
load("data/out/colony.RData")
load("data/out/country_list.RData")
load("data/out/sipri_tiv.RData")
baci_aggregated = readRDS("data/out/baci_aggregated.rds")

thrshld1 = 0
thrshld2 = 100 
year = 2001
present = (EX[, year - 1949] == 1)
n = sum(present)

mat1 = (sipri_tiv[[year - 1949]][present, present] > thrshld1) * 1
mat2 = (baci_aggregated[[year - 1994]][present, present] > thrshld2) * 1
diag = diag(1, n, n)
colnames(diag) = country_list$V1[present]
rownames(diag) = country_list$V1[present]
mat = rbind(cbind(mat2, diag), cbind(diag, mat1))
  
net = graph_from_adjacency_matrix(mat, mode = "directed", diag = FALSE)
net = set_vertex_attr(net, "lvl", index = V(net), c(rep(1, n), rep(2, n)))
vertex_attr_names(net)


pdf(file = "figures/plot_multilevel_2001.pdf", paper = "a4r")
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



# http://mr.schochastics.net/netVizR.html (for igraph / ggraph)
# http://blog.schochastics.net/post/visualizing-multilevel-networks-with-graphlayouts/
data("multilvl_ex", package = "graphlayouts")
vertex_attr_names(multilvl_ex)
vertex_attr(multilvl_ex, "lvl", index = V(multilvl_ex))

pdf(file = "figures/test.pdf", paper = "a4r")
cols2 <- c("#3A5FCD", "#CD00CD", "#EE30A7", "#EE6363", 
           "#CD2626", "#458B00", "#EEB422", "#EE7600")

xy <- layout_as_multilevel(multilvl_ex,type = "all", alpha = 25, beta = 45)
ggraph(multilvl_ex, "manual", x = xy[, 1], y = xy[, 2]) +
  geom_edge_link0(
    aes(filter = (node1.lvl == 1 & node2.lvl == 1)),
    edge_colour = "firebrick3",
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
    edge_colour = "goldenrod3",
    edge_width = 0.3,
    alpha = 0.5
  ) +
  geom_node_point(aes(shape = as.factor(lvl)), fill = "grey25", size = 3) +
  scale_shape_manual(values = c(21, 22)) +
  theme_graph() +
  coord_cartesian(clip = "off", expand = TRUE) +
  theme(legend.position = "none")

dev.off()