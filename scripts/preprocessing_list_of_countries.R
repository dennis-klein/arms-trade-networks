rm(list = ls(all.names = TRUE))

source("utils/utils.R")
path = data_path.get()


# load data
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/EX.RData"))


# select countries included in analysis
include = rowSums(EX[, (2000:2016)-1949]) == 17
country_list = country_list[include, ]
sum(include)


sink("figures/list of included countries.txt")
cat(paste(country_list$V1, collapse = "\n \n"), sep = "\n")
sink()
