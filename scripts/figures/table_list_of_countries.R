rm(list = ls(all.names = TRUE))

source("utils/utils.R")
path = data_path.get()


# load data
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/EX.RData"))


# select countries included in analysis
include = rowSums(EX[, (1995:2018)-1949]) == length(1995:2018)
country_list = country_list[include, ]
sum(include)


sink("figures/list_of_included_countries.txt")
cat(paste(country_list$V1, collapse = "\n \n"), sep = "\n")
sink()
