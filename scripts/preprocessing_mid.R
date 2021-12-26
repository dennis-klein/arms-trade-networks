library(data.table)
library(countrycode)


rm(list = ls(all.names = TRUE))

source("utils/utils.R")
path = data_path.get()


# load data
load(file.path(path, "out/country_list.RData"))
load(file.path(path, "out/EX.RData"))


# load mid data
midb = fread(file.path(path, "raw/Dyadic MIDs 4.02/dyadic_mid_4.02.csv"))


# helper function to melt data
help_func <- function(i) return(c(midb$strtyr[i]:midb$endyear[i]))
list_years <- lapply(1:nrow(midb), FUN = help_func)

midb[, year := list_years]
midb = midb[,.(year = unlist(year)), by = setdiff(names(midb), 'year')]
midb = midb[year %in% c(1950:2014)]

midb = unique(midb[, .(statea, stateb, year)])[, conflict := 1]

midb[, from_id:=match(statea, country_list$COW)]
midb[, to_id:=match(stateb, country_list$COW)]

unique(midb$statea[is.na(midb$from_id)]) 
unique(midb$stateb[is.na(midb$to_id)]) 

midb[statea == 260, from_id:= 72]
midb[stateb == 260, to_id:= 72]

midb = midb[!(is.na(midb$from_id) | is.na(midb$to_id))] # Delete all observations where the id is not known 
tmp = list()


# fill matrices
for (i in 1:length(1950:2014)){
  
  tmp[[i]] = matrix(0, 257, 257)
  colnames(tmp[[i]]) = country_list$V1
  rownames(tmp[[i]]) = country_list$V1
  
  yr = i + 1949
  tmp_dfl = midb[year == yr] 
  tmp[[i]][cbind(tmp_dfl$from_id, tmp_dfl$to_id)] <- tmp_dfl$conflict
}


# save
saveRDS(tmp, file = file.path(path, "out/midb_conflict.rds"))





