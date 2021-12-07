# utils

# get local path to data repository
data_path.get <- function() {
  r <- readLines("data/data_path.txt", warn = F)
  return(r)
}

