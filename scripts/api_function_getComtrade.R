library(rjson)
string <- "http://comtrade.un.org/data/cache/partnerAreas.json"
reporters <- fromJSON(file = string)
reporters <- as.data.frame(t(sapply(reporters$results, rbind)))


# https://comtrade.un.org/data/Doc/api/ex/r
get.Comtrade <- function(url = "http://comtrade.un.org/api/get?",
                         maxrec = 50000,
                         type = "C",
                         freq = "A",
                         px = "HS",
                         ps = "now",
                         r,
                         p,
                         rg = "all",
                         cc = "TOTAL",
                         fmt = "json") {
  string <- paste(url,
    "max=", maxrec, "&" # maximum no. of records returned
    , "type=", type, "&" # type of trade (c=commodities)
    , "freq=", freq, "&" # frequency
    , "px=", px, "&" # classification
    , "ps=", ps, "&" # time period
    , "r=", r, "&" # reporting area
    , "p=", p, "&" # partner country
    , "rg=", rg, "&" # trade flow
    , "cc=", cc, "&" # classification code
    , "fmt=", fmt # Format
    ,
    sep = ""
  )

  if (fmt == "csv") {
    raw.data <- read.csv(string, header = TRUE)
    return(list(validation = NULL, data = raw.data))
  } else {
    if (fmt == "json") {
      raw.data <- fromJSON(file = string)
      data <- raw.data$dataset
      validation <- unlist(raw.data$validation, recursive = TRUE)
      ndata <- NULL
      if (length(data) > 0) {
        var.names <- names(data[[1]])
        data <- as.data.frame(t(sapply(data, rbind)))
        ndata <- NULL
        for (i in 1:ncol(data)) {
          data[sapply(data[, i], is.null), i] <- NA
          ndata <- cbind(ndata, unlist(data[, i]))
        }
        ndata <- as.data.frame(ndata)
        colnames(ndata) <- var.names
      }
      return(list(validation = validation, data = ndata))
    }
  }
}




s <- get.Comtrade(r = "842", p = "124,484", fmt = "csv")
