BugulaData <- read.csv("data-raw/BugulaData.csv", header=TRUE, stringsAsFactors=FALSE)

devtools::use_data(BugulaData, overwrite = TRUE)
