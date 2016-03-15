PelagicData <- read.csv("data-raw/PelagicData.csv", header=TRUE, stringsAsFactors=FALSE)

devtools::use_data(PelagicData, overwrite = TRUE)
