CormorantData <- read.csv("data-raw/CormorantData.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

devtools::use_data(CormorantData, overwrite = TRUE)
