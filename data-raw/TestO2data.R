TestO2data <- read.csv("data-raw/TestO2data.csv", header=TRUE, stringsAsFactors=FALSE)

devtools::use_data(TestO2data, overwrite = TRUE)
