taxes <- read.csv("data-raw/taxes.csv", header=TRUE, stringsAsFactors=FALSE)

devtools::use_data(taxes, overwrite = TRUE)
