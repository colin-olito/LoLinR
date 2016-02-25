thinned_toad_data <- read.csv("data-raw/thinned_toad_data.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

devtools::use_data(thinned_toad_data, overwrite = TRUE)
