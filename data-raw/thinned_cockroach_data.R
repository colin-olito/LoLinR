thinned_cockroach_data <- read.csv("data-raw/thinned_cockroach_data.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

devtools::use_data(thinned_cockroach_data, overwrite = TRUE)
