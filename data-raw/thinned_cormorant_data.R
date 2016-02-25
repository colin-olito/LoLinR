thinned_cormorant_data <- read.csv("data-raw/thinned_cormorant_data.csv", header=TRUE, stringsAsFactors=FALSE, row.names=1)

devtools::use_data(thinned_cormorant_data, overwrite = TRUE)
