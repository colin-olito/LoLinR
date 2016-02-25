BenchOut <- read.table("data-raw/BenchOut.txt", header=TRUE, stringsAsFactors=FALSE)

devtools::use_data(BenchOut, overwrite = TRUE)
