## ---- echo = FALSE-------------------------------------------------------
library(httr)
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)

## ------------------------------------------------------------------------
library(BaMORC)

## ----cache=TRUE----------------------------------------------------------
## Arguments:
sequence = paste(RefDB_data$carbonDat[[1]]$AA, collapse = "")
secondary_structure = paste(RefDB_data$carbonDat[[1]]$SS, collapse = "")
chemical_shifts_input = RefDB_data$carbonDat[[1]][, c(4,5)]
from= -5
to = 5

## Running bamorc() function:
bamorc(sequence, secondary_structure, chemical_shifts_input, from=-5, to=5)

## ----cache=TRUE----------------------------------------------------------
## Arguments:
sequence = "RPAFCLEPPYAGPGKARIIRYFYNAAAGAAQAFVYGGVRAKRNNFASAADALAACAAA"
sample_data_generator(input_type = "ssc_sample") # this will generate a temperary sample file and later will be removed.
file_path = "./bpti_HNcoCACB.txt" # temperary sample file path.

## Running unassigned_bamorc() function:
unassigned_bamorc(peakList_file_loc = file_path, sequence = sequence, secondary_structure = NULL, from = -5, to = 5)

## Delete the temperary sample file.
unlink("./bpti_HNcoCACB.txt")

## ------------------------------------------------------------------------
input_type = "ws" 
sample_data_generator(input_type = input_type)
head(read_raw_file("sample_input_ws.txt", delim = "ws"))
unlink("sample_input_ws.txt")

## ------------------------------------------------------------------------
input_type = "csv" 
sample_data_generator(input_type = input_type)
head(read_raw_file("sample_input.csv", delim = "comma"))
unlink("sample_input.csv")

## ------------------------------------------------------------------------
input_type = "sc" 
sample_data_generator(input_type = input_type)
head(read_raw_file("sample_input_sc.txt", delim = "semicolon"))
unlink("sample_input_sc.txt")

