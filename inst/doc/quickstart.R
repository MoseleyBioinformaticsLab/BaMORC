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

## ------------------------------------------------------------------------
## Download a BMRB file
library(BMRBr)
bmrb_download(4020, output_dir = "./")

## Read in BMRB file and procec
file_path = "bmr4020.str"
head(read_NMRSTAR_file(file_path))

## Delete downloaded BMRB file
unlink("./bmr4020.str")

## ------------------------------------------------------------------------
id = 4022
output <- read_DB_File(id)
head(output[[1]])
head(output[[2]])

## ------------------------------------------------------------------------
protein_sequence <- "MQVWPIEGIKKFETLSYLPPLTVEDLLKQI"
secondary_structure <- jpred_fetcher(protein_sequence)

secondary_structure

## ------------------------------------------------------------------------


## ------------------------------------------------------------------------


## ------------------------------------------------------------------------
## Arguments:
sequence = paste(RefDB_data$carbonDat[[1]]$AA, collapse = "")
secondary_structure = paste(RefDB_data$carbonDat[[1]]$SS, collapse = "")

## Function:
calculate_RCF(sequence, secondary_structure)

## ------------------------------------------------------------------------
## chemicalShifts and aaFreq are predefined sample variables for demo purpose within the BaMORC Package.
calculate_MSE(step_ca=1, step_cb=1, dat_cacb=chemicalShifts[, c(3,4)], aa_Freq=aaFreq)

