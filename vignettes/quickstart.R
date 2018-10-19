## ---- echo = FALSE-------------------------------------------------------
library(knitr)
opts_chunk$set(tidy.opts=list(width.cutoff=60),tidy=TRUE)
knitr::opts_chunk$set(comment = "#>", collapse = TRUE)

## ------------------------------------------------------------------------
library(BaMORC)

## ------------------------------------------------------------------------
## Arguments:
sequence = paste(RefDB_data$carbonDat[[1]]$AA, collapse = "")
secondary_structure = paste(RefDB_data$carbonDat[[1]]$SS, collapse = "")
chemical_shifts_input = RefDB_data$carbonDat[[1]][, c(4,5)]
from= -5
to = 5

## Running bamorc() function:
bamorc(sequence = sequence, secondary_structure = secondary_structure, 
       chemical_shifts_input = chemical_shifts_input, from = -5, to = 5)

## ---- eval= FALSE--------------------------------------------------------
#  ## Arguments:enerate a temperary sample NMR spectra file and later will be
#  ## removed.
#  sequence = "RPAFCLEPPYAGPGKARIIRYFYNAAAGAAQAFVYGGVRAKRNNFASAADALAACAAA"
#  sample_data_generator(input_type = "ssc_sample")
#  file_path = "./bpti_HNcoCACB.txt" # temperary sample file path.
#  
#  ## Running unassigned_bamorc() function:
#  # unassigned_bamorc will take a while to run
#  unassigned_bamorc(peakList_file_loc = file_path, sequence = sequence,
#                   secondary_structure = NULL, from = -5, to = 5)
#  
#  ## Delete the temperary sample file.
#  unlink("./bpti_HNcoCACB.txt")

## ------------------------------------------------------------------------
## Arguments:enerate a temperary sample NMR spectra file and later will be 
## removed.
input_type = "ws" 
sample_data_generator(input_type = input_type)

## Running reading function
head(read_raw_file(file_path = "sample_input_ws.txt", delim = "ws"))
unlink("sample_input_ws.txt")

## ------------------------------------------------------------------------
## Arguments:enerate a temperary sample NMR spectra file and later will be 
## removed.
input_type = "csv" 
sample_data_generator(input_type = input_type)

## Running reading function
head(read_raw_file(file_path = "sample_input.csv", delim = "comma"))
unlink("sample_input.csv")

## ------------------------------------------------------------------------
## Arguments:enerate a temperary sample NMR spectra file and later will be 
## removed.
input_type = "sc" 
sample_data_generator(input_type = input_type)

## Running reading function
head(read_raw_file(file_path = "sample_input_sc.txt", delim = "semicolon"))
unlink("sample_input_sc.txt")

## ---- eval= FALSE--------------------------------------------------------
#  ## Download a BMRB file
#  library(BMRBr)
#  bmrb_download(4020, output_dir = "./")
#  
#  ## Read in BMRB file and procec
#  file_path = "bmr4020.str"
#  head(read_nmrstar_file(file_path = file_path))
#  
#  ## Delete downloaded BMRB file
#  unlink("./bmr4020.str")

## ------------------------------------------------------------------------
id = 4022
output <- read_db_file(id = id)
head(output[[1]])
head(output[[2]])

## ---- eval= FALSE--------------------------------------------------------
#  protein_sequence <- "MQVWPIEGIKKFETLSYLPPLTVEDLLKQI"
#  secondary_structure <- jpred_fetcher(protein_sequence = protein_sequence)
#  
#  secondary_structure

## ------------------------------------------------------------------------
calculate_AA_Prob(chi_squared_stat = c(0.05, 0.1, 0.5), df = 2)

## ------------------------------------------------------------------------
calculate_chi_squared_stat(cacb_pair = c(54,45))

## ------------------------------------------------------------------------
## Arguments:
sequence = paste(RefDB_data$carbonDat[[1]]$AA, collapse = "")
secondary_structure = paste(RefDB_data$carbonDat[[1]]$SS, collapse = "")

## Function:
calculate_RCF(sequence = sequence, secondary_structure = secondary_structure)

## ------------------------------------------------------------------------
## chemicalShifts and aaFreq are predefined sample variables for demo purpose 
## within the BaMORC Package.
calculate_MSE(step_ca = 1, step_cb = 1, dat_cacb = chemicalShifts[, c(3,4)], 
              aa_Freq = aaFreq)

