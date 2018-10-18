#' RefDB object
#'
#' Collection of data from Re-referenced Protein Chemical shift Database (RefDB) and their access functions.
#'

#'
#' \describe{
#'     \item{$AminoAcid}{data set of 19 amino acid single-letter naming convention (excluding glycine)}
#'     \item{$carbonDat}{data set (list) of all 1557 RefDB carbon 13 raw data pre-processed in data frame format}
#'     \item{$carbonDat_narm}{data set (list) of all 1557 RefDB carbon 13 raw data pre-processed in data frame format with NA removed}
#'     \item{$carbonDat_rmGU}{data set (list) of all 1557 RefDB carbon 13 raw data pre-processed in data frame format with NA, glycine and undetermined secondary structure elements removed}
#'     \item{$DBID}{data set of all 1557 RefDB IDs}
#'     \item{$getData(index = NA, ID = NA, type = "raw")}{function that return single dataset, (index or ID, only one is allowed)
#'        \itemize{
#'          \item index.          the index of the dataset
#'          \item ID.             the RefDB ID number of the dataset
#'          \item type            what kind of data will be fetched. "raw": returns raw data; "datatable": return data in data frame format; "rmNA": return data with NA removed; "rmGU: return data with NA, glycine and undetermined secondary structure elements removed
#'        }}
#'     \item{$getFreq(index = NA, ID = NA)}{function that return amino acid relative cumulative frequency for each data set (index or ID, only one is allowed)}
#'     \item{$getSecStr(index = NA, ID = NA)}{function that return secondary structure information for each data set (index or ID, only one is allowed)}
#'     \item{$getSequence(index = NA, ID = NA)}{function that return protein sequence information for each data set (index or ID, only one is allowed)}
#'     \item{$rawData}{the data set contains all the raw data from RefDB}
#'
#' }
#'
#' @docType class
#' @name RefDB_data
#' @usage RefDB_data
#' @format A class contains 1557 13C datasets and access functions.
"RefDB_data"


#' Overlapping Matrix for 57 amino acid and secondary structure combinations.
#'
#' A dataset containing the pre-calculated overlapping matrix 57 amino acid and secondary structure combinations:
#'
#' \itemize{
#'   \item column names or row names.
#'     \itemize{
#'       \item First letter. amino acid typings in single-letter format, not include glycine ("A", "C", "D", "E", "F", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
#'       \item second letter. secondary structure in single-letter format ("B", "H", "C")
#'     }
#' }
#'
#' @docType data
#' @keywords datasets
#' @name AA57OLMatrix
#' @usage AA57OLMatrix
#' @format A matrix with 57 by 57 elements.
'AA57OLMatrix'


#' Overlapping weights for 57 amino acid and secondary structure combinations.
#'
#' A vector containing the pre-calculated overlapping weights 57 amino acid and secondary structure combinations:
#'
#' \itemize{
#'   \item names.
#'     \itemize{
#'       \item First letter. amino acid typings in single-letter format, not include glycine ("A", "C", "D", "E", "F", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
#'       \item second letter. secondary structure in single-letter format ("B", "H", "C")
#'     }
#' }
#'
#' @docType data
#' @keywords datasets
#' @name AA57OLWeights
#' @usage AA57OLWeights
#' @format A vector of 57 elements.
'AA57OLWeights'

#' Three-letter amino acid naming convention with first letter capitalized.
#'
#' A vector containing amino acid names with first letter capitalized (no glycine), they are:
#' "Ala", "Cys", "Asp", "Glu", "Phe", "His", "Ile", "Lys", "Leu", "Met", "Asn", "Pro", "Gln", "Arg", "Ser", "Thr", "Val", "Trp", "Tyr".
#'
#' @docType data
#' @keywords datasets
#' @name aaCodes3Letter1stCap
#' @usage aaCodes3Letter1stCap
#' @format A vector of 19 elements.
'aaCodes3Letter1stCap'

#' Three-letter amino acid naming convention with all letters capitalized.
#'
#' A vector containing amino acid names with all letters capitalized (no glycine), they are:
#' "ALA", "CYS", "ASP", "GLU", "PHE", "HIS", "ILE", "LYS", "LEU", "MET", "ASN", "PRO", "GLN", "ARG", "SER", "THR", "VAL", "TRP", "TYR".
#'
#' @docType data
#' @keywords datasets
#' @name aaCodes3LetterAllCap.v
#' @usage aaCodes3LetterAllCap.v
#' @format A vector of 19 elements.
'aaCodes3LetterAllCap'

#' Single-letter amino acid naming convention.
#'
#' A vector containing amino acid names with all letters capitalized, (no glycine, but includes the oxidized cystine), they are:
#' "A", "B", "C", "D", "E", "F", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y".
#'
#' @docType data
#' @keywords datasets
#' @name aaCodes1Letter20
#' @usage aaCodes1Letter20
#' @format A vector of 20 elements.
'aaCodes1Letter20'

#' Single-letter amino acid naming convention.
#'
#' A vector containing amino acid names with all letters capitalized, (no glycine), they are:
#' "A", "C", "D", "E", "F", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y".
#'
#' @docType data
#' @keywords datasets
#' @name aaCodes1Letter19
#' @usage aaCodes1Letter19
#' @format A vector of 19 elements.
'aaCodes1Letter19'

#' Average covariance across three secondary structures for all amino acid typings (including oxidized cystine).
#'
#' A dataset containing average covariance across three secondary structure for all amino acid typings (including oxidized cystine):
#' \itemize{
#'   \item AA. amino acid typing names
#'   \item AvgCov. average covariance value across three secondary structures
#' }
#'
#' @docType data
#' @keywords datasets
#' @name AvgCov.t
#' @usage AvgCov.t
#' @format A data frame with 20 rows and 2 columns.
'AvgCov.t'

#' Mean chemical shift values for alpha carbon .
#'
#' A dataset containing mean chemical shift values for alpha carbon for each amino acid across three secondary structure for all amino acid typings (including oxidized cystine and glycine):
#' \itemize{
#'   \item Residue. amino acid typing name (single-letter convention)
#'   \item C. coil
#'   \item H. helix
#'   \item B. beta strand
#'   \item A. average of three secondary structures.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name CAMuTable
#' @usage CAMuTable
#' @format A data frame with 21 rows and 5 columns.
'CAMuTable'

#' Standard deviation of chemical shift values for alpha carbon .
#'
#' A dataset containing Standard deviation of chemical shift values for alpha carbon for each amino acid across three secondary structure for all amino acid typings (including oxidized cystine and glycine):
#' \itemize{
#'   \item Residue. amino acid typing name (single-letter convention)
#'   \item C. coil
#'   \item H. helix
#'   \item B. beta strand
#'   \item A. average of three secondary structures.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name CASdTable
#' @usage CASdTable
#' @format A data frame with 21 rows and 5 columns.
'CASdTable'

#' Mean chemical shift values for beta carbon .
#'
#' A dataset containing mean chemical shift values for beta carbon for each amino acid across three secondary structure for all amino acid typings (including oxidized cystine and glycine):
#' \itemize{
#'   \item Residue. amino acid typing name (single-letter convention)
#'   \item C. coil
#'   \item H. helix
#'   \item B. beta strand
#'   \item A. average of three secondary structures.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name CBMuTable
#' @usage CBMuTable
#' @format A data frame with 21 rows and 5 columns.
'CBMuTable'

#' Standard deviation of chemical shift values for beta carbon .
#'
#' A dataset containing Standard deviation of chemical shift values for beta carbon for each amino acid across three secondary structure for all amino acid typings (including oxidized cystine and glycine):
#' \itemize{
#'   \item Residue. amino acid typing (single-letter convention)
#'   \item C. coil
#'   \item H. helix
#'   \item B. beta strand
#'   \item A. average of three secondary structures.
#' }
#'
#' @docType data
#' @keywords datasets
#' @name CBSdTable
#' @usage CBSdTable
#' @format A data frame with 21 rows and 5 columns.
'CBSdTable'

#' Covariance values of chemical shifts of alpha and beta carbons.
#'
#' A dataset containing covariance values of chemical shifts of alpha and beta carbons for each amino acid across three secondary structure for all amino acid typings (including oxidized cystine but not glycine):
#' \itemize{
#'   \item AA.SS. amino acid typing and secondary structure combination (single-letter convention)
#'   \item Value. covariance value

#' }
#'
#' @docType data
#' @keywords datasets
#' @name CarbonCov.t
#' @usage CarbonCov.t
#' @format A data frame with 79 rows and 2 columns.
'CarbonCov.t'

#' RefDB ID included in the BaMORC package.
#'
#' A vector containing all the RefDB ID included in the BaMORC package.
#'
#' @docType data
#' @keywords datasets
#' @name ID
#' @usage ID
#' @format A vector of 1557 elements.
'ID'

#' Statistics of chemical shifts values of alpha carbons from RefDB.
#'
#' A dataset containing statistics of chemical shifts of alpha carbons for each amino acid from RefDB data:
#' \itemize{
#'   \item Residue. amino acid typing (single-letter convention)
#'   \item Coil.mu. mean value of alpha carbon in coil
#'   \item Coil.sd. mean value of alpha carbon in coil
#'   \item Helix.mu. mean value of alpha carbon in helix
#'   \item Helix.sd. standard deviation value of alpha carbon in helix
#'   \item Beta.mu. mean value of alpha carbon in coil
#'   \item Beta.sd. standard deviation value of alpha carbon in coil
#'   \item Avg.mu. average of mean values of alpha carbon across three secondary structures
#'   \item Avg.sd. average of standard deviation values of alpha carbon across three secondary structures
#' }
#'
#' @docType data
#' @keywords datasets
#' @name RefDB.StatCA
#' @usage RefDB.StatCA
#' @format A data frame with 21 rows and 9 columns.
'RefDB.StatCA'

#' Statistics of chemical shifts values of beta carbons from RefDB.
#'
#' A dataset containing statistics of chemical shifts of beta carbons for each amino acid from RefDB data:
#' \itemize{
#'   \item Residue. amino acid typing (single-letter convention)
#'   \item Coil.mu. mean value of beta carbon in coil
#'   \item Coil.sd. mean value of beta carbon in coil
#'   \item Helix.mu. mean value of beta carbon in helix
#'   \item Helix.sd. standard deviation value of beta carbon in helix
#'   \item Beta.mu. mean value of beta carbon in coil
#'   \item Beta.sd. standard deviation value of beta carbon in coil
#'   \item Avg.mu. average of mean values of beta carbon across three secondary structures
#'   \item Avg.sd. average of standard deviation values of beta carbon across three secondary structures
#' }
#'
#' @docType data
#' @keywords datasets
#' @name RefDB.StatCB
#' @usage RefDB.StatCB
#' @format A data frame with 21 rows and 9 columns.
'RefDB.StatCB'

#' All the amino acids and secondary structures combinations for easy access.
#'
#' A dataset containing all the amino acids and secondary structures combinations:
#' A-H", "A-B", "A-C", "R-H", "R-B", "R-C", "N-H", "N-B", "N-C", "D-H", "D-B", "D-C", "B-H", "B-B", "B-C", "C-H", "C-B", "C-C", "Q-H", "Q-B", "Q-C", "E-H", "E-B", "E-C", "H-H", "H-B", "H-C", "I-H", "I-B", "I-C", "L-H", "L-B", "L-C", "K-H", "K-B", "K-C", "M-H", "M-B", "M-C", "F-H", "F-B", "F-C", "P-H", "P-B", "P-C", "S-H", "S-B", "S-C", "T-H", "T-B", "T-C", "Y-H", "Y-B", "Y-C", "W-H", "W-B", "W-C", "V-H", "V-B", "V-C".
#' @docType data
#' @keywords datasets
#' @name cname
#' @usage cname
#' @format A vector with 60 elements.
'cname'

#' inverseMatrices.
#'
#' A data set containing all the pre-calculated inverse matrices and the access functions.
#' \describe{
#'  \item{$getInvMatrix(name)}{get the inverse matrix basing on the name, which is the amino acid and secondary structure combination (single-letter naming convention)}
#'  \item{$getNames()}{get all the names of the inverse matrices}
<<<<<<< HEAD
#'  \item{$Matrices}{the data set of all the pre-calculated inverse matrices}
#' }
#' @docType class
#' @name inverseMatrices
#' @format An object of inverse matrices with their accessing functions.
'inverseMatrices'
=======
#'  \item{$Matrices}{the data set of all the pre-calculated inversed matricies}
#' }
#' @docType class
#' @name inversedMatrices
#' @format An object of inversed matrices with their accessing functions.
'inversedMatrices'
>>>>>>> 3f51f9f9c7f032a01412e344093c3dcc19241810


#' Pre-defined sample chemical shifts data.
#'
#' Sample data for testing \code{calculate_MSE()} function.
#' \itemize{
#'   \item AA. amino acid typing (single-letter convention)
#'   \item SS. secondary structure (single-letter convention)
#'   \item CA. chemical shift value of alpha carbon
#'   \item CB. chemical shift value of beta carbon
#' }
#'
#' @docType data
#' @keywords datasets
#' @name chemicalShifts
#' @usage chemicalShifts
#' @format A data frame with 55 rows and 4 columns.
'chemicalShifts'

#' Pre-defined amino acid frequency data.
#'
#' Sample data for testing \code{calculate_MSE()} function.
#' \itemize{
#'   \item AA_SS. amino acid typing and secondary structure naming combination (single-letter convention)
#'   \item Freq. relative cumulated frequency
#' }
#'
#' @docType data
#' @keywords datasets
#' @name aaFreq
#' @usage aaFreq
#' @format A data frame with 57 rows and 2 columns.
'aaFreq'



