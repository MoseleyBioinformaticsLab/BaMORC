% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{class}
\name{RefDB_data}
\alias{RefDB_data}
\title{RefDB object}
\format{A class contains 1557 13C datasets and access functions.}
\usage{
RefDB_data
}
\description{
Collection of data from Re-referenced Protein Chemical shift Database (RefDB) and their access functions.
}
\details{
\describe{
    \item{$AminoAcid}{data set of 19 amino acid single-letter naming convention (excluding glycine)}
    \item{$carbonDat}{data set (list) of all 1557 RefDB carbon 13 raw data pre-processed in data frame format}
    \item{$carbonDat_narm}{data set (list) of all 1557 RefDB carbon 13 raw data pre-processed in data frame format with NA removed}
    \item{$carbonDat_rmGU}{data set (list) of all 1557 RefDB carbon 13 raw data pre-processed in data frame format with NA, glycine and undetermined secondary structure elements removed}
    \item{$DBID}{data set of all 1557 RefDB IDs}
    \item{$getData(index = NA, ID = NA, type = "raw")}{function that return single dataset, (index or ID, only one is allowed)
       \itemize{
         \item index.          the index of the dataset
         \item ID.             the RefDB ID number of the dataset
         \item type            what kind of data will be fetched. "raw": returns raw data; "datatable": return data in data frame format; "rmNA": return data with NA removed; "rmGU: return data with NA, glycine and undetermined secondary structure elements removed
       }}
    \item{$getFreq(index = NA, ID = NA)}{function that return amino acid relative cumulative frequency for each data set (index or ID, only one is allowed)}
    \item{$getSecStr(index = NA, ID = NA)}{function that return secondary structure information for each data set (index or ID, only one is allowed)}
    \item{$getSequence(index = NA, ID = NA)}{function that return protein sequence information for each data set (index or ID, only one is allowed)}
    \item{$rawData}{the data set contains all the raw data from RefDB}

}
}
\keyword{datasets}
