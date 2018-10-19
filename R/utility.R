#' Extracts data from a protein NMR experimental peak list.
#' \code{read_raw_file()} function reads in a user provided protein NMR experimental peak list. It currently supports file format in csv, txt with deliminator of comma, whitespace or semicolon.
#' Note: please don't leave space between sequence and chemical shifts data, otherwise it will report error.
#' @param file_path     File path where input chemical shifts file is located
#' @param delim         Delimiter for parsing file
#' @param assigned      Flag tell whether the input chemical shifts file is already assigned or not
#'
#' @return              A list contains protein sequence and chemical shift table.
#' @export              read_raw_file
#'
#' @examples
#' input_type = "ws" 
#' sample_data_generator(input_type = input_type)
#' head(read_raw_file("sample_input_ws.txt", delim="ws"))
#' unlink("sample_input_ws.txt")
#'
#' input_type = "csv"
#' sample_data_generator(input_type = input_type)
#' head(read_raw_file("sample_input.csv", delim="comma"))
#' unlink("sample_input.csv")
#'
#' input_type = "sc"
#' sample_data_generator(input_type = input_type)
#' head(read_raw_file("sample_input_sc.txt", delim="semicolon"))
#' unlink("sample_input_sc.txt")
#'

read_raw_file <- function(file_path, delim="comma", assigned=FALSE){
        skip_line = ifelse(assigned, 0, 1)

        raw_doc <- switch(delim,
                          comma     = readr::read_csv(file_path, col_names = FALSE, skip = skip_line),
                          ws        = readr::read_table2(file_path, col_names = FALSE, skip = skip_line),
                          semicolon = readr::read_delim(file_path, ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = skip_line),
                          tab       = readr::read_delim(file_path, "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = skip_line)
        )
        if(assigned==FALSE){
                con <- file(file_path, "r")
                seq <- readLines(con, n = 1)
                close(con)
                return(list(seq, raw_doc))
        }
        else{
                return(raw_doc)
        }



}



#' Extracts data from BMRB STAR 3.0 file.
#' \code{read_nmrstar_file()} parses BMRB STAR 3.0 file. It will extract sequence information and chemical shifts for both alpha and beta carbons.
#'
#' @param file_path     File path where input chemical shifts file is located
#' @importFrom BMRBr bmrb_download
#' @importFrom RBMRB export_star_data
#' @importFrom RBMRB fetch_entry_chemical_shifts
#' @importFrom tidyr spread
#' @importFrom stringr str_replace_all
#' @return Protein sequence and chemical shifts dataframe.
#' @export read_nmrstar_file
#'
#' @examples
#' ## Download a BMRB file
#' library(BMRBr)
#' \dontrun{bmrb_download(id_list = "4020", output_dir = "./", verbose = F)}
#'
#' ## Read in BMRB file and procec
#' file_path = "bmr4020.str"
#' \dontrun{head(read_nmrstar_file(file_path))}
#'
#' ## Delete downloaded BMRB file
#' unlink("./bmr4020.str")
#'
read_nmrstar_file <- function(file_path){

        token <- RBMRB::export_star_data(file_path)
        raw_data <- RBMRB::fetch_entry_chemical_shifts(token)[,c("Seq_ID", "Comp_ID", "Atom_ID", "Val")]
        raw_data_carbon <- raw_data[which(raw_data$Atom_ID=="CA" | raw_data$Atom_ID=="CB"),]
        dat <- tidyr::spread(raw_data_carbon, "Atom_ID", "Val")
        dat <- dat[order(as.numeric(dat$Seq_ID)), ]
        cs_df <- data.frame(dat[,c("CA", "CB")])
        seq <- dat[,"Comp_ID"]
        seq <- stringr::str_replace_all(seq, c(
                "ALA"="A", "ARG"="R", "ASN"="N", "ASP"="D",
                "CYS"="C", "GLU"="E", "GLN"="Q", "GLY"="G",
                "HIS"="H", "ILE"="I", "LEU"="L", "LYS"="K",
                "MET"="M", "PHE"="F", "PRO"="P", "SER"="S",
                "THR"="T", "TRP"="W", "TYR"="Y", "VAL"="V"
        ))
        return(list(seq, cs_df))
}

# Extracts data from exisitng database entry.
#' \code{read_db_file()} reads in data from existing database that included in the BaMORC package. This database was extracted from RefDB database.
#'
#' @param id BMRB or RefDB entry ID.
#' @importFrom RBMRB fetch_entry_chemical_shifts
#' @importFrom tidyr spread
#' @importFrom stringr str_replace_all
#' @return Protein sequence, secondary structure information and chemical shifts dataframe.
#' @export read_db_file
#'
#' @examples
#' id = 4022
#' head(read_db_file(id))

read_db_file <- function(id){
        raw_data <- RBMRB::fetch_entry_chemical_shifts(id)[,c("Seq_ID", "Comp_ID", "Atom_ID", "Val")]
        raw_data_carbon <- raw_data[which(raw_data$Atom_ID=="CA" | raw_data$Atom_ID=="CB"),]

        dat <- tidyr::spread(raw_data_carbon, "Atom_ID", "Val")
        dat <- dat[order(as.numeric(dat$Seq_ID)), ]
        cs_df <- data.frame(dat[,c("CA", "CB")])
        seq <- dat[,"Comp_ID"]
        seq <- stringr::str_replace_all(seq, c(
                "ALA"="A", "ARG"="R", "ASN"="N", "ASP"="D",
                "CYS"="C", "GLU"="E", "GLN"="Q", "GLY"="G",
                "HIS"="H", "ILE"="I", "LEU"="L", "LYS"="K",
                "MET"="M", "PHE"="F", "PRO"="P", "SER"="S",
                "THR"="T", "TRP"="W", "TYR"="Y", "VAL"="V"
        ))
        return(list(seq, cs_df))
}

#' Using JPred Mass-submission scheduler program to submit protein sequence and return secondary structure results.
#'
#' @param protein_sequence protein sequence
#' @importFrom utils installed.packages
#' @importFrom utils untar
#' @importFrom utils read.csv2
#' @importFrom utils download.file
#' @return protein secondary structure information
#' @export jpred_fetcher
#'
#' @examples
#' protein_sequence <- "MQVWPIEGIKKFETLSYLPPLTVEDLLKQI"
#' \dontrun{secondary_structure <- jpred_fetcher(protein_sequence)}

jpred_fetcher <- function(protein_sequence){
        if(is.null(protein_sequence) | protein_sequence==""){
                stop("Please provide protein sequence.")
        }
        host = "http://www.compbio.dundee.ac.uk/jpred4/cgi-bin/rest/job"
        
        query <- paste("skipPDB=on", "format=seq", paste0(">query\n", protein_sequence), sep = "\u00a3\u20ac\u00a3\u20ac")
        response <- httr::POST(host, body = query, httr::add_headers("Content-type" = "text/txt"))
        
        if(response$status_code == 202 & grepl(pattern = "created jpred job", tolower(httr::content(response, "text")))){
                result_url <- httr::headers(response)$location
                jobid <- stringr::str_match(string = result_url, pattern = "(jp_.*)$")[2]
        }
        else{
                stop(paste(message(httr::content(response, "text")), response$status_code))
        }
        
        job_url <- paste(host, "id", jobid, sep = "/")
        result_url <- httr::headers(response)$location

        
        for(i in c(0:60)) {
                # check whether it's done or not
                if (grepl(pattern = "finished", tolower(httr::content(httr::GET(job_url), "text")))){
                        
                        # Download results
                        archive_url <- paste("http://www.compbio.dundee.ac.uk/jpred4", "results", jobid, paste0(jobid, ".tar.gz"), sep = "/") # where file located on-line
                        dir.create("temp_jpred")
                        jpred_saved_loc <- file.path(paste0("temp_jpred/", jobid, ".tar.gz")) # declare where to save
                        download.file(archive_url, jpred_saved_loc)
                        jpred_target_file <- paste0(jobid, ".jnet") # declare what file to get in the tar arch
                        untar(jpred_saved_loc, files=jpred_target_file) # save the file in root
                        
                        # Remove downloaded tar file
                        unlink("temp_jpred", recursive=TRUE)
                        
                        # Need to process the Secondary Structure result
                        secondary_structure <- read.csv2(jpred_target_file, header = F)[[1]][1]
                        secondary_structure <- gsub("jnetpred:|,", "", secondary_structure)
                        secondary_structure <- gsub("E", "B", secondary_structure)
                        secondary_structure <- gsub("-", "C", secondary_structure)
                        
                        # Remove jnet file after get the secondary structure.
                        garbage_file <- list.files(pattern = ".jnet")
                        file.remove(garbage_file)
                        
                        message(paste("The predicted secondary structure is:", secondary_structure))
                        
                        return(secondary_structure)
                }
                # if not done, wait
                else {
                        Sys.sleep(60)
                }
        }
        
        
}
