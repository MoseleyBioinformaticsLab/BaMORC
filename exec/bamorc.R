#!/usr/bin/env Rscript
"
Usage:
        bamorc.R assigned   (--table=<csv> | --bmrb=<bmrb> | --id=<id>)  [--ppm_range=<range>] [--output=<output_filename>] [--delimiter=<delim>] [--report=<report_filename>]

        bamorc.R unassigned (--table=<csv>) (--seq=<sequence>) [--ppm_range=<range>] [--output=<output_filename>] [--delimiter=<delim>] [--report=<report_filename>] [--ssc=<path>]

        bamorc.R valid_ids

        bamorc.R -h | --help
        bamorc.R -v | --version

Options:
        -t, --table=<csv>                       Input file path
        -s, --seq=<sequence>                    Sequence of Protein of interest, if protein sequence file is not provided.
        -d, --delimiter=<delim>                 Delimiter option can be 'comma', 'tab' and 'whitespace'.
        -p, --ppm_range=<range>                 The ppm search range for reference correction value. [default: '-5,5']
        -i, --id=<id>                           RefDB or BMRB ID.
        -o, --output=<output_filename>          Filename of output of BaMORC result in csv format.
        -g, --ssc=<path>                        Spin system creater. [default: moseleybioinformaticslab/ssc]
        -h, --help                              Show this help message.
        -v, --version                           Show api version
        --bmrb                                  filepath fo BMRB file

        --report=<report_filename>
        --grouped
" -> doc

# checking the first line for table. Check all the values. table 3rd 4rd column numerical.
# unassigned one file only, then need to parse the sequence, generate error if sequence is not available. If sequence in both file and sequence, then check whether they are same.
# for the unassigned, out put the original file with error correct, not ssc one.

#library(docopt)
#library(BaMORC)
#library(BMRBr)

BaMORC_version <- as.character(packageVersion("BaMORC"))
cmdargs        <- docopt::docopt(doc, version = BaMORC_version)


#' BaMORC CLI
#'
#' @param cmdargs List of command-line arguments for BaMORC.R
#'
#' @return BaMORC analytic results
#'
#' @examples

cli <- function(cmdargs){
        message("")
        message("")
        message("")
        message(" *************************************** ")
        message("                                         ")
        message(" * Thanks for choosing BaMORC package! * ")
        message("                                         ")
        message(" *************************************** ")

        if(cmdargs$valid_ids) {
                message(BaMORC::ID)
        }

        if(cmdargs$assigned) {
                if(length(which(c(cmdargs$table, cmdargs$bmrb, cmdargs$id) != FALSE)) != 1){
                        stop("Please privide only one required argument for --table, bmrb, or id!")
                }

                if(cmdargs$ppm_range){
                        ppm_range <- as.numeric(strsplit(cmdargs$ppm_range, ",")[[1]])
                        from <- ppm_range[1]
                        to <- ppm_range[2]
                        if(from >= to){
                                stop("Please correct the searching range for reference correction value.")
                        }
                }


                if(cmdargs$table){
                        input_data <- BaMOARC::read_file(cmdargs$table,
                                                         delim=cmdargs$delimiter,
                                                         assigned=TRUE)
                        corr_value <- BaMOARC::BaMORC(sequence=input_data[,1],
                                                      secondaryStructure=input_data[,2],
                                                      chemicalShifts=input_data[,c(3,4)],
                                                      from=from,
                                                      to=to)
                        output_data <- input_data[,c(3,4)] + corr_value
                }

                else if(cmdargs$bmrb){
                        input_data <- BaMOARC::read_NMRSTAR_file(cmdargs$bmrb)
                        corr_value <- BaMOARC::BaMORC_BMRB(sequence=input_data[[1]],
                                                           chemicalShifts=input_data[[2]][,c(2,3)],
                                                           from=from,
                                                           to=to)
                        output_data<- input_data[[2]][,c(2,3)] + corr_value
                }

                else if(cmdargs$id){
                        if(!(cmdargs$id %in% BaMORC::ID)){
                                stop("please make sure the the dataset ID is correct. To show all the IDs, use 'bamorc.R valid_ids'.")
                        }
                        input <- BaMORC::read_DB_File(cmdargs$id)
                        corr_value <- BaMORC::BaMORC(sequence = paste(input[,1], collapse = ''),
                                                     secondaryStructure = paste(input[,2], collapse = ''),
                                                     chemicalShifts_input = input[,c(3,4)],
                                                     from = from,
                                                     to = to)
                        output_data <- input_data[,c(3,4)] + corr_value
                }

                output_file <- ifelse(cmdargs$report, cmdargs$report, "./BaMORC_output.csv")
                write.csv(output_data, output_file)
                message(paste("The reference correction value is:", corr_value))
                message("The corrected NMR dataset is")
                message(output_data)

        }
        else if(cmdargs$unassigned){
                if(cmdargs$ppm_range){
                        ppm_range <- as.numeric(strsplit(cmdargs$ppm_range, ",")[[1]])
                        from <- ppm_range[1]
                        to <- ppm_range[2]
                        if(from >= to){
                                stop("Please correct the searching range for reference correction value.")
                        }
                }

                corr_value <- BaMORC::UnassignedBaMORC(peakList_file_loc = cmdargs$table,
                                                       sequence = cmdargs$seq,
                                                       secondaryStructure = NULL,
                                                       from = from,
                                                       to = to,
                                                       ssc = cmdargs$ssc,
                                                       para = cmdargs$para)
                message(paste("The reference correction value is:", corr_value))
        }

}

cli(cmdargs)

# To find CLI executable:
# system.file("exec", "bamorkapi", package = "BaMORC")

