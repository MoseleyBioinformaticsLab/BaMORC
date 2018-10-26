#'
#' Calculates the referencing correction value for unassigned protein NMR peaklists.
#' \code{unassigned_bamorc()} will analyze unassigned protein NMR spectra, first groups the peaklist via SSC, then estimates the secondary structure via JPred, finally using BaMORC core function to calculate the reference correction value.
#'
#' @param peakList_file_loc NMR HNCACB file location
#' @param sequence sequence string of protein of interest
#' @param secondary_structure string of secondary structure if available
#' @param from the lower bound of the optimization
#' @param to the upper bound of the optimization
#' @param ssc location of ssc docker image
#' @param para parameter input for ssc function (no need to change)
#' @importFrom stats na.omit
#' @return Reference correction value.
#' @export unassigned_bamorc
#'
#' @examples
#' sequence = "RPAFCLEPPYAGPGKARIIRYFYNAAAGAAQAFVYGGVRAKRNNFASAADALAACAAA"
#' peakList_file_loc = system.file("extdata", "bpti_HNcoCACB.txt", package = "BaMORC")

#' \dontrun{unassigned_bamorc(peakList_file_loc, sequence, secondary_structure=NULL,
#' from=-5, to=5, ssc="moseleybioinformaticslab/ssc",
#' para="--plformat=sparky --stype=HNcoCACB --dims=H,N,CA/CB --rdims=H,N")}
#' # Expected result should be around (due to randomness): 0.0007890328

unassigned_bamorc <- function(peakList_file_loc, sequence, secondary_structure=NULL, from=-5, to=5, ssc="moseleybioinformaticslab/ssc", para="--plformat=sparky --stype=HNcoCACB --dims=H,N,CA/CB --rdims=H,N"){

        secondary_structure <- jpred_fetcher(sequence)

        # ssc
        ## Check docker is installed or not
        if(is.null(system("docker -v"))){
                stop("Please install docker first!")

        }
        message("Running grouping algorithm")
        ## Check ssc is installed or not
        if(is.null(system("docker images -q moseleybioinformaticslab/ssc"))){
                stop("Please install ssc following the installation!")
        }


        # return the results
        peakList_file_loc <- normalizePath(peakList_file_loc)
        #ssc_cmd <- "docker run -v /Users/bill/ssc-master/tests/peaklists/solution_nmr/sparky/bpti_HNcoCACB.txt:/ssc/test.txt -t moseleybioinformaticslab/ssc group --plpath=/ssc/test.txt --plformat=sparky --stype=HNcoCACB --dims=H,N,CA/CB --rdims=H,N --view --result=/ssc/result/"
        ssc_cmd <- paste0("docker run -v ", peakList_file_loc, ":/ssc/test.txt -t ", ssc, " group --plpath=/ssc/test.txt ", para, " --view --result=/ssc/result/")
        message(ssc_cmd)

        # system(command = ssc_cmd)
        group_output <- jsonlite::fromJSON(system(command = ssc_cmd, intern=T))$peaks
        group_output[length(group_output)] <- NULL # remove the last one
        ref_cor_dat <- sapply(group_output, function(x){data.frame(x$dimensions)[3,]})
        indx <- which(lapply(ref_cor_dat, function(x) length(x)) != 2)
        ref_cor_dat[indx] <- NULL
        df <- data.frame(matrix(unlist(ref_cor_dat), nrow=length(ref_cor_dat), byrow=T))


        # Unassigned BaMORC
        RCF <- BaMORC::calculate_rcf(sequence, secondary_structure)

        AA_SS <- RCF[[1]]
        actual_RCF <- RCF[[2]]

        aa_freq <- data.frame(AA_SS=BaMORC::cname, Freq = sapply(BaMORC::cname, function(x){
                freq = 0.0
                if(x %in% actual_RCF$AA_SS){
                        freq = actual_RCF[actual_RCF$AA_SS==x,]$Freq
                }
                return(freq)
        }))
        aa_freq <- aa_freq[-grep("B-", BaMORC::cname),]
        # AA <- strsplit(sequence, "")[[1]]
        # SS <- strsplit(secondary_structure, "")[[1]]

        chemical_shifts <- data.frame(df)
        if(any(is.na(chemical_shifts))){
                chemical_shifts <- na.omit(chemical_shifts)
        }

        fn <- function(x){
                mse <- BaMORC::calculate_MSE(x, x, chemical_shifts, aa_freq)
                return(mse)
        }
        message("Running reference correction algorithm.")
        #t = Sys.time()
        est_ras50 <- DEoptim::DEoptim(fn, lower=from, upper=to, control=list(trace=FALSE, itermax=5, CR=1, c=1))
        #results <- outer(c(seq(from, to, length.out = 50)), c(seq(from, to, length.out = 50)), Vectorize(fn))
        #message(Sys.time() - t)

        return(est_ras50$optim$bestval)

}
