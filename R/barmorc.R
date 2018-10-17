#' Calculates the referencing correction value.
#'
#' The core package function that calculates the carbon-13 reference correction using an input protein sequence with associated secondary structure information along with a table of alpha and beta carbon chemical shift  pairs. The output of this function is the correction value that should be added to the input carbon chemical shifts.
#'
#'
#' @name bamorc
#' @param sequence string of sequence
#' @param secondary_structure string of secondary structure (optional)
#' @param chemical_shifts_input table of alpha and beta carbon chemical shift pairs in data.frame
#' @param from the lower bound of the optimization search window
#' @param to the upper bound of the optimization search window
#' @importFrom stats na.omit
#' @return Carbon-13 reference correction value that should be applied (added) to the input carbon chemical shifts data.
#'
#' @export bamorc
#'
#' @examples
#' sequence <- paste(RefDB_data$carbonDat[[1]]$AA,collapse = "")
#' secondary_structure <- paste(RefDB_data$carbonDat[[1]]$SS,collapse = "")
#' chemical_shifts_input <- RefDB_data$carbonDat[[1]][,c(4,5)]
#' from=-5
#' to=5
#'
#' bamorc(sequence, secondary_structure, chemical_shifts_input, from=-5, to=5)
#'
#' # Expected output
#' # [1] 0.0142443
#'
#'
#' sequence <- paste(BaMORC::RefDB_data$carbonDat[[1]]$AA,collapse = "")
#' chemical_shifts_input <- BaMORC::RefDB_data$carbonDat[[1]][,c(4,5)]
#' from=-5
#' to=5
#'
#'
#' \dontrun{bamorc(sequence=sequence, chemical_shifts_input=chemical_shifts_input, from=-5, to=5)}
#' # Expected output
#' # [1] 0.009805279

bamorc <- function(sequence, secondary_structure=NULL, chemical_shifts_input, from=-5, to=5){
        # calculating relative cummulateive frequency

        if(is.null(secondary_structure)){
                secondary_structure <-  BaMORC::jpred_fetcher(sequence)
        }

        secondary_structure <- gsub("U", "C", secondary_structure)
        temp_result <- calculate_RCF(sequence, secondary_structure)

        AA_SS <- temp_result[[2]]$AA_SS
        actual_RCF <- temp_result[[2]]$Freq
        AA <- strsplit(sequence, "")[[1]]
        SS <- strsplit(secondary_structure, "")[[1]]


        chemical_shifts <- data.frame(AA=AA, SS=SS, CA=chemical_shifts_input[,1], CA=chemical_shifts_input[,2])
        # Remove NA and remove G
        if("G" %in% chemical_shifts$AA){
                chemical_shifts <- chemical_shifts[-which(AA=="G"),]
        }

        if(any(is.na(chemical_shifts))){
                chemical_shifts <- na.omit(chemical_shifts)
        }
        label_AA_SS <- paste(chemical_shifts$AA, chemical_shifts$SS, sep = "-")
        fn <- function(x){
                # Debug
                ca_x=0
                cb_x=0
                i=1
                ca_x <- x[1]
                cb_x <- x[2]
                InvDet_v <- sapply(AA_SS, function(aass){
                        det(BaMORC::inversedMatrices$getInvMatrix(aass))
                })
                est_prob <- unlist(lapply(c(1:nrow(chemical_shifts)), function(i){
                        cacb_cs <- unlist(c(chemical_shifts[i,3] + ca_x, chemical_shifts[i,4] + cb_x))
                        aa <- as.character(unlist(chemical_shifts[i,1]))
                        ss <- as.character(unlist(chemical_shifts[i,2]))
                        inv_matrices <- BaMORC::inversedMatrices$getInvMatrix(paste(aa, ss, sep = "-"))

                        chi_str <- t(as.matrix(cacb_cs) - as.matrix(c(BaMORC::CAMuTable[BaMORC::CAMuTable$Residue==aa, ss], BaMORC::CBMuTable[BaMORC::CBMuTable$Residue==aa, ss]))) %*%
                                inv_matrices %*%
                                (as.matrix(cacb_cs) - as.matrix(c(BaMORC::CAMuTable[BaMORC::CAMuTable$Residue==aa, ss], BaMORC::CBMuTable[BaMORC::CBMuTable$Residue==aa, ss])))
                        den <-  calculate_AA_Prob(chi_str)
                        if(aa == "C"){
                                aa <- "B"
                                ss <- unlist(chemical_shifts[i,2])
                                inv_matrices <- BaMORC::inversedMatrices$getInvMatrix(paste(aa, ss, sep = "-"))
                                chi_str <- t(as.matrix(cacb_cs) - as.matrix(c(BaMORC::CAMuTable[BaMORC::CAMuTable$Residue==aa, ss], BaMORC::CBMuTable[BaMORC::CBMuTable$Residue==aa, ss]))) %*%
                                        inv_matrices %*%
                                        (as.matrix(cacb_cs) - as.matrix(c(BaMORC::CAMuTable[BaMORC::CAMuTable$Residue==aa, ss], BaMORC::CBMuTable[BaMORC::CBMuTable$Residue==aa, ss])))
                                if(calculate_AA_Prob(chi_str) > den) {
                                        den <- calculate_AA_Prob(chi_str)
                                        InvDet_v[i] <- det(inv_matrices)
                                }

                        }

                        return(den)
                }))
                est_prob <- data.frame(AA_SS=label_AA_SS, Prob=est_prob)
                est_prob_reduced <- sapply(AA_SS, function(x){
                        return(sum(est_prob[which(est_prob$AA_SS==x),]$Prob))
                })
                mse <- sum((est_prob_reduced - actual_RCF)^2)/length(est_prob)
                return(mse)

        }

        est.ras50 <- DEoptim::DEoptim(fn, lower=c(from,from), upper=c(to,to), control=list(storepopfrom=1, trace=FALSE, itermax=50, CR=1))
        #results <- outer(c(seq(from, to, length.out = 50)), c(seq(from, to, length.out = 50)), Vectorize(fn))
        # print(Sys.time() - t)

        return(est.ras50$optim$bestval)
}



#' Calculates an amino acid typing probability.
#'
#' Function returns the probability (density) for a certain type of amino acid based on a chi-squared statistics wtih 2 degrees of freedom.
#'
#' @param chi_squared_stat a single or a vector of chi-squared statistics
#' @param df degrees of freedom, default is 2
#' @importFrom stats dchisq
#' @return Input can be a single value or a vector of values, the output will be probability density for each value.
#' @examples
#' # Find density for a chi square parameter with 3 degrees of freedom
#' calculate_AA_Prob(0.314, df=3)
#' # Find density for a list of (chi square statistics) with 2 degrees of freedom
#' calculate_AA_Prob(c(0.05, 0.1, 0.5), 2)
#' @export calculate_AA_Prob

calculate_AA_Prob <- function(chi_squared_stat, df=2){
        return(dchisq(chi_squared_stat, df=df))
}


#' Calculates a chi squared statistic(s).
#'
#' \code{calculate_chi_squared_stat} Given a pair of C_alpha and C_beta chemical shifts, this function will return a list of caculated chisquare statistics based on the combination of amino acid typings and secondary structures.
#'
#' @param cacb_pair A pair of carbon chemical shifts \code{c(Ca, Cb)}
#'
#' @return A list of chisquare statisitcs basing the combination of amino acid typings and secondary structures.
#' @examples
#' calculate_chi_squared_stat(c(54,45))
#' @export calculate_chi_squared_stat
#'
#'
calculate_chi_squared_stat <- function(cacb_pair) {
        chiSquaredStat.v <- unlist(lapply(BaMORC::inversedMatrices$getNames(), function(name){
                aa <- strsplit(name, "-")[[1]][1]
                ss <- strsplit(name, "-")[[1]][2]
                chi_str <- t(as.matrix(cacb_pair) - as.matrix(c(BaMORC::CAMuTable[BaMORC::CAMuTable$Residue == aa, ss], BaMORC::CBMuTable[BaMORC::CBMuTable$Residue == aa, ss]))) %*% BaMORC::inversedMatrices$getInvMatrix(name) %*% (as.matrix(cacb_pair) - as.matrix(c(BaMORC::CAMuTable[BaMORC::CAMuTable$Residue == aa, ss], BaMORC::CBMuTable[BaMORC::CBMuTable$Residue == aa, ss])))
                cacb_pair <- rev(cacb_pair)
                chi_str_rev <- t(as.matrix(cacb_pair) - as.matrix(c(BaMORC::CAMuTable[BaMORC::CAMuTable$Residue == aa, ss], BaMORC::CBMuTable[BaMORC::CBMuTable$Residue == aa, ss]))) %*% BaMORC::inversedMatrices$getInvMatrix(name) %*% (as.matrix(cacb_pair) - as.matrix(c(BaMORC::CAMuTable[BaMORC::CAMuTable$Residue == aa, ss], BaMORC::CBMuTable[BaMORC::CBMuTable$Residue == aa, ss])))
                return(min(chi_str, chi_str_rev))
        }))
        #print(chi_str.v)

        chiSquaredStat <- data.table::data.table(cbind(BaMORC::aaCodes1Letter20, matrix(chiSquaredStat.v, ncol = 3, nrow = 20, byrow = T)))
        colnames(chiSquaredStat) <- c("AA", "H", "B", "C")
        return(chiSquaredStat)
}


#' Calculates the relative cummulative frequency for amino acid and secondary structure.
#'
#' This function calculates the relative cummulative freqeucny of amino acid and secondary structure combination.
#'
#' @param sequence String of protein sequence with one letter convention
#' @param secondary_structure String of protein secondary structure with single letter convention
#'
#' @return Relative cummulative frequency.
#' @export calculate_RCF
#'
#' @examples
#' sequence = paste(RefDB_data$carbonDat[[1]]$AA, collapse = "")
#' secondary_structure = paste(RefDB_data$carbonDat[[1]]$SS, collapse = "")
#' relativeCummulativeFrequency = calculate_RCF(sequence, secondary_structure)

calculate_RCF <- function(sequence, secondary_structure){
        # Check whether same length
        if(nchar(secondary_structure) !=nchar(sequence)){
                stop("Sequence and secondary structure must have the same length!")
        }
        # Convert all to upper case
        sequence <- toupper(sequence)
        secondary_structure <- toupper(secondary_structure)

        # Check whether other characters in the sequence or secondary structure.

        standard_AA_code = c("G", BaMORC::aaCodes1Letter19)
        aa <- strsplit(sequence, "")[[1]]
        ss <- strsplit(secondary_structure, "")[[1]]

        # check whether the canonical symbels are used.
        if(!(all(aa %in% standard_AA_code) & all(ss %in% c("B", "H", "C")))){
                stop("Unexpected amino acid or secondary structure symbels included in the input, only 20 standard amino acids and 3 secondary structure one letter symbels are accepted.")
        }

        # remove glycine.
        glycin_Index <- grep("G", aa)
        if(length(glycin_Index) > 0){
                aa <- aa[-glycin_Index]
                ss <- ss[-glycin_Index]
        }

        AA_SS <- paste(aa, ss, sep = "-")
        relative_cumulative_freq <- data.frame(table(AA_SS))

        relative_cumulative_freq$Freq <- relative_cumulative_freq$Freq/sum(relative_cumulative_freq$Freq)
        return(list(AA_SS, relative_cumulative_freq))
}


#' Calculates mean squared error
#'
#' This function will return a mean squared error between estimated density of amino acid typing and secondary structure combination based on the given dataset and reference correction values for the alpha and beta carbons. The estimated amino acid typing density is based on the BaMORC method.
#'
#' @param step_ca Potential correction value for alpha carbon.
#' @param step_cb Potential correction value for beta carbon.
#' @param dat_cacb Chemical shift data fram of alpha and beta carbons.
#' @param aa_Freq Actual amino acid typing and secondary structure frequency calculated basing on provided protein sequence.
#'
#' @return Mean squared error.
#' @examples
#' # chemicalShifts and aaFreq are predefined sample variables for demo purpose.
#'
#' calculate_MSE(step_ca=1, step_cb=1, dat_cacb=chemicalShifts[, c(3,4)], aa_Freq=aaFreq)
#'
#' @export calculate_MSE

calculate_MSE <- function(step_ca, step_cb, dat_cacb, aa_Freq) {

        sum_AA_Prob_v <- apply(
                do.call(cbind, lapply(c(1:nrow(dat_cacb)), function(i){

                        cacb_cs <- unlist(c(dat_cacb[i,1] + step_ca, dat_cacb[i,2] + step_cb))
                        chiStar_v <- as.numeric(as.vector(t(BaMORC::calculate_chi_squared_stat(cacb_cs)[,c(2:4)])))
                        density_v <- BaMORC::calculate_AA_Prob(chiStar_v)
                        names(density_v) <- BaMORC::cname
                        # remove B (cyctine)
                        density_v[grep("C-", BaMORC::cname)] <- density_v[grep("C-", BaMORC::cname)] + density_v[grep("B-", BaMORC::cname)]
                        density_v <- density_v[-grep("B-", BaMORC::cname)]
                        return(density_v)
                })),
                1,
                sum)*BaMORC::AA57OLWeights

        sum_AA_Prob_v <- sum_AA_Prob_v/sum(sum_AA_Prob_v)

        actual_freq_v <- aa_Freq[,2]
        actual_freq_v <- actual_freq_v %*% as.matrix(BaMORC::AA57OLMatrix)
        actual_freq_v <- actual_freq_v*BaMORC::AA57OLWeights
        # actual_freq_v <- actual_freq_v/sum(actual_freq_v)

        mse <- sum((sum_AA_Prob_v- actual_freq_v)^2)/length(sum_AA_Prob_v)

        return(mse)
}
