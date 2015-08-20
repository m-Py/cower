
#' Power analysis for two independent correlations
#'
#' This function determines the power for detecting a significant difference between two independent Pearson correlation coefficients given n for each group. The power is determined for a Fischer z-test. Currently, it does not behave like functions from package `pwr`. It can only be used to determine the power given two r; not n, effect (q) or significance level.
#'
#' @param r1 Value of first Pearson correlation coefficient
#' @param r2 Value of the other Pearson correlation coefficient
#' @param n1 Sample size of group 1
#' @param n2 Sample size of group 2
#' @param sig.level The employed alpha level
#' @param alternative Must be either 'greater', 'less', or 'two.sided'. 'greater' or 'less' will result in one sided tests (to be used if a directional hypothesis exists), two.sided in a two sided test.
#'
#' @return Single value \code{Vector}. The power to detect a difference between the two given correlation coefficients.
#'
#' @references
#'   Fisher RA. \emph{Statistical Methods for Research Workers}. Edinburgh, Scotland: Oliver and Boyd; 1925. Available: http://psychclassics.yorku.ca.
#'
#'   Cohen, J. (1969). \emph{Statistical power analysis for the behavioural sciences}. New York: Academic Press.
#'
#' @author Martin Papenberg \email{martin.papenberg@@hhu.de}
#' @export


power.indep.cor <- function(r1, r2, n1 = NULL, n2 = NULL, power = NULL, sig.level = 0.05, alternative = "two.sided") {
   
   # add some error handling
   if (alternative != "two.sided" && alternative != "less" && alternative != "greater") {
      stop("alternative must be 'two.sided', 'less', or 'greater'.")
   }
   if ( (is.null(n1) && is.null(power)) || (is.null(n2) && is.null(power)) ) {
      stop("either n1 and n2 or power must be given, not both")
   }
   if ( is.null(n1) && is.null(n2) && is.null(power) ) {
      stop("either n1 and n2 or power must be given, not both")
   }

   effSizeQ         <- atanh(r1) - atanh(r2)           # difference of two Fischer-z-transformed r -> effect size q
   
   # prepare which type of analysis to do: determine power or n?
   if (!is.null(power)) {
      wantedPower      <- power
      startPoint <- 3
      i <- 0 # for approximating power
      n1 <- startPoint + 2^i
      n2 <- startPoint + 2^i
      status <- "increase" # start increasing n in power approximation
      steps  <- 0 # how exact should approximation be
   }
   else { 
      wantedPower <- Inf # not relevant if power is to be determined
   }
   
   currentPower     <- 0                               # approximate power   
   bestFitFound     <- FALSE
   
   # loop until power is reached (if n is searched for)
   while (bestFitFound == FALSE) {
      s                <- sqrt( 1/(n1 - 3) + 1/(n2 - 3) ) # see G-Power 3.1 manual
      distributionH1   <- effSizeQ/s                      # mean of the effect size distribution assuming an effect
      
      #### determine power for given N
      if (alternative  == "less") {
         criticalValue  <- qnorm(sig.level)
         pwr            <- 1 - pnorm(criticalValue, mean = distributionH1, sd = 1)
      }
      else if (alternative == "greater") {
         criticalValue  <- qnorm(1-sig.level)
         pwr            <- 1 - pnorm(criticalValue, mean = distributionH1, sd = 1)
      }
      else if (alternative == "two.sided") {
         if (r1 > r2) {
            criticalValue <- qnorm(1-(sig.level/2))
            pwr           <- 1 - pnorm(criticalValue, mean = distributionH1, sd = 1)
         }
         else if (r1 < r2) {
            criticalValue <- qnorm(sig.level/2)
            pwr           <- pnorm(criticalValue, mean = distributionH1, sd = 1)
         }
      }
      if (is.null(power)) { 
         break # need not loop if power is to be determined
      } 
      
      ### CASE: determine n, continue looping until good fit for power is found
      currentPower <- pwr
      i  <- i+1
      if (steps == 10 && currentPower >= power) { # make sure that approximation is rather exact
         bestFitFound <- TRUE
         print(paste("n per condition", n1))
         return(pwr)
      }
      if (currentPower < wantedPower) { # increase n
         if (status == "decrease") {
            i          <- 0
            startPoint <- n1
            steps      <- steps + 1
         }
         status <- "increase"
         n1 <- startPoint + 2^i
         n2 <- startPoint + 2^i
      }
      else if (currentPower >= wantedPower) {
         if (status == "increase") {
            i          <- 0
            startPoint <- n1
            steps      <- steps + 1
         }
         status <- "decrease"
         n1         <- startPoint - 2^i
         n2         <- startPoint - 2^i
      }
   }
   
   print(c(n1, n2))
   return(pwr)
}
