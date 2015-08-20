
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


power.indep.cor <- function(r1, r2, n1, n2, sig.level = 0.05, alternative = "two.sided") {
   
   # add some error handling
   if (alternative != "two.sided" && alternative != "less" && alternative != "greater") {
      stop("alternative must be 'two.sided', 'less', or 'greater'.")
   }

   effSizeQ         <- atanh(r1) - atanh(r2)           # difference of two Fischer-z-transformed r -> effect size q
   s                <- sqrt( 1/(n1 - 3) + 1/(n2 - 3) ) # see G-Power 3.1 manual
   distributionH1   <- effSizeQ/s                      # mean of the effect size distribution assuming an effect
   
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
   return(pwr)
}
