
#' Power analysis for two independent correlations
#'
#' This function can be used to conduct power analyses for the
#' comparison of independent correlation coefficients, as tested via
#' Fisher's z-test. It can (1) determine the power for detecting a
#' significant difference between 2 correlations given N ('post hoc'),
#' and (2) determine the required N to achieve a desired power
#' ('a-priori').
#'
#' @param r1 Value of first Pearson correlation coefficient, must be
#'     passed.
#' @param r2 Value of second correlation coefficient, must be passed.
#' @param n1 Sample size of group 1; pass only for 'post hoc' analysis.
#' @param n2 Sample size of group 2; pass only for 'post hoc' analysis.
#' @param power The statistical power; pass only for 'a-priori'
#'     analysis.
#' @param sig.level The employed alpha level, defaults to 0.05.
#' @param alternative Must be 'greater', 'less', or
#'     'two.sided'. 'greater' (r1 > r2) or 'less' (r1 < r2) will result
#'     in one sided tests (to be used if a directional hypothesis
#'     exists). Default test is two.sided.
#'
#' @return A \code{list} that contains all passed and computed values.
#'   \item{r1}{passed correlation coefficient r1}
#'   \item{r2}{passed correlation coefficient r2}
#'   \item{q}{The effect size. Represents the difference between
#'      the z-transformed values of r1 and r2}
#'   \item{n1}{sample size in group 1, either passed ('post hoc'), or computed ('a-priori')}
#'   \item{n2}{sample size in group 1, either passed ('post hoc'), or computed ('a-priori')}
#'   \item{power}{Power to detect a significant difference between
#'      the two correlation coefficients r1 and r2, either passed
#'      ('a-priori') or computed ('post hoc')}
#'   \item{sig.level}{The significance level that was employed for the power analysis}
#'   \item{hypothesis}{Passed argument 'alternative'}
#'
#' @references
#'
#' Fisher RA. \emph{Statistical Methods for Research
#'     Workers}. Edinburgh, Scotland: Oliver and Boyd; 1925. Available:
#'     http://psychclassics.yorku.ca.
#'
#'    Cohen, J. (1969). \emph{Statistical power analysis for the
#'    behavioural sciences}. New York: Academic Press.
#'
#' @author Martin Papenberg \email{martin.papenberg@@hhu.de}
#' 
#' @examples
#'
#' ## A priori analysis
#' power.indep.cor(r1 = 0.4, r2 = 0.3, power = 0.8)
#'
#' ## Post-hoc analysis (one-sided)
#' power.indep.cor(r1 = 0.4, r2 = 0.3, n1 = 450, n2 = 450, alternative = "greater")
#' 
#' @export


power.indep.cor <- function(r1, r2, n1 = NULL, n2 = NULL, power = NULL,
                            sig.level = 0.05, alternative = "two.sided") {
   
    ## add some error handling
    if (alternative != "two.sided" && alternative != "less" && alternative != "greater") {
        stop("Error: alternative must be 'two.sided', 'less', or 'greater'.")
    }
    if ( (!is.null(n1) && !is.null(power)) || (!is.null(n2) && !is.null(power)) ) {
        stop("Error: either n1 and n2 or power must be given, but n and power were passed.")
    }
    if ( is.null(n1) && is.null(n2) && is.null(power) ) {
        stop("Error: either n1 and n2 or power must be given, but none was passed.")
    }
    
    ## Difference of two Fischer-z-transformed r -> effect size q
    effSizeQ <- atanh(r1) - atanh(r2)
   
    ## Post-hoc or a-priori analysis?
    if (is.null(power)) { # post-hoc power analysis
        pwr <- get_power(n1, n2, r1, r2, effSizeQ, alternative, sig.level)
    } else { # a priori power analysis
        values <- a_priori(power, r1, r2, effSizeQ, alternative, sig.level)
        n1  <- values[1]
        n2  <- values[2]
        pwr <- values[3]
    }

    ## return all values in one list
    rtn         <- list()
    rtn[["r1"]] <- r1
    rtn[["r2"]] <- r2
    rtn[["q"]]  <- abs(effSizeQ)
    rtn[["n1"]] <- n1
    rtn[["n2"]] <- n2
    rtn[["power"]] <- pwr
    rtn[["sig.level"]] <- sig.level
    rtn[["hypothesis"]] <- alternative
    return(rtn)
}


a_priori <- function(power, r1, r2, effSizeQ, alternative, sig.level) {

    ## create some book keeping variables
    currentPower <- 0 
    wantedPower  <- power
    startPoint <- 3
    i <- 0 # used to vary n
    ## start n = 4 per condition
    n1 <- startPoint + 2^i # approximate the required n in steps that are a power of 2
    n2 <- startPoint + 2^i
    status <- "increase" # start increasing n in power approximation
    steps  <- 0          # keep track about exactness of approximation
    ## Initialize variable that tracks whether power was achieved using given N:
    bestFitFound     <- FALSE

    ## loop until required power is achieved (if N is searched for)
    while (bestFitFound == FALSE) { 
        ## determine power for current N
        currentPower <- get_power(n1, n2, r1, r2, effSizeQ, alternative, sig.level)
        i  <- i + 1
        # make sure that approximation is rather exact by using many
        # steps
        if (steps >= 10 && currentPower >= wantedPower) {
            bestFitFound <- TRUE
            break ## best N was found!
        }
        if (currentPower < wantedPower) {
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
    return(c(n1, n2, currentPower))
}

## Function that determines achieved power for correlation coefficients
## and n.
get_power <- function(n1, n2, r1, r2, effSizeQ, alternative, sig.level) {

    ## Some error handling
    if (r1 == r2) {
        stop("r1 has the same value as r2, power cannot be calculated")
    } else if (r1 > r2 & alternative == "less") {
        stop("r1 is greater than r2, but argument alternative is 'less'; no power can be calculated")
    } else if (r1 < r2 & alternative == "greater") {
        stop("r1 is less than r2, but argument alternative is 'greater'; no power can be calculated")
    }
    
    s              <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3)) # see G-Power 3.1 manual
    distributionH1 <- effSizeQ / s # mean of the effect size distribution
    if (alternative  == "less") {
        criticalValue <- qnorm(sig.level)
        pwr           <- pnorm(criticalValue, mean = distributionH1, sd = 1)
    }
    else if (alternative == "greater") {
        criticalValue <- qnorm(1-sig.level)
        pwr           <- 1 - pnorm(criticalValue, mean = distributionH1, sd = 1)
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
