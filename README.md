# cower

`cower` is an R-package to conduct power analyses on the comparison of
correlation coefficients. Currently, power analyses for the comparisons
of independent correlations -- as tested via Fisher's z-test -- are
available. Results have been tested against <a href
="http://www.gpower.hhu.de/" target="_blank">G-Power 3.1</a>.

## General purpose

Studies that compare independent correlation coefficients require very
large sample sizes to obtain reasonable power (e.g.: *r* = 0.5 versus
*r* = 0.4 requires 787 participants per group for a power of .80 in a
one-sided test). Therefore, it is often preferable to consider a
comparison of dependent correlations if that is feasible. However, this
is not always possible; it then crucial to know the power to detect a
significant difference in independent correlations. `cower` is used to
conduct such power analyses.

## Types of power analysis that supported

- "post-hoc" - determine the power for two given correlation
  coefficients and given sample size
- "a priori" - specify a desired power and two correlation coefficients
  to determine the required sample size

## Installation

As `cower` is not available from CRAN, you can install it directly from
this GitHub repository. To do so, you need the `devtools` package. Then
run the following commands:

```R

library("devtools") # if not available: install.packages("devtools")
install_github("m-Py/cower")

# load the package via 
library("cower")

```

## "A priori" power analysis

To compute the number of participants needed to obtain a certain power,
we can use the function `power.indep.cor`. We specify two hypothesized
population correlation coefficients and the desired power:

```R
power.indep.cor(r1 = 0.4, r2 = 0.3, power = .8)

$r1
[1] 0.4

$r2
[1] 0.3

$q
[1] 0.1141293

$n1
[1] 1209

$n2
[1] 1209

$power
[1] 0.8002746

$sig.level
[1] 0.05

$hypothesis
[1] "two.sided"

```

By default, the power for a two-sided test is computed, and an alpha
level of .5 is adapted. The alpha-level can be changed using the
parameter `sig.level` and the sidedness can be changed using the
parameter `alternative` (for a one-sided test, set `alternative` to
"less" or "greater", depending on whether r1 is smaller or greater than
r2).

## "Post-hoc" power analysis

To determine the achieved power for the comparisons of two given
population correlation coefficients and sample size, we can use
`power.indep.cor` the following way. Here we do not specify the
parameter `power` (which is to be computed), but instead specify two
sample sizes `n1` and `n2`.


```R
power.indep.cor(r1 = 0.4, r2 = 0.3, n1 = 450, n2 = 350)

$r1
[1] 0.4

$r2
[1] 0.3

$q
[1] 0.1141293

$n1
[1] 450

$n2
[1] 350

$power
[1] 0.3576306

$sig.level
[1] 0.05

$hypothesis
[1] "two.sided"

```

Again, we could also adjust the alpha-level and the sidedness of the
hypothesis test.

## Questions and suggestions

If you have any questions or suggestions (which are greatly
appreciated), just open an issue at Github or contact me via
martin.papenberg at hhu.de.
