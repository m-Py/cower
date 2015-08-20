# cower

cower is an r-package for conducting power analyses on the comparison of correlation coefficients. Currently, power analysis for independent correlation coefficients, as tested via Fisher's test, is available. Results are tested against <a href ="http://www.gpower.hhu.de/" target="_blank">G-Power 3.1</a> (and are still being tested).

### General usefulness

Studies, in which independent correlation coefficients are compared require very large sample sizes for a reasonable significance test (e.g.: r = 0.5 vs. r = 0.4 requires 787 participants per group for a power of .8 in a one-sided test). It is hence always sensible to consider a within-subject design if correlation coefficients are to be compared. Since this is not always possible, it is however crucial to know the power to detect a significant difference in such a comparison. cower can be used to conduct such power analyses.

### Types of power analysis that supported

- "post-hoc" - determine the power for two given correlation coefficients and two given n
- "a priori" - specify a desired power and two correlation coefficients, to determine the required n
