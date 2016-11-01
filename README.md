# permGS
permGS is an R package implementing permutational group-sequential tests for time-to-event data based on
the two-sample log-rank test statistic. It supports exact permutation test when the censoring distributions 
are equal in the treatment and the control group and approximate imputation-permutation methods when the 
censoring distributions are different. One- and two-sided testing is possible.

## Installation
```R
# install.packages("devtools")
devtools::install_github("mbrueckner/permGS")
```
