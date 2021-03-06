#' permGS
#'
#' This package implements permutational group-sequential tests for time-to-event data based on (weighted) log-rank
#' test statistics. It supports exact permutation test when the censoring distributions are equal in the treatment
#' and the control group and the approximate imputation-permutation methods of Heinze et al. (2003) and Wang et al. (2010)
#' and  when the censoring distributions are different. Permutations can be stratified, i.e. only patients within the
#' same stratum are treated as exchangeable. Rejection boundaries are monotone and finite even when only a random
#' subset of all permutations is used. One- and Two-sided testing possible.
#' 
#' @author Matthias Brueckner \email{m.bruckner@@lancaster.ac.uk}, Franz Koenig \email{Franz.Koenig@@meduniwien.ac.at}, Martin Posch \email{martin.posch@@meduniwien.ac.at}
#'
#' @references
#' Brueckner, M., Koenig, F. and Posch, M. Group-sequential permutation tests for time-to-event data.
#' 
#' Heinze, G., Gnant, M. and Schemper, M. Exact Log-Rank Tests for Unequal Follow-Up. Biometrics, 59(4), December 2003.
#'
#' Wang, R., Lagakos, S.~W. and Gray, R.~J. Testing and interval estimation for two-sample survival comparisons with small sample sizes and unequal censoring. Biostatistics, 11(4), 676--692, January 2010.
#'
#' Kelly, P., Zhou, Y., Whitehead, N. J., Stallard, N. and Bowman, C. Sequentially testing for a gene–drug interaction in a genomewide analysis. Statistics in Medicine, 27(11), 2022--2034, May 2008.
#'
#' @examples
#' ## IPZ method based on logrank test with 1000 restricted random permutations
#' x <- createPermGS(1000, TRUE, "IPZ", type="logrank")
#'
#' T <- rexp(100) ## event times
#' R <- runif(100, 0, 12)  ## recruitment times
#' Z <- rbinom(100, 1, 0.5)  ## treatment assignment
#' C <- rexp(100) ## drop-out times
#'
#' ## two-stage design
#' t1 <- 9  ## calendar time of interim analysis
#' t2 <- 18  ## calendar time of final analysis
#'
#' ## Stage 1
#' data.t1 <- data.frame(time=pmin(T, C, max(0, (t1-R))), status=(T<=pmin(C, t1-R)), trt=Z)
#' data.t1 <- data.t1[R <= t1,] 
#' x <- nextStage(x, 0.00153, Surv(time, status) ~ trt, data.t1)
#' summary(x)
#'
#' if(!x$results$reject[1]) { ## Stage 2
#'    data.t2 <- data.frame(time=pmin(T, C, max(0, (t2-R))), status=(T<=pmin(C, t2-R)), trt=Z)
#'    data.t2 <- data.t2[R <= t2,]
#'    data.t2$strata <- rep.int(c(1,2), c(nrow(data.t1), nrow(data.t2)-nrow(data.t1)))
#'    x <- nextStage(x, alpha=0.025, Surv(time, status) ~ trt + strata(strata), data.t2)           
#'    summary(x)
#' }
#'
#' @name permGS
#' @docType package
#' @import survival stats
#' @importFrom utils combn
#' @importFrom coin logrank_trafo
NULL
