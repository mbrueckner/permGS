% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permLR.R
\name{permLR}
\alias{permLR}
\title{Convenience function which calls createPermGS and nextStage to perform fixed sample size permutation test without imputation}
\usage{
permLR(formula, data, B = 1000, alpha = 0.05, pool = TRUE,
  type = c("logrank", "Gehan-Breslow", "Tarone-Ware", "Prentice",
  "Prentice-Marek", "Andersen-Borgan-Gill-Keiding", "Fleming-Harrington",
  "Self"))
}
\arguments{
\item{formula}{a formula object, as used by \code{\link{coxph}}, left hand side must be a 'Surv' object, right hand side must only consist of a factor (treatment indicator) and optionally a special strata() term identifying the permutation strata}

\item{data}{a data.frame or list containing the variables in "formula", by default "formula" is evaluated in the parent frame}

\item{B}{number of random permutations (default: 1000)}

\item{alpha}{significance level (default: 0.05)}

\item{pool}{if TRUE impute event times from Kaplan-Meier estimator calculated from pooled data}

\item{type}{logrank weights to be used with coin::logrank_trafo}
}
\value{
An object of class permGS
}
\description{
Convenience function which calls createPermGS and nextStage to perform fixed sample size permutation test without imputation
}
\examples{
## Two-sided permutation test
T <- rexp(100) ## event times
Z <- rbinom(100, 1, 0.5)  ## treatment assignment
C <- rexp(100) ## drop-out times
data <- data.frame(time=pmin(T,C), status=T<=C, Z=Z)
x <- permLR(Surv(time, status) ~ Z, data, alpha=c(0.025, 0.025))
summary(x)
}
