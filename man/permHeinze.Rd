% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permLR.R
\name{permHeinze}
\alias{permHeinze}
\title{Convenience function which calls createPermGS and nextStage to perform fixed sample size permutation test with Heinze method}
\usage{
permHeinze(formula, data, B = 1000, alpha = 0.05, pool = TRUE,
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
Convenience function which calls createPermGS and nextStage to perform fixed sample size permutation test with Heinze method
}
\examples{
T <- rexp(30) ## event times
Z <- rbinom(30, 1, 0.5)  ## treatment assignment
C <- rexp(30) ## drop-out times
data <- data.frame(time=pmin(T,C), status=T<=C, Z=Z)
x <- permHeinze(Surv(time, status) ~ Z, data)
summary(x)
}
