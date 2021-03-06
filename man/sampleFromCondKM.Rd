% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{sampleFromCondKM}
\alias{sampleFromCondKM}
\title{sampleFromCondKM}
\usage{
sampleFromCondKM(U, fit, tmax = NULL, dv = 1, f = NULL)
}
\arguments{
\item{U}{vector of observed times}

\item{fit}{Kaplan-Meier fit as returned by survfit}

\item{tmax}{largest observation of the pooled sample}

\item{dv}{1 if imputing events, 0 if imputing censoring times}

\item{f}{interpolated Kaplan-Meier estimate}
}
\value{
Random sample of survival times drawn from conditional distribution of T given T > U
}
\description{
Sample from conditional distribution estimated by Kaplan-Meier estimator. Imputed values > tmax are right-censored.
}
