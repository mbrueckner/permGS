#' impute.IPT
#'
#' Impute data according to IPT method. Output is supposed to be passed to permute.IPT
#'
#' @param data matrix as returned by as.matrix(generateData(param))
#'
#' @return matrix containing imputed survival and censoring times (columns 1 and 2), and original treatment indicator (column 3)
#' @references
#' Wang, R., Lagakos, S.~W. and Gray, R.~J. Testing and interval estimation for two-sample survival comparisons with small sample sizes and unequal censoring. Biostatistics, 11(4), 676--692, January 2010.
#'
impute.IPT <- function(data) {
    ## KM estimator from pooled data
    time <- data[,1]
    status <- data[,2]

    g0 <- data[,3] == 0
    g1 <- data[,3] == 1
    tmax <- max(time)

    ## split data set by treatment groups
    data1 <- data[g0,]
    data2 <- data[g1,]

    ## extract variables because indexing is not allowed in Surv function
    time1 <- data1[,1]
    status1 <- data1[,2]

    time2 <- data2[,1]
    status2 <- data2[,2]

    ## pooled KM
    fitS <- survfit(Surv(time, status) ~ 1)

    ## KM for censoring time in each group
    fit1 <- survfit(Surv(time1, 1-status1) ~ 1)
    fit2 <- survfit(Surv(time2, 1-status2) ~ 1)

    T <- time
    C <- time
    delta <- as.logical(status)
    Tdelta <- delta

    ## only impute survival times for censored obs.
    if(!all(delta)) {
        tmp <- sampleFromCondKM(T[!delta], fitS, tmax, 1)
        T[!delta] <- tmp[1,]
        Tdelta[!delta] <- tmp[2,]
    }

    ## only impute censoring times for uncensored obs.
    v1 <- delta & g0
    v2 <- delta & g1
    if(any(v1)) C[v1] <- sampleFromCondKM(C[v1], fit1, tmax, 0)[1,]
    if(any(v2)) C[v2] <- sampleFromCondKM(C[v2], fit2, tmax, 0)[1,]

    matrix(c(T, C, data[,3], Tdelta), nrow=nrow(data), ncol=4)
}

#' permute.IPT
#'
#' Permute survival times after imputation (IPT)
#'
#' @param data matrix as returned by impute.IPT
#' @param pp vector of permuted indices
#' @param index not used
#'
#' @return matrix with time, status, trt columns
#' @references
#' Wang, R., Lagakos, S.~W. and Gray, R.~J. Testing and interval estimation for two-sample survival comparisons with small sample sizes and unequal censoring. Biostatistics, 11(4), 676--692, January 2010.
#' 
permute.IPT <- function(data, pp, index=TRUE) {
    pT <- data[pp, 1]
    tmp <- matrix(c(pmin(pT, data[,2]), (pT <= data[,2]), data[,3]), nrow=nrow(data), ncol=3)

    tmp[,2] <- (pT <= data[,2]) * data[pp, 4]

    ## same rule as in permImpHeinze
    ## eq <- pT == data[,2]
    ## if(any(eq)) tmp[eq, 2] <- data[pp, 4][eq]

    tmp
}
