#' impute.heinze
#'
#' Impute data according to Heinze et al. method. Output is supposed to be passed to permute.heinze
#'
#' @param data matrix as returned by as.matrix(generateData(param))
#'
#' @return list containing Kaplan-Meier estimators of censoring and survival distributions and the original data
#' @references
#' Heinze, G., Gnant, M. and Schemper, M. Exact Log-Rank Tests for Unequal Follow-Up. Biometrics, 59(4), December 2003.
#' 
impute.heinze <- function(data) {    
    time <- data[,1]
    status <- data[,2]

    g0 <- data[,3] == 0
    g1 <- data[,3] == 1

    ## split data set by treatment groups
    data1 <- data[g0,]
    data2 <- data[g1,]

    tmax <- max(time)
    ##tmax1 <- max(data1[,1])
    ##tmax2 <- max(data2[,1])

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

    ## interpolation of pooled survival KM
    fS <- approxfun(fitS$time, fitS$surv, method="constant", yleft=1, rule=2, f=0)

    ## interpolation of censoring KM in group 1
    f1 <- approxfun(fit1$time, fit1$surv, method="constant", yleft=1, rule=2, f=0)

    ## interpolation of censoring KM in group 2
    f2 <- approxfun(fit2$time, fit2$surv, method="constant", yleft=1, rule=2, f=0)

    list(fS=fS, f1=f1, f2=f2, fit1=fit1, fit2=fit2, tmax=tmax, g0=g0, g1=g1, data=data)
}

#' permute.heinze
#'
#' Perform single imputation and permutation step
#'
#' @param imp list as returned by impute.heinze
#' @param pp vector of permuted indices
#' @param index not used
#'
#' @return matrix with time, status, trt columns
#' @references
#' Heinze, G., Gnant, M. and Schemper, M. Exact Log-Rank Tests for Unequal Follow-Up. Biometrics, 59(4), December 2003.
#' 
permute.heinze <- function(imp, pp, index=TRUE) {
    ## permute rows
    pdata <- imp$data[pp, ]
    
    T <- pdata[,1]
    C <- imp$data[,1]
    pdelta <- as.logical(pdata[,2])
    
    ## only impute survival times for censored obs.
    if(!all(pdelta)) {
        tmp <- sampleFromCondKM(T[!pdelta], imp$fitS, imp$tmax, 1, imp$fS)
        T[!pdelta] <- tmp[1,]
        pdelta[!pdelta] <- tmp[2,]
    }
    
    ## only impute censoring times for uncensored obs.
    v1 <- imp$data[,2] & imp$g0
    v2 <- imp$data[,2] & imp$g1
    
    ## 1-f1(U) = fit1$surv[v1]
    ## 1-f2(U) = fit2$surv[v2]
    if(any(v1)) C[v1] <- sampleFromCondKM(C[v1], imp$fit1, imp$tmax, 0, imp$f1)[1,]
    if(any(v2)) C[v2] <- sampleFromCondKM(C[v2], imp$fit2, imp$tmax, 0, imp$f2)[1,]
    
    pY <- pmin(T, C)
    pdelta <- (T <= C) * pdelta

    matrix(c(pY, pdelta, imp$data[,3]), ncol=3, nrow=length(pY))
}
