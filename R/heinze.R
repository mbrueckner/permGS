sampleKM <- function(n, fit) {
    wp <- -diff(c(1, fit$surv))
    sample(fit$time, size=n, replace=TRUE, prob=wp)
}

sampleFromKM.Heinze <- function(n, fit, start=0, tmax=NULL, dv=1) {
    if(is.null(tmax)) tmax <- fit$time[length(fit$time)]
    sapply(runif(n, start, 1), function(v) {
        if(v > max(1-fit$surv)) c(tmax, 0)
        else c(min(fit$time[fit$surv <= 1-v]), dv)
    })
}

sampleFromCondKM.Heinze <- function(U, fit, f, tmax=NULL, dv=1) {
    n <- length(U)
    ##f <- approxfun(fit$time, fit$surv, method="constant", yleft=1, rule=2, f=0)
    sampleFromKM.Heinze(n, fit, 1-f(U), tmax, dv)
}

## validate.sampleKM <- function(M) {
##     par <- param.heinze
##     par$ftype <- 2
##     par$maxTime <- 60
##     par$n <- c(50, 50)

##     f <- function() {
##         data <- generateData(par)
##         fit <- survfit(Surv(time, status) ~ 1, data=data)
##         x <- sampleKM(1000, fit)
##         y <- sampleFromKM.Heinze(1000, fit)[1,]

##         c(median(x), median(y))
##     }

##     replicate(M, f())
## }

## validate.cond <- function(M, u=0) {
##     par <- param.heinze
##     par$ftype <- 2
##     par$maxTime <- 60
##     par$n <- c(50, 50)
##     par$cens$trt <- Exp(0.04)
##     par$cens$ctrl <- Exp(0.04)

##     f <- function() {
##         data <- generateData(par)
##         fit <- survfit(Surv(time, 1-status) ~ 1, data=data)

##         f <- approxfun(fit$time, fit$surv, method="constant", yleft=1, rule=2, f=0)
##         x <- sampleFromCondKM.Heinze(rep.int(u, 100), fit, f)[1,]
##         y <- sampleFromKM.Heinze(100, fit)[1,]
##         c(median(x), median(y))
##     }

##     replicate(M, f())
## }

permImpHeinze <- function(B, data, alpha=0.025) {

    n <- nrow(data)
    data <- as.matrix(data)

    obsS <- sum(logrank_trafo(data[,1], data[,2]) * data[,3])

    ## KM estimator from pooled data
    time <- data[,1]
    status <- data[,2]

    g0 <- data[,3] == 0
    g1 <- data[,3] == 1

    ## split data set by treatment groups
    data1 <- data[g0,]
    data2 <- data[g1,]

    tmax <- max(time)
    tmax1 <- max(data1[,1])
    tmax2 <- max(data2[,1])

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

    fS <- approxfun(fitS$time, fitS$surv, method="constant", yleft=1, rule=2, f=0)
    f1 <- approxfun(fit1$time, fit1$surv, method="constant", yleft=1, rule=2, f=0)
    f2 <- approxfun(fit2$time, fit2$surv, method="constant", yleft=1, rule=2, f=0)

    S <- sapply(1:B, function(i) {
        ## permute rows
        pdata <- data[sample.int(n), ]

        T <- pdata[,1]
        C <- data[,1]
        pdelta <- as.logical(pdata[,2])

        ##only impute survival times for censored obs.
        if(!all(pdelta)) {
            ##T[!pdelta] <- -log(runif(sum(!pdelta))) / 0.04 + T[!delta]
            tmp <- sampleFromCondKM.Heinze(T[!pdelta], fitS, fS, tmax, 1)
            T[!pdelta] <- tmp[1,]
            pdelta[!pdelta] <- tmp[2,]
        }

        ## only impute censoring times for uncensored obs.
        v1 <- data[,2] & g0
        v2 <- data[,2] & g1

        ## 1-f1(U) = fit1$surv[v1]
        ## 1-f2(U) = fit2$surv[v2]
        if(any(v1)) {
            ##C[v1] <- pmin(-log(runif(sum(v1))) / 0.04 + C[v1], 60-pdata[v1,4])
            C[v1] <- sampleFromCondKM.Heinze(C[v1], fit1, f1, tmax, 0)[1,]
        }
        if(any(v2)) {
            ##C[v2] <- pmin(-log(runif(sum(v2))) / 0.04 + C[v2], 60-pdata[v2,4])
            C[v2] <- sampleFromCondKM.Heinze(C[v2], fit2, f2, tmax, 0)[1,]
        }

        pY <- pmin(T, C)
        pdelta <- (T <= C) * pdelta

        sum(logrank_trafo(pY, pdelta) * data[,3])
    })

    p <- (sum(S >= obsS) + 1) / (B + 1)

    cv <- quantile(c(obsS, S), probs=1-alpha, names=FALSE) ##sort(S)[B*(1-alpha)]

    sdS <- sd(S[is.finite(S)])

    list(p=p, obsS=obsS, Z=obsS/sdS, cv=cv, std.cv=cv/sdS, stop=(obsS > cv))
}
