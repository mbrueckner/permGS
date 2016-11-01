#' permLR
#'
#' One-sided exact / approximate permutation and asymptotic log-rank test
#'
#' This function performs a standard exact or approximate permutation test
#' which is only valid under the extended null hypothesis of equal survival
#' AND censoring distributions.
#'
#' @param B number of random permutations (only used if type="approximate")
#' @param formula formula specifying data model
#' @param data data.frame or list containing the variables in "formula"
#' @param type if type="exact" performs complete enumeration of all permutations, if type="approximate" draw random permutations, if type="asymptotic" perform asymptotic log-rank test
#'
#' @return A list containing the exact or approximate permutation p-value and the observed test statistic
#' @examples
#' T <- rexp(20)
#' C <- rexp(20)
#' data <- data.frame(time=pmin(T, C), status=(T<=C), trt=rbinom(20, 1, 0.5))
#'
#' # Approximate permutation test using 1000 random permutations
#' x <- permLR(1000, Surv(time, status) ~ trt, data, "approximate")
#'
#' print(paste("Approximate permutation p-value:", x$p))
#'
#' # Exact permutation test
#' y <- permLR(0, Surv(time, status) ~ trt, data, "exact")
#' print(paste("Exact permutation p-value:", y$p))
#'
#' @export
permLR <- function(B, formula, data, type="exact") {
    data <- parseFormula(formula, data)
    
    trt <- data$trt

    ## calculate logrank scores
    x <- logrank_trafo(data$time, data$status)

    ## centered logrank statistic
    obsS <- sum(x * trt)

    if(sum(trt) == 0) return(list(p=1, obsS=0))

    if(type == "exact") { ## complete enumeration of all permutations
        co <- utils::combn(1:length(trt), sum(trt))
        S <- vapply(1:ncol(co), function(i) sum(x[co[,i]]), NA_real_)
        p <- mean(S >= obsS)
    } else if(type == "approximate") { ## Monte-Carlo sampling of permutation distribution
        S <- vapply(1:B, function(b) sum(x * sample(trt)), NA_real_)
        p <- (sum(S >= obsS) + 1) / (B + 1)
    } else if(type == "asymptotic") {
        fit <- survdiff(Surv(time, status) ~ trt, data=data)       
        Z <- -sqrt(fit$chisq) * sign(fit$obs[2] - fit$exp[2])
        p <- 1 - pnorm(Z)
    } else stop(paste("Unknown type:", type))

    list(p=p, obsS=obsS)
}

#' shuffleBlock
#' Permute block preserving group sizes, randomization blocks
#' @param block vector of row indices to be permuted
#' @param strata factor defining strata with block
#' @return random permutation of each stratum within block
shuffleBlock <- function(block, strata=0) {
    len <- length(block)
    urb <- unique(strata)
    if(length(urb) == 1) block[sample.int(len, len)]
    else do.call(c, lapply(urb, function(j) {
        ## equivalent to sample(block[which(rnd.block == j)]) but faster
        spln <- block[which(strata == j)]
        spln[sample.int(length(spln), length(spln))]
    }))
}

#' impute.IPT
#'
#' Impute data according to IPT method. Output is supposed to be passed to permute.IPT
#'
#' @param data matrix as returned by as.matrix(generateData(param))
#'
#' @return matrix containing imputed survival and censoring times (columns 1 and 2), and original treatment indicator (column 3)
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

#' impute.IPZ
#'
#' Impute data according to IPZ method. Output is supposed to be passed to permute.IPZ
#'
#' @param data matrix as returned by as.matrix(generateData(param))
#'
#' @return original data with 4 new columns (V1 and V2) containing the imputed observations
#'
impute.IPZ <- function(data) {
    time <- data[,1]
    status <- data[,2]

    tmax <- max(time)
    
    ## KM estimator from pooled data
    fitS <- survfit(Surv(time, status) ~ 1)

    ## split data set by treatment groups
    data1 <- data[data[,3] == 0,]
    data2 <- data[data[,3] == 1,]

    ## extract variable because indexing is not allowed in Surv function
    time1 <- data1[,1]
    status1 <- data1[,2]

    time2 <- data2[,1]
    status2 <- data2[,2]

    ## KM for censoring time in each group
    fit0 <- survfit(Surv(time1, 1-status1) ~ 1)
    fit1 <- survfit(Surv(time2, 1-status2) ~ 1)

    f <- function(work.data, trt.level, fitS, fitK) {
        U <- work.data[,1]
        delta <- as.logical(work.data[,2])
        Tdelta <- delta
        T <- U

        if(!all(delta)) {
            tmp <- sampleFromCondKM(U[!delta], fitS, tmax, 1)
            T[!delta] <- tmp[1,]
            Tdelta[!delta] <- tmp[2,]
        }

        n <- length(U)

        C <- sampleFromKM(n, fitK, 0, tmax, 0)[1,]

        sel1 <- delta & (U <= C)
        sel2 <- delta & (U > C)
        sel3 <- !delta & (U > C)
        sel4 <- !delta & (T <= C)
        sel5 <- !delta & (U < C) & (C < T)

        time <- data[,1]
        status <- data[,2]

        time1 <- numeric(n)
        status1 <- logical(n)

        if(any(sel1)) {
            time1[sel1] <- U[sel1]
            status1[sel1] <- TRUE
        }

        s <- sel2 | sel3 | sel5
        if(any(s)) {
            time1[s] <- C[s]
            status1[s] <- FALSE
        }

        if(any(sel4)) {
            time1[sel4] <- T[sel4]
            status1[sel4] <- Tdelta[sel4] ##TRUE
        }

        time[data[,3] == trt.level] <- time1
        status[data[,3] == trt.level] <- status1

        matrix(c(time, status), nrow=length(time), ncol=2)
    }

    V1 <- f(data2, 1, fitS, fit0)

    V2 <- f(data1, 0, fitS, fit1)

    ## columns: time, status, trt, entry, id, block, rnd.block, time1, status1, time2, status2
    cbind(data[,1:3], V1, V2)
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

#' permute.IPZ
#'
#' Permute treatment assignment after imputation (IPZ)
#'
#' @param data matrix as returned by impute.IPT
#' @param pZ vector of permuted indices if index is TRUE, else binary vector of treatment assignments
#' @param index indicates if pZ is a vector of indices or a binary vector of treatment assignments
#'
#' @return matrix with time, status, Z columns
#'
permute.IPZ <- function(data, pZ, index=FALSE) {
    Z <- data[,3]
    if(index) pZ <- data[pZ, 3]

    time <- data[,1]
    status <- data[,2]

    if(length(pZ) != length(Z)) browser()
    
    a <- pZ > Z
    b <- pZ < Z
    ## pZ == Z: keep as is
    
    time[a] <- data[a, 6]
    status[a] <- data[a, 7]

    time[b] <- data[b, 4]
    status[b] <- data[b, 5]
    
    matrix(c(time, status, pZ), nrow=nrow(data), ncol=3)
}

#' createPermGS
#' 
#' Create permGS object representing a permutational group-sequential trial.
#'
#' @param B number of random permutations
#' @param restricted if TRUE only permute within strata
#' @param IP.method imputation/permuation method IPZ, IPT or none (default: IPZ)
#' @param imputeData user-supplied imputation function (ignored if IP.method is given)
#' @param permuteData user-supplied permutation function (ignore if IP.method is given)
#' @return object of class permGS
#' @examples
#' ## standard permutation test (no imputation, free permutations)
#' x <- createPermGS(1000, FALSE, "none")
#' ## imputation using IPT method, restricted permutations
#' y <- createPermGS(1000, TRUE, "IPT")
#' @export
createPermGS <- function(B=1000, restricted=TRUE, IP.method="IPZ", imputeData=NULL, permuteData=NULL) {
    if(!is.null(IP.method)) {
        if(IP.method == "IPT") {
            imputeData <- impute.IPT
            permuteData <- permute.IPT
        } else if(IP.method == "IPZ") {
            imputeData <- impute.IPZ
            permuteData <- permute.IPZ
        } else if(IP.method == "none") {
            imputeData <- function(data) data
            permuteData <- NULL
        } else stop(paste("Unknown imputation method:", IP.method))
    } else IP.method <- "user"
    
    x <- list(IP.method=IP.method, imputeData=imputeData, permuteData=permuteData,
              B=B, restricted=restricted, S=matrix(nrow=B, ncol=0), perms=matrix(nrow=0, ncol=B), strata=NULL,
              results=data.frame(cv.l=double(), cv.u=double(), obs=double(), reject=logical(), alpha=double(), p=double()))
    class(x) <- c("permGS", "list")
    x
}

#' summary of permGS object
#' @param object permGS object as returned by \code{\link{createPermGS}}
#' @param ... additional parameters (currently unused)
#' @method summary permGS
#' @return nothing
#' @method summary permGS
#' @export
summary.permGS <- function(object, ...) {
    print(paste("Permutations:", object$B))
    print(paste("Restricted:", object$restricted))
    print(paste("Imputation method:", object$IP.method))
    print(paste("Stage:", nrow(object$results)))
    print(object$results)
}

#' nextStage
#'
#' Imputation permutation group-sequential log-rank test.
#' Random permutations of a block a reused in all later stages. This automatically
#' results in blockwise permutations.
#'
#' @param pgs.obj permGS object as returned by \code{\link{createPermGS}}
#' @param alpha alpha at current stage
#' @param time event times
#' @param status censoring indicator
#' @param trt treatment indicator
#' @param strata stratum indicator
#' @return An updated permGS object.
#' @examples
#' ## Two-stage design with one-sided O'Brien-Fleming boundaries using IPZ method
#' x <- createPermGS(1000, TRUE, "IPZ")
#'
#' t1 <- 9  ## calendar time of interim analysis
#' t2 <- 18  ## calendar time of final analysis
#'
#' T <- rexp(100) ## event times
#' R <- runif(100, 0, 12)  ## recruitment times
#' Z <- rbinom(100, 1, 0.5)  ## treatment assignment
#' C <- rexp(100) ## drop-out times
#' 
#' ## Stage 1 data
#' data.t1 <- data.frame(time=pmin(T, C, max(0, (t1-R))), status=(T<=pmin(C, t1-R)), trt=Z)
#' data.t1 <- data.t1[R <= t1,]
#'
#' ## Stage 2 data
#' data.t2 <- data.frame(time=pmin(T, C, max(0, (t2-R))), status=(T<=pmin(C, t2-R)), trt=Z)
#' data.t2 <- data.t2[R <= t2,] 
#' x <- nextStage(x, 0.00153, data.t1$time, data.t1$status, data.t1$trt, rep.int(1, nrow(data.t1)))
#' summary(x)
#'
#' if(!x$results$reject[1]) {
#'    x <- nextStage(x, alpha=0.025, time=data.t2$time, status=data.t2$status, trt=data.t2$trt,
#'           strata=rep.int(c(1,2), c(nrow(data.t1), nrow(data.t2)-nrow(data.t1))))
#'    summary(x)
#' }
#' @export
nextStage <- function(pgs.obj, alpha, time, status, trt, strata) { ##formula, data, entry=NULL) {
    if(any(pgs.obj$results$reject)) return(pgs.obj)

    data <- data.frame(time=time, status=status, trt=trt, strata=strata)
    ##data <- parseFormula(formula, data)

    n <- nrow(data)
    strata <- data$strata
    
    B <- pgs.obj$B
    index <- 1:B
    stage <- nrow(pgs.obj$results) + 1

    cv.l <- pgs.obj$results$cv.l
    cv.u <- pgs.obj$results$cv.u
    S <- pgs.obj$S
    permuteData <- pgs.obj$permuteData

    ## perform imputation step
    data <- as.matrix(data)
    imp.data <- pgs.obj$imputeData(data)

    ## observed scores
    scores <- logrank_trafo(data[,1], data[,2])

    rejected <- logical(B)

    if(pgs.obj$restricted) {
        if(stage > 1) {
            rejected <- vapply(index, function(i) any(S[i,] > cv.u[stage-1] | S[i,] < cv.l[stage-1]), FALSE)
        }
        ## We only need to shuffle the one new block of the current stage B times
        ## perms = matrix of blockwise permuted 1:n vectors (1 perm = 1 column), B columns
        new.strata <- setdiff(unique(strata), pgs.obj$strata)
        if(length(new.strata) > 0) {
            sb <- strata %in% new.strata
            spln <- which(sb)
            srb <- strata[sb]
            pgs.obj$perms <- rbind(pgs.obj$perms, sapply(index, function(i) shuffleBlock(spln, srb)))
        }
        if(nrow(data) != nrow(pgs.obj$perms)) stop("nextStage: New patients but no new strata!")
        pgs.obj$strata <- unique(strata)
    } else {
        spln <- 1:nrow(data)
        srb <- rep.int(1, nrow(data))
        pgs.obj$perms <- sapply(index, function(i) shuffleBlock(spln, srb))
    }
        
    perms <- pgs.obj$perms
    newS <- numeric(B)
    newS[rejected] <- Inf

    ## calculate linear rank statistic for each permutation (scores are calculated once for observed data)
    if(pgs.obj$IP.method == "none") {
        newS[!rejected] <- vapply(index[!rejected], function(i) sum(scores * data[perms[,i], 3]), NA_real_)
        d <- sum(data[,2])
    } else { ## scores need to be re-calculated for every permutation        
        tmp <- vapply(index[!rejected], function(i) {
            pdata <- permuteData(imp.data, perms[,i], TRUE)
            c(sum(pdata[,2]), sum(logrank_trafo(pdata[,1], pdata[,2]) * pdata[,3]))
        }, c(NA_real_, NA_real_))

        d <- mean(tmp[1,])
        newS[!rejected] <- tmp[2,]
    }

    ## observed linear rank statistic
    obsS <- sum(scores * data[,3])
    
    ## critical value
    if(length(alpha) == 2) {
        alpha.l <- alpha[1]
        alpha.u <- alpha[2]
    } else {
        alpha.l <- 0  ## cv.l=-Inf
        alpha.u <- alpha
    }    
    cv.u <- quantile(c(obsS, newS), probs=1-alpha.u, names=FALSE)
    cv.l <- quantile(c(obsS, newS), probs=alpha.l, names=FALSE)

    ## pvalue
    p <- (sum(newS >= obsS) + 1)/(B + 1)

    reject <- (obsS > cv.u) | (obsS < cv.l)
    
    pgs.obj$results <- rbind(pgs.obj$results, data.frame(cv.u=cv.u, cv.l=cv.l, obs=obsS, reject=reject, alpha=alpha, p=p))
    pgs.obj$S <- cbind(pgs.obj$S, newS)
    colnames(pgs.obj$S) <- NULL
    
        pgs.obj
}
