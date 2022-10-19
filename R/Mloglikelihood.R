################################################################################
#  Minus the log-likelihood                                                    #
################################################################################
#                                                                              #
#  This function computes minus the logarithm of the likelihood function       #
#                                                                              #
#  Its parameters are                                                          #
#   - p      : the parameters vector, in the form                              #
#              c( frailty distribution parameter(s),                           #
#                 baseline hazard parameter(s),                                #
#                 regression parameter(s) )                                    #
#   - obs    : the observed data, in the form                                  #
#              list( time   = event/censoring times,                           #
#                   [trunc  = left truncation times, ]                         #
#                    event   = event indicators,                               #
#                    x       = covariate data.frame, intercept included        #
#                    cluster = cluster ID vector,                              #
#                    ncl     = number of clusters,                             #
#                    di      = vector giving the numbers of events per cluster #
#   - dist   : the baseline hazard distribution name                           #
#   - frailty: the frailty distribution name                                   #
#   - correct  : (only for possta) the correction to use in case of many       #
#                events per cluster to get finite likelihood values.           #
#                When correct!=0 the likelihood is divided by                  #
#                10^(#clusters * correct) for computation,                     #
#                but the value of the log-likelihood in the output             #
#                is the re-adjusted value.                                     #
#   - transform: should the parameters tranformed to their parameter space or  #
#                are they assumed to be already on their scale?                #
#                The first case (TRUE, the default) is for Mll optimization,   #
#                The secand case (FALSE) is used to straightforwardly compute  #
#                the Hessian matrix on the risght parameter scale              #
#                                                                              #
################################################################################
#                                                                              #
#   Date:              December 19, 2011                                       #
#   Last modification: January 31, 2017                                        #
#                                                                              #
################################################################################

Mloglikelihood <- function(p,
                           obs,
                           dist,
                           frailty,
                           correct,
                           transform = TRUE) { 
    # ---- Assign the number of frailty parameters 'obs$nFpar' --------------- #
    # ---- and compute Sigma for the Positive Stable frailty ----------------- #
    
    if (frailty %in% c("gamma", "ingau")) {
        theta <- ifelse(transform, exp(p[1]), p[1])
    } else if (frailty == "lognormal") {
        sigma2 <- ifelse(transform, exp(p[1]), p[1])
    } else if (frailty == "possta") {
        nu <- ifelse(transform, exp(-exp(p[1])), p[1])
        D <- max(obs$dqi)
        Omega <- Omega(D, correct = correct, nu = nu)
    }
    
    
    # ---- Baseline hazard --------------------------------------------------- #
    if (frailty == 'none') obs$nFpar <- 0
    
    # baseline parameters
    if (dist %in% c("weibull", "inweibull", "frechet")) {
        if (transform) {
            pars <- cbind(rho    = exp(p[obs$nFpar + 1:obs$nstr]),
                          lambda = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
        } else {
            pars <- cbind(rho    = p[obs$nFpar + 1:obs$nstr],
                          lambda = p[obs$nFpar + obs$nstr + 1:obs$nstr])
        }
        beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
    } else if (dist == "exponential") {
        if (transform) {
            pars <- cbind(lambda = exp(p[obs$nFpar + 1:obs$nstr]))
        } else {
            pars <- cbind(lambda = p[obs$nFpar + 1:obs$nstr])
        }
        beta <- p[-(1:(obs$nFpar + obs$nstr))]
    } else if (dist == "gompertz") {
        if (transform) {
            pars <- cbind(gamma  = exp(p[obs$nFpar + 1:obs$nstr]),
                          lambda = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
        } else {
            pars <- cbind(gamma  = p[obs$nFpar + 1:obs$nstr],
                          lambda = p[obs$nFpar + obs$nstr + 1:obs$nstr])
        }  
        beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
    } else if (dist == "lognormal") {
        if (transform) {
            pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                          sigma = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
        } else {
            pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                          sigma = p[obs$nFpar + obs$nstr + 1:obs$nstr])
        }
        beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
    } else if (dist == "loglogistic") {
        if (transform) {
            pars <- cbind(alpha = p[obs$nFpar + 1:obs$nstr],
                          kappa = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]))
        } else  {
            pars <- cbind(alpha = p[obs$nFpar + 1:obs$nstr],
                          kappa = p[obs$nFpar + obs$nstr + 1:obs$nstr])
        }
        beta <- p[-(1:(obs$nFpar + 2 * obs$nstr))]
    } else if (dist == "logskewnormal") {
        if (transform) {
            pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                          sigma = exp(p[obs$nFpar + obs$nstr + 1:obs$nstr]),
                          alpha = exp(p[obs$nFpar + 2 * obs$nstr + 1:obs$nstr]))
        } else {
            pars <- cbind(mu    = p[obs$nFpar + 1:obs$nstr],
                          sigma = p[obs$nFpar + obs$nstr + 1:obs$nstr],
                          alpha = p[obs$nFpar + 2 * obs$nstr + 1:obs$nstr])
        }
        beta <- p[-(1:(obs$nFpar + 3 * obs$nstr))]
    }
    rownames(pars) <- levels(as.factor(obs$strata))
    
    # baseline: from string to the associated function
    dist <- eval(parse(text = dist))
    
    
    # ---- Cumulative Hazard by cluster and by strata ------------------------- #
    
    cumhaz <- NULL
    cumhaz <- matrix(unlist(
        sapply(levels(as.factor(obs$strata)),
               function(x) {t(
                   cbind(dist(pars[x, ], obs$time[obs$strata == x], what = "H"
                   ) * exp(as.matrix(obs$x)[
                       obs$strata == x, -1, drop = FALSE] %*% as.matrix(beta)),
                   obs$cluster[obs$strata == x]))
               })), ncol = 2, byrow = TRUE)
    cumhaz <- aggregate(cumhaz[, 1], by = list(cumhaz[, 2]), 
                        FUN = sum)[, 2, drop = FALSE]
    ### NO FRAILTY
    if (frailty == "none") cumhaz <- sum(cumhaz)
    
    # Possible truncation
    if (!is.null(obs$trunc)) {
        cumhazT <- matrix(unlist(
            sapply(levels(as.factor(obs$strata)),
                   function(x) {t(
                       cbind(dist(pars[x, ], obs$trunc[obs$strata == x], what = "H"
                       ) * exp(as.matrix(obs$x)[
                           obs$strata == x, -1, drop = FALSE] %*% as.matrix(beta)),
                       obs$cluster[obs$strata == x]))
                   })), ncol = 2, byrow = TRUE)
        cumhazT <- aggregate(cumhazT[, 1], by = list(cumhazT[, 2]), 
                             FUN = sum)[, 2, drop = FALSE]
        ### NO FRAILTY
        if (frailty == "none") cumhazT <- sum(cumhazT)
    }
    
    # ---- log-hazard by cluster --------------------------------------------- #
    loghaz <- NULL
    if (frailty != "none")  {
        loghaz <- matrix(unlist(
            sapply(levels(as.factor(obs$strata)),
                   function(x) {
                       t(cbind(obs$event[obs$strata == x] * (
                           dist(pars[x, ], obs$time[obs$strata == x],
                                what = "lh") + 
                               as.matrix(obs$x)[
                                   obs$strata == x, -1, drop = FALSE] %*% 
                               as.matrix(beta)),
                           obs$cluster[obs$strata == x]))
                   })), ncol = 2, byrow = TRUE)
        loghaz <- aggregate(loghaz[, 1], by = list(loghaz[, 2]), FUN = sum)[
            , 2, drop = FALSE]
    } else {
        loghaz <- sum(apply(cbind(rownames(pars), pars), 1,
                            function(x) {
                                sum(obs$event[obs$strata == x[1]] * (
                                    dist(as.numeric(x[-1]), 
                                         obs$time[obs$strata == x[1]],
                                         what = "lh") + 
                                        as.matrix(obs$x[
                                            obs$strata == x[1], -1, drop = FALSE]
                                        ) %*% as.matrix(beta)))
                            }))
    }
    
    
    # ---- log[ (-1)^di L^(di)(cumhaz) ]-------------------------------------- #
    logSurv <- NULL
    if (frailty == "gamma") {
        logSurv <- mapply(fr.gamma, 
                          k = obs$di, s = as.numeric(cumhaz[[1]]), 
                          theta = rep(theta, obs$ncl), 
                          what = "logLT") 
    } else if (frailty == "ingau") {
        logSurv <- mapply(fr.ingau, 
                          k = obs$di, s = as.numeric(cumhaz[[1]]), 
                          theta = rep(theta, obs$ncl), 
                          what = "logLT") 
    } else if (frailty == "possta") {
        logSurv <- sapply(1:obs$ncl, 
                          function(x) fr.possta(k = obs$di[x], 
                                                s = as.numeric(cumhaz[[1]])[x], 
                                                nu = nu, Omega = Omega, 
                                                what = "logLT",
                                                correct = correct))
    } else if (frailty == "lognormal") {
        logSurv <- mapply(fr.lognormal, 
                          k = obs$di, s = as.numeric(cumhaz[[1]]), 
                          sigma2 = rep(sigma2, obs$ncl), 
                          what = "logLT")
    } else if (frailty == "none") {
        logSurv <- mapply(fr.none, s = cumhaz, what = "logLT")
    }
    
    ### Possible left truncation
    if (!is.null(obs$trunc)) {
        logSurvT <- NULL
        if (frailty == "gamma") {
            logSurvT <- mapply(fr.gamma, 
                               k = 0, s = as.numeric(cumhazT[[1]]), 
                               theta = rep(theta, obs$ncl), 
                               what = "logLT") 
        } else if (frailty == "ingau") {
            logSurvT <- mapply(fr.ingau, 
                               k = 0, s = as.numeric(cumhazT[[1]]), 
                               theta = rep(theta, obs$ncl), 
                               what = "logLT") 
        } else if (frailty == "possta") {
            logSurvT <- sapply(1:obs$ncl, 
                               function(x) fr.possta(
                                   k = 0, 
                                   s = as.numeric(cumhazT[[1]])[x], 
                                   nu = nu, Omega = Omega, 
                                   what = "logLT",
                                   correct = correct))
        } else if (frailty == "lognormal") {
            logSurvT <- mapply(fr.lognormal, 
                               k = 0, s = as.numeric(cumhazT[[1]]), 
                               sigma2 = rep(sigma2, obs$ncl), 
                               what = "logLT") 
        } else if (frailty == "none") {
            logSurvT <- mapply(fr.none, s = cumhazT, what = "logLT")
        }
    }
    
    
    # ---- Minus the log likelihood ------------------------------------------ #
    Mloglik <- -sum(as.numeric(loghaz[[1]]) + logSurv)
    if (!is.null(obs$trunc)) {
        Mloglik <- Mloglik + sum(logSurvT)
    }
    attr(Mloglik, "cumhaz") <- as.numeric(cumhaz[[1]])
    if (!is.null(obs$trunc)) {
        attr(Mloglik, "cumhazT") <- as.numeric(cumhazT[[1]])
    } else {
        attr(Mloglik, "cumhazT") <- NULL
    }
    attr(Mloglik, "loghaz") <- as.numeric(loghaz[[1]])
    attr(Mloglik, "logSurv") <- logSurv
    if (!is.null(obs$trunc)) {
        attr(Mloglik, "logSurvT") <- logSurvT
    }
    return(Mloglik)
}
################################################################################
################################################################################



################################################################################
# the same as Mloglikelihood, without attributes, to be passed to optimx()     #
################################################################################
optMloglikelihood <- function(p, obs, dist, frailty, correct) {
    res <- Mloglikelihood(p = p, obs = obs, dist = dist, 
                          frailty = frailty, correct = correct)
    as.numeric(res)}
################################################################################
