
#----- Data for Mloglikelihood() and loglikelihood_i ------------------------------------------#
obsdata_creator <- function(formula,data,cluster,strata,frailty,dist){
obsdata <- NULL
#time
if (length(formula[[2]]) == 3) {# --> without left truncation
  obsdata$time <- eval(
    formula[[2]][[2]], 
    envir = data
  )
  obsdata$event <- eval(
    formula[[2]][[3]],
    envir = data
  )
} else if (length(formula[[2]]) == 4) {# --> with left truncation
  obsdata$trunc <- eval(
    formula[[2]][[2]],
    envir = data
  )    
  obsdata$time <- eval(
    formula[[2]][[3]] ,
    envir = data
  )
  obsdata$event <- eval(
    formula[[2]][[4]], 
    envir = data
  )
}
if (!all(levels(as.factor(obsdata$event)) %in% 0:1)) {
  stop(paste("The status indicator 'event' in the Surv object",
             "in the left-hand side of the formula object",
             "must be either 0 (no event) or 1 (event)."))
}

#covariates (an intercept is automatically added)
obsdata$x <- as.data.frame(model.matrix(formula, data=data))

#cluster
if (is.null(cluster)) {
  if (frailty != "none") {
    stop(paste("if you specify a frailty distribution,\n",
               "then you have to specify the cluster variable as well"))
  } else {
    obsdata$cluster <- rep(1, nrow(data))
  }
  #number of clusters
  obsdata$ncl <- 1
  #number of events in each cluster
  obsdata$di <- sum(obsdata$event)
} else {
  if (! cluster %in% names(data)) {
    stop(paste0("object '", cluster, "' not found"))
  }
  obsdata$cluster <- eval(#parse(text=paste0("data$", 
    as.name(cluster),
    envir = data #))
  )
  #number of clusters
  obsdata$ncl <- length(levels(as.factor(obsdata$cluster)))
  #number of events in each cluster
  obsdata$di <- aggregate(obsdata$event,
                          by=list(obsdata$cluster), 
                          FUN=sum)[,, drop=FALSE]
  cnames <- obsdata$di[,1]
  obsdata$di <- as.vector(obsdata$di[,2])
  names(obsdata$di) <- cnames
}

#strata
if (is.null(strata)) {
  obsdata$strata <- rep(1, length(obsdata$time))
  #number of strata
  obsdata$nstr <- 1
  #number of events in each stratum
  obsdata$dq <- sum(obsdata$event)
} else {
  if (!strata %in% names(data)) {
    stop(paste0("object '", strata, "' not found"))
  }
  obsdata$strata <- eval(
    as.name(strata), 
    envir = data
  )
  #number of strata
  obsdata$nstr <- length(levels(as.factor(obsdata$strata)))
  #number of events in each stratum
  obsdata$dq <- aggregate(obsdata$event, 
                          by = list(obsdata$strata), 
                          FUN = sum)[, , drop = FALSE]
  snames <- obsdata$dq[,1]
  obsdata$dq <- as.vector(obsdata$dq[,2])
  names(obsdata$dq) <- snames
}

#cluster+strata
if (!is.null(cluster) && !is.null(strata)) {
  #number of events in each cluster for each stratum
  obsdata$dqi <- xtabs(x~Group.1+Group.2, data = aggregate(
    obsdata$event, 
    by = list(obsdata$cluster, obsdata$strata), 
    FUN = sum))
  dimnames(obsdata$dqi) <- list(cluster = dimnames(obsdata$dqi)[[1]], 
                                strata  = dimnames(obsdata$dqi)[[2]])
} else if (!is.null(cluster)) {
  obsdata$dqi <- obsdata$di
} else if (!is.null(strata)) {
  obsdata$dqi <- obsdata$dq
} else {
  obsdata$dqi <- sum(obsdata$event)
}



#----- Dimensions ---------------------------------------------------------#
#nFpar: number of heterogeneity parameters
if (frailty == "none") {
  nFpar <- 0
} else if (frailty %in% c("gamma", "ingau", "possta", "lognormal")) {
  nFpar <- 1
}
obsdata$nFpar <- nFpar

#nBpar: number of parameters in the baseline hazard
if (dist == "exponential") {
  nBpar <- 1
} else if (dist %in% c("weibull", "inweibull", "frechet", "gompertz",
                       "lognormal", "loglogistic")) {
  nBpar <- 2
} else if (dist %in% c("logskewnormal")) {
  nBpar <- 3
}
obsdata$nBpar <- nBpar

#nRpar: number of regression parameters
nRpar <- ncol(obsdata$x) - 1
obsdata$nRpar <- nRpar  
obsdata
}

# This function computes, for each observation, the logarithm of the likelihood function according to Equation (1) in \cite{Munda2012}.
# It returns a vector of length n (log density for each obs). Parameters in imput must have the same structure of Mloglikelihood function
# in the parfm R pkg. Note: I deleted all the code related to left truncation as we do not consider it (yet) in the mixture framework
loglikelihood_i <- function(p,
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
  
  # baseline: from string to the associated function (baselines R script)
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
  cumhaz_i <- cumhaz[,1]
  cumhaz <- aggregate(cumhaz[, 1], by = list(cumhaz[, 2]), 
                      FUN = sum)[, 2, drop = FALSE]
  ### NO FRAILTY
  if (frailty == "none") cumhaz <- sum(cumhaz)
  
  
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
    loghaz_i <- loghaz[,1] 
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
    loghaz_i <- c(apply(cbind(rownames(pars), pars), 1,
                        function(x) {
                          obs$event[obs$strata == x[1]] * (
                            dist(as.numeric(x[-1]), 
                                 obs$time[obs$strata == x[1]],
                                 what = "lh") + 
                              as.matrix(obs$x[
                                obs$strata == x[1], -1, drop = FALSE]
                              ) %*% as.matrix(beta))
                        }))
  }
  
  
  # ---- log[ (-1)^di L^(di)(cumhaz) ]-------------------------------------- #
  logSurv <- NULL
  if (frailty == "gamma") {
    # logSurv <- mapply(fr.gamma, 
    #                   k = obs$di, s = as.numeric(cumhaz[[1]]), 
    #                   theta = rep(theta, obs$ncl), 
    #                   what = "logLT") 
    logSurv_i <- sapply(1:length(cumhaz_i), function(i) fr.gamma(
                      k = obs$event[i], s = cumhaz_i[i], 
                      theta = theta, 
                      what = "logLT")) 
  } else if (frailty == "ingau") {
    # logSurv <- mapply(fr.ingau, 
    #                   k = obs$di, s = as.numeric(cumhaz[[1]]), 
    #                   theta = rep(theta, obs$ncl), 
    #                   what = "logLT") 
    logSurv_i <- sapply(1:length(cumhaz_i), function(i) fr.ingau(
      k = obs$event[i], s = cumhaz_i[i], 
      theta = theta, 
      what = "logLT")) 
  } else if (frailty == "possta") {
    # logSurv <- sapply(1:obs$ncl, 
    #                   function(x) fr.possta(k = obs$di[x], 
    #                                         s = as.numeric(cumhaz[[1]])[x], 
    #                                         nu = nu, Omega = Omega, 
    #                                         what = "logLT",
    #                                         correct = correct))
    logSurv_i <- sapply(1:length(cumhaz_i), function(i) fr.possta(
      k = obs$event[i], s = cumhaz_i[i], 
      nu = nu, Omega = Omega, 
      what = "logLT", correct=correct)) 
  } else if (frailty == "lognormal") {
    # logSurv <- mapply(fr.lognormal, 
    #                   k = obs$di, s = as.numeric(cumhaz[[1]]), 
    #                   sigma2 = rep(sigma2, obs$ncl), 
    #                   what = "logLT")
    logSurv_i <- sapply(1:length(cumhaz_i), function(i) fr.lognormal(
      k = obs$event[i], s = cumhaz_i[i], 
      sigma2=sigma2, 
      what = "logLT")) 
  } else if (frailty == "none") {
    logSurv_i <- mapply(fr.none, s = cumhaz_i, what = "logLT")
  }
  
  # ---- Log likelihood (vector of length n) ------------------------------------------ #
  loglik_i <- loghaz_i + logSurv_i
  attr(loglik_i, "parameters") <- pars
  loglik_i
}


paramaters_transformator <-
  function(estim_par,
           frailty,
           dist,
           obsdata) {
    
    nBpar <- obsdata$nBpar
    nRpar <- obsdata$nRpar
    nFpar <- obsdata$nFpar  
    
  #heterogeneity parameter
  if (frailty %in% c("gamma", "ingau")) {
    theta <- exp(estim_par[1:nFpar])
    sigma2 <- NULL
    nu <- NULL
  } else if (frailty == "lognormal") {
    theta <- NULL
    sigma2 <- exp(estim_par[1:nFpar])
    nu <- NULL
  } else if (frailty == "possta") {
    theta <- NULL
    sigma2 <- NULL
    nu <- exp(-exp(estim_par[1:nFpar]))
  } else if (frailty == "none") {
    theta <- NULL
    sigma2 <- NULL
    nu <- NULL
  }
  
  #baseline hazard parameter(s)
  if (dist == "exponential") {
    lambda <- exp(estim_par[nFpar + 1:obsdata$nstr])
    ESTIMATE <- c(lambda = lambda)
  } else if (dist %in% c("weibull", "inweibull", "frechet")) {
    rho <- exp(estim_par[nFpar + 1:obsdata$nstr])
    lambda <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(rho = rho, lambda = lambda)
  } else if (dist == "gompertz") {
    gamma <- exp(estim_par[nFpar + 1:obsdata$nstr])
    lambda <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(gamma = gamma, lambda = lambda)
  } else if (dist == "lognormal") {
    mu <- estim_par[nFpar + 1:obsdata$nstr]
    sigma <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(mu = mu, sigma = sigma)
  } else if (dist == "loglogistic") {
    alpha <- estim_par[nFpar + 1:obsdata$nstr]
    kappa <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    ESTIMATE <- c(alpha = alpha, kappa = kappa)
  } else if (dist == "logskewnormal") {
    xi <- estim_par[nFpar + 1:obsdata$nstr]
    omega <- exp(estim_par[nFpar + obsdata$nstr + 1:obsdata$nstr])
    alpha <- estim_par[nFpar + 2 * obsdata$nstr + 1:obsdata$nstr]
    ESTIMATE <- c(xi = xi, omega = omega, alpha = alpha)
  }
  
  #regression parameter(s)
  if (nRpar == 0) {
    beta <- NULL
  } else {
    beta <- estim_par[-(1:(nFpar + nBpar * obsdata$nstr))]
    names(beta) <- paste("beta", names(obsdata$x), sep=".")[-1]
  }
  
  #all together
  ESTIMATE <- c(theta = theta,
                sigma2 = sigma2,
                nu = nu,
                ESTIMATE,
                beta = beta)
  ESTIMATE
}