################################################################################
#  Parametric frailty models fitting                                           #
################################################################################
#                                                                              #
#  This function is the core of the package,                                   #
#   performing estimation and returning results.                               #
#                                                                              #
#  It is also the front-end of the package,                                    #
#   being the only that will be visible to the user.                           #
#                                                                              #
#  Its parameters are                                                          #
#   - formula  : a formula object, with the response                           #
#                on the left of a ~ operator,                                  #
#                and the terms on the right.                                   #
#                The response must be a survival object                        #
#                as returned by the Surv() function;                           #
#   - cluster  : the name of the variable in data containing cluster IDs;      #
#   - strata   : the name of the variable in data containing strata IDs;       #
#   - data     : a data.frame containing all the variables;                    #
#   - inip     : initial values for the baseline hazard & reg parameters;      #
#   - iniFpar  : initial value(s) for the heterogeneity parameter(s);          #
#   - dist     : the baseline hazard;                                          #
#   - frailty  : the frailty distribution;                                     #
#   - method   : the optimisation method (See optim());                        #
#   - maxit    : the mamumum number of iterations (See optim());               #
#   - Fparscale: the scaling value for all the frailty parameter(s) in optim() #
#                Optimisation is performed on Fpar/Fparscale                   #
#   - showtime : is the execution time displayed?                              #
#   - correct  : (only for possta) the correction to use in case of many       #
#                events per cluster to get finite likelihood values.           #
#                When correct!=0 the likelihood is divided by                  #
#                10 ^ (#clusters * correct) for computation,                   #
#                but the value of the log-likelihood in the output             #
#                is the re-adjusted value.                                     #
#                                                                              #
################################################################################
#                                                                              #
#   Date:               December 21, 2011                                      #
#   Last modification : January  31, 2017                                      #
#                                                                              #
################################################################################

parfm <- function(formula,
                  cluster   = NULL,
                  strata    = NULL,
                  data,
                  inip      = NULL,
                  iniFpar   = NULL,
                  dist      =  c("weibull", "inweibull", "frechet", "exponential", 
                                 "gompertz", "loglogistic", "lognormal",
                                 "logskewnormal"),
                  frailty   = c("none", "gamma", "ingau", "possta",
                                "lognormal", "loglogistic"),
                  method    = "nlminb",
                  maxit     = 500,
                  Fparscale = 1,
                  showtime  = FALSE,
                  correct   = 0){
    if (missing(data)) {
        data <- eval(parse(text = paste("data.frame(", 
                                      paste(all.vars(formula), collapse = ", "),
                                      ")")))
    }
    
    #----- Check the baseline hazard and the frailty distribution ---------------#
    dist <- tolower(dist)
    dist <- match.arg(dist)
    
    frailty <- tolower(frailty)
    frailty <- match.arg(frailty)
    if (frailty == "none" &&  !is.null(cluster)) {
        warning(paste0("With frailty='none' the cluster variable '",
                       cluster, "' is not used!"))
    }
    if (frailty == "none" &&  !is.null(iniFpar)) {
        warning("With frailty='none' the argument 'iniFpar' is not used!")
    }
    
    #----- 'Correct' is useless except for frailty="possta" -------------------#
    if (frailty == "possta") {  #Do not exaggerate when setting 'correct' !
        if (10 ^ correct == Inf || 10 ^ -correct == 0) {
            stop("'correct' is too large!")
        }
        if (10 ^ correct == 0 || 10 ^ -correct == Inf) {
            stop("'correct' is too small!")
        }
    } else if (correct != 0) {
        warning(paste0("'correct' has no effect when 'frailty = ", frailty, "'"))
    }
    
    
    obsdata <-
      obsdata_creator(
        formula = formula,
        data = data,
        cluster = cluster,
        strata = strata,
        frailty = frailty,
        dist = dist
      )
    
    nBpar <- obsdata$nBpar
    nRpar <- obsdata$nRpar
    nFpar <- obsdata$nFpar
    
    #----- Initial parameters -------------------------------------------------#
    if (!is.null(inip)) {
        #if they are specified, then 
        # (1) check the dimension,
        # (2) check whether they lie in their parameter space, and
        # (3) reparametrise them so that they take on values on the whole real line
        
        # (1)
        if (length(inip) != nBpar * obsdata$nstr + nRpar) {
            stop(paste("number of initial parameters 'inip' must be", 
                       nBpar * obsdata$nstr + nRpar))
        }
        p.init <- inip
        
        # (2)-(3)
        if (dist %in% c("exponential", "weibull", "inweibull", "frechet", "gompertz")) {
            # 1st initial par: log(lambda), log(rho), or log(gamma)
            if (any(p.init[1:obsdata$nstr] <= 0)) {
                stop(paste("with that baseline, the 1st parameter has to be > 0"))
            }
            p.init[1:obsdata$nstr] <- log(p.init[1:obsdata$nstr]) 
        }
        if (dist %in% c("weibull", "inweibull", "frechet", "gompertz", 
                        "lognormal", "loglogistic", "logskewnormal")) {
            #2nd initial par: log(lambda), log(lambda), 
            #                 log(sigma), log(kappa), or log(omega)
            if (any(p.init[obsdata$nstr + 1:obsdata$nstr] <= 0)) {
                stop(paste("with that baseline, the 2nd parameter has to be > 0"))
            }
            p.init[obsdata$nstr + 1:obsdata$nstr] <- 
                log(p.init[obsdata$nstr + 1:obsdata$nstr]) 
        }
    } else {
        # sink('NUL')
        inires <- optimx(par = rep(0, nRpar + nBpar),
                         fn = optMloglikelihood, method = method,
                         obs = obsdata, dist = dist, frailty = 'none',
                         correct = correct,
                         hessian = FALSE,
                         control = list(maxit = maxit,
                                        starttests = FALSE,
                                        dowarn = FALSE))
        # sink()
        p.init <- inires[1:(nRpar + nBpar)]
        # rm(inires)
    }
    
    # --- frailty parameters initialisation --- #
    if (frailty == "none") {
        pars <- NULL
    } else if (frailty %in% c("gamma", "ingau", "lognormal")) {
        if (is.null(iniFpar)) {
            iniFpar <- 1
        } else if (iniFpar <= 0) {
            stop("initial heterogeneity parameter (theta) has to be > 0")
        }
        pars <- log(iniFpar)
    } else if (frailty == "possta") {
        if (is.null(iniFpar)) {
            iniFpar <- 0.5
        } else if (iniFpar <= 0 || iniFpar >= 1) {
            stop("initial heterogeneity parameter (nu) must lie in (0, 1)")
        }
        pars <- log(-log(iniFpar))
    }
    
    pars <- c(pars, unlist(p.init))
    res <- NULL
    
    #--------------------------------------------------------------------------#
    #----- Minimise Mloglikelihood() ------------------------------------------#
    #--------------------------------------------------------------------------#
    # sink('NUL')
    todo <- expression({
        res <- optimx(par = pars, fn = optMloglikelihood, method = method,
                      obs = obsdata, dist = dist, frailty = frailty,
                      correct = correct,
                      hessian = FALSE,
                      control = list(maxit = maxit,
                                     starttests = FALSE,
                                     dowarn = FALSE))
        })
    if (showtime) {
        extime <- system.time(eval(todo))[1]
    } else {
        eval(todo)
        extime <- NULL
    }
    # sink()
    #--------------------------------------------------------------------------#
    
    if (res$convcode > 0) {
        warning("optimisation procedure did not converge,
              conv = ", bquote(.(res$convergence)), ": see ?optimx for details")
    }
    it <- res$niter   #number of iterations
    lL <- -res$value     # value of the marginal loglikelihood
    if (frailty == "possta") {
        lL <- lL + correct * log(10) * obsdata$ncl
    }
    
    
    #--------------------------------------------------------------------------#
    #----- Recover the estimates ----------------------------------------------#
    #--------------------------------------------------------------------------#
    estim_par <- as.numeric(res[1:(nFpar + nBpar  * obsdata$nstr + nRpar)])
    
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
    #--------------------------------------------------------------------------#
    
    
    #--------------------------------------------------------------------------#
    #----- Recover the standard errors ----------------------------------------#
    #--------------------------------------------------------------------------#
    resHessian <- # attr(res, 'details')[1, 'nhatend'][[1]]
        optimHess(par = ESTIMATE, fn = Mloglikelihood, 
                  obs = obsdata, dist = dist, frailty = frailty,
                  correct = correct, transform = FALSE)
    
    var <- try(diag(solve(resHessian)), silent=TRUE)
    if (class(var) == "try-error" | any(is.nan(var))) {
        warning(var[1])
        STDERR <- rep(NA, nFpar + nBpar * obsdata$nstr + nRpar)
        PVAL <- rep(NA, nFpar + nBpar * obsdata$nstr + nRpar)
    } else {
        if (any(var <= 0)) {
            warning(paste("negative variances have been replaced by NAs\n",
                          "Please, try other initial values",
                          "or another optimisation method"))
        }
        
        # heterogeneity (frailty distribution) parameter(s)
        if (frailty %in% c("gamma", "ingau")) {
            seTheta <- sapply(1:nFpar, function(x){
                ifelse(var[x] > 0, sqrt(var[x]), NA)
            })
            seSigma2 <- seNu <- NULL
        } else if (frailty == "lognormal") {
            seSigma2 <- sapply(1:nFpar, function(x){
                ifelse(var[x] > 0, sqrt(var[x]), NA)
            })
            seTheta <- seNu <- NULL
        } else if (frailty == "possta") {
            seNu <- sapply(1:nFpar, function(x){
                ifelse(var[x] > 0, sqrt(var[x]), NA)
            })
            seTheta <- seSigma2 <- NULL
        }
        
        # baseline hazard parameter(s)
        if (dist == "exponential") {
            seLambda <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
            })
            STDERR <- c(seLambda = seLambda)
        } else if (dist %in% c("weibull", "inweibull", "frechet")) {
            seRho <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
            })
            seLambda <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + obsdata$nstr + x] > 0, 
                       sqrt(var[nFpar + obsdata$nstr + x]), NA)
            })
            STDERR <- c(seRho = seRho, seLambda = seLambda)
        } else if (dist == "gompertz") {
            seGamma <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
            })
            seLambda <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + obsdata$nstr + x] > 0,
                       sqrt(var[nFpar + obsdata$nstr + x]), NA)
            })
            STDERR <- c(seGamma = seGamma, seLambda = seLambda)
        } else if (dist == "lognormal") {
            seMu <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
            })
            seSigma <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + obsdata$nstr + x] > 0,
                       sqrt(var[nFpar + obsdata$nstr + x]), NA)
            })
            STDERR <- c(seMu=seMu, seSigma=seSigma)
        } else if (dist == "loglogistic") {
            seAlpha <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
            })
            seKappa <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + obsdata$nstr + x] > 0,
                       sqrt(var[nFpar + obsdata$nstr + x]), NA)
            })
            STDERR <- c(seAlpha=seAlpha, seKappa=seKappa)
        } else if (dist == "logskewnormal") {
            seXi <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + x] > 0, sqrt(var[nFpar + x]), NA)
            })
            seOmega <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + obsdata$nstr + x] > 0,
                       sqrt(var[nFpar + obsdata$nstr + x]), NA)
            })
            seAlpha <- sapply(1:obsdata$nstr, function(x){
                ifelse(var[nFpar + 2 * obsdata$nstr + x] > 0,
                       sqrt(var[nFpar + 2 * obsdata$nstr + x]), NA)
            })
            STDERR <- c(seXi    = seXi, 
                        seOmega = seOmega, 
                        seAlpha = seAlpha)
        }
        
        #regression parameter(s)
        if (nRpar == 0) {
            seBeta <- NULL
        } else {
            seBeta <- numeric(nRpar)
            varBeta <- var[-(1:(nFpar + nBpar * obsdata$nstr))]
            for (i in 1:nRpar) {
                seBeta[i] <- ifelse(varBeta[i] > 0, sqrt(varBeta[i]), NA)
            }
            PVAL <- c(rep(NA, nFpar + nBpar * obsdata$nstr), 
                      2 * pnorm(q = -abs(beta / seBeta)))
        }
        
        #all together
        STDERR <- c(STDERR, se.beta = seBeta)
        if (frailty != "none") {
            STDERR <- c(se.theta  = seTheta,
                        se.sigma2 = seSigma2,
                        se.nu     = seNu,
                        STDERR)
        }
    }
    
    
    
    #--------------------------------------------------------------------------#
    #----- Output -------------------------------------------------------------#
    #--------------------------------------------------------------------------#
    resmodel <- cbind(ESTIMATE = ESTIMATE, 
                      SE       = STDERR)
    rownames(resmodel) <- gsub("beta.","", rownames(resmodel))
    
    if (nRpar > 0)
        resmodel <- cbind(resmodel, "p-val" = PVAL)
    
    class(resmodel) <- c("parfm", class(resmodel))
    
    ### - Checks - #############################################################
    Call <- match.call()
    if (!match("formula", names(Call), nomatch = 0))
        stop("A formula argument is required")
    
    Terms <- terms(formula, data = data)
    ###################################################### - End of Checks - ###
    
    attributes(resmodel) <- c(attributes(resmodel), list(
        call        = Call,
        convergence = res$convergence,
        it          = it,
        extime      = extime,
        nobs        = nrow(data),
        shared      = (nrow(data) > obsdata$ncl),
        loglik      = lL,
        dist        = dist,
        cumhaz      = attributes(Mloglikelihood(p       = estim_par,
                                                obs     = obsdata, 
                                                dist    = dist, 
                                                frailty = frailty,
                                                correct = correct))$cumhaz,
        cumhazT     = attributes(Mloglikelihood(p       = estim_par,
                                                obs     = obsdata,
                                                dist    = dist, 
                                                frailty = frailty,
                                                correct = correct))$cumhazT,
        di          = obsdata$di,
        dq          = obsdata$dq,
        dqi         = obsdata$dqi,
        frailty     = frailty,
        clustname   = cluster,
        stratname   = strata,
        correct     = correct,
        formula     = as.character(Call[match("formula", names(Call), 
                                              nomatch = 0)]),
        terms       = attr(Terms, "term.labels"),
        FisherI     = resHessian,
        estim_par=estim_par
    ))
    if (frailty != "none") {
        names(attr(resmodel, "cumhaz")) <-
            names(attr(resmodel, "di")) <- 
            unique(obsdata$cluster)
    }
    if (showtime)
        cat("\nExecution time:", extime, "second(s) \n")
    
    return(resmodel)
}