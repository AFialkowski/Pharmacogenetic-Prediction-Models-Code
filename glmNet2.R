# This function is a modified version of glmNet() from Nengjun Yi's BhGLM R package.  It provides the 
# option to scale y or not (default = TRUE).

# Required packages: glmnet, BhGLM

# Reference: Nengjun Yi (2018). BhGLM: Bayesian hierarchical GLMs and survival models, with applications 
# to Genetics and Epidemiology. R package version 1.1.0. https://www.soph.uab.edu/ssg/software/bhglm, 
# https://github.com/nyiuab/BhGLM.

glmNet2 <- function (x, y, family = c("gaussian", "binomial", "poisson", 
    "cox"), weights = rep(1, nrow(x)), offset = NULL, alpha = c(1, 
    0.5, 0), lambda, penalty.factor = rep(1, ncol(x)), nfolds = 10, 
    ncv = 10, verbose = TRUE, scale.y = TRUE) 
{
    require(glmnet)
    start.time <- Sys.time()
    call <- match.call()
    if (any(is.na(x)) | any(is.na(y))) {
        a <- apply(cbind(y, x), 1, function(z) !any(is.na(z)))
        y <- y[a]
        x <- x[a, ]
    }
    x <- as.matrix(x)
    if (NROW(x) != NROW(y)) 
        stop("nobs of 'x' and 'y' are different")
    if (is.null(colnames(x))) 
        colnames(x) <- paste("x", 1:ncol(x), sep = "")
    family <- family[1]
    if (family == "cox") 
        if (!is.Surv(y)) 
            stop("'y' should be a 'Surv' object")
    if (family == "gaussian" & scale.y == TRUE) {
        y <- (y - mean(y))/sd(y)
    }
    alpha <- alpha[1]
    if (length(penalty.factor) != ncol(x)) 
        stop("give each predictor a penalty")
    names(penalty.factor) <- colnames(x)
    penalty.factor <- ncol(x) * penalty.factor/sum(penalty.factor)
    if (nfolds == NROW(x)) 
        ncv <- 1
    if (missing(lambda)) {
        lam <- rep(NA, ncv)
        for (i in 1:ncv) {
            cv.f <- cv.glmnet(x = x, y = y, family = family, 
                weights = weights, offset = offset, alpha = alpha, 
                nfolds = nfolds, grouped = TRUE, penalty.factor = penalty.factor, 
                standardize = FALSE)
            lam[i] <- cv.f$lambda.min
            if (verbose & ncv > 1) {
                pre <- rep("\b", ncv + 1)
                cat(pre, i, "/", ncv, sep = "")
                flush.console()
            }
        }
        lambda <- mean(lam)
    }
    f <- glmnet(x = x, y = y, family = family, weights = weights, 
        offset = offset, alpha = alpha, lambda = lambda, penalty.factor = penalty.factor, 
        standardize = FALSE)
    f$coefficients <- as.numeric(coef(f))
    names(f$coefficients) <- rownames(coef(f))
    f$linear.predictors <- predict(f, newx = x, type = "link", 
        offset = offset)
    if (family == "gaussian") 
        f$dispersion <- bglm(y ~ f$linear.predictors - 1, weights = weights, 
            start = 1, prior = "de", prior.mean = 1, prior.scale = 0, 
            verbose = FALSE)$dispersion
    f$x <- x
    f$y <- y
    f$family <- family
    f$alpha <- alpha
    f$penalty.factor <- penalty.factor
    if (alpha == 1) 
        f$prior.scale <- mean(1/(f$lambda * penalty.factor * 
            nrow(x)))
    if (alpha == 0) 
        f$prior.sd <- mean(sqrt(1/(f$lambda * penalty.factor * 
            nrow(x))))
    f$weights <- weights
    f$offset <- offset
    f$aic <- deviance(f) + 2 * f$df
    f$call <- call
    class(f) <- c("glmnet", "glmNet", "GLM")
    if (family == "cox") 
        class(f) <- c("glmnet", "glmNet", "COXPH")
    stop.time <- Sys.time()
    minutes <- round(difftime(stop.time, start.time, units = "min"), 
        3)
    if (verbose) {
        cat("\n")
        cat("Computational time:", minutes, "minutes \n")
    }
    return(f)
}