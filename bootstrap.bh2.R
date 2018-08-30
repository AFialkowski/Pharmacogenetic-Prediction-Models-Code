# This function is a scaled down version of another bootstrapping function.  This was done to 
# reduce computational time.  This function works for type = "data" and bca.ci = FALSE only:
# does case resampling (select n observations with replacement, each with equal probability)

# The elements of the returned list are:
# coefs.est: observed (coefficient from original model), mean, sd, median, 
#            b.l (lower bound of bootstrap percentile CI), b.h (upper bound of bootstrap percentile CI), pvalue,
#            lower_norm (lower bound of normal CI), upper_norm (upper bound of normal CI)
# coefs: the bootstrapped coefficient estimates, with the 1st column equal to the coefficients from the original model
# coefs.p: proportion of times each coefficient is non-zero in the bootstrap samples (for glmNet or bmlasso models).
#
# Reference: Fox, J. and Weisberg, S. (2011). An R Companion to Applied Regression. 
#            Sage, Thousand Oaks, CA, second edition.  
#            https://socserv.socsci.mcmaster.ca/jfox/Books/Companion/appendix/Appendix-Bootstrapping.pdf.

bootstrap.bh2 <- function(object, nbs = 100, verbose = FALSE, alpha = 0.05, 
                          boot.type = c("data", "residual", "smooth"), 
                          bca.ci = FALSE, seed = 1234, scale.y = FALSE) 
{
  start.time <- Sys.time()
  set.seed(seed)
  require(glmnet)
  x.obj <- object$x
  y.obj <- object$y
  n <- NROW(y.obj)
  out <- list()
  measures <- NULL
  coefs <- object$coefficients
  lp <- NULL
  if (boot.type == "data") {
    for (k in 1:nbs) {
      obs <- sample(1:n, size = n, replace = TRUE)
      fit <- update(object, x = x.obj[obs, ], y = y.obj[obs], 
        weights = object$weights[obs], offset = object$offset[obs], 
        alpha = object$alpha, lambda = object$lambda, verbose = FALSE, 
        scale.y = scale.y)
      coefs <- cbind(coefs, fit$coefficients)
      if (verbose) {
        pre <- rep("\b", nbs)
        cat(pre, k, "/", nbs, sep = "")
        flush.console()
      }
    }
  }
  Mean <- apply(coefs, 1, mean)
  Median <- apply(coefs, 1, median)
  Sd <- apply(coefs, 1, sd)
  b.l <- apply(coefs, 1, quantile, probs = 0.5 * alpha)
  b.h <- apply(coefs, 1, quantile, probs = 1 - 0.5 * alpha)
  tvalue <- Mean/(Sd + 1e-10)
  pvalue <- 2 * pnorm(-abs(tvalue))
  coefs.est <- as.data.frame(cbind(coefs[, 1], Mean, Sd, Median, b.l, b.h, 
                                   pvalue))
  colnames(coefs.est) <- c("original", "mean", "sd", "median", 
                           paste(0.5 * alpha * 100, "%", sep = ""), 
                           paste(100 - 0.5 * alpha * 100, "%", sep = ""), 
                           "pvalue")
  coefs.est$lower_norm <- 2 * coefs.est$original - coefs.est$mean - 
    qnorm(1 - 0.5 * alpha) * coefs.est$sd
  coefs.est$upper_norm <- 2 * coefs.est$original - coefs.est$mean + 
    qnorm(1 - 0.5 * alpha) * coefs.est$sd
  colnames(coefs)[1] <- "original"
  for(i in 2:ncol(coefs)) {
    colnames(coefs)[i] <- paste("boot", (i - 1), sep = "")
  }
  out$coefs.est <- coefs.est
  out$coefs <- coefs
  nonzero <- ifelse(coefs == 0, 0, 1)
  nonzero <- apply(nonzero, 1, mean)
  out$coefs.p <- nonzero
  stop.time <- Sys.time()
  out$Time <- round(difftime(stop.time, start.time, units = "min"), 
                3)
  if (verbose) {
    cat("\n")
    cat("Bootstrap time:", Time, "minutes \n")
  }
  out
}
