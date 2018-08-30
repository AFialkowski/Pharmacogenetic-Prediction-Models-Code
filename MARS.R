
# Due to the long computational times, MARS models were run in 180 jobs with
# 500 datasets each.

runID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library("earth")

# treatment group:
# number of treatment groups
tx <- c(rep(2, 60), rep(3, 60), rep(5, 60))
# sample size per tx group
n_tx <- c(rep(100, 20), rep(250, 20), rep(500, 20), 
          rep(100, 20), rep(250, 20), rep(500, 20), 
          rep(100, 20), rep(250, 20), rep(500, 20))

if (runID %in% seq(1, 161, 20)) rep <- 1:500
if (runID %in% seq(2, 162, 20)) rep <- 501:1000
if (runID %in% seq(3, 163, 20)) rep <- 1001:1500
if (runID %in% seq(4, 164, 20)) rep <- 1501:2000
if (runID %in% seq(5, 165, 20)) rep <- 2001:2500
if (runID %in% seq(6, 166, 20)) rep <- 2501:3000
if (runID %in% seq(7, 167, 20)) rep <- 3001:3500
if (runID %in% seq(8, 168, 20)) rep <- 3501:4000
if (runID %in% seq(9, 169, 20)) rep <- 4001:4500
if (runID %in% seq(10, 170, 20)) rep <- 4501:5000
if (runID %in% seq(11, 171, 20)) rep <- 5001:5500
if (runID %in% seq(12, 172, 20)) rep <- 5501:6000
if (runID %in% seq(13, 173, 20)) rep <- 6001:6500
if (runID %in% seq(14, 174, 20)) rep <- 6501:7000
if (runID %in% seq(15, 175, 20)) rep <- 7001:7500
if (runID %in% seq(16, 176, 20)) rep <- 7501:8000
if (runID %in% seq(17, 177, 20)) rep <- 8001:8500
if (runID %in% seq(18, 178, 20)) rep <- 8501:9000
if (runID %in% seq(19, 179, 20)) rep <- 9001:9500
if (runID %in% seq(20, 180, 20)) rep <- 9501:10000

if (min(rep) %in% 1:2001) rep2 <- 1:2500
if (min(rep) %in% 2501:4501) rep2 <- 2501:5000
if (min(rep) %in% 5001:7001) rep2 <- 5001:7500
if (min(rep) %in% 7501:9501) rep2 <- 7501:10000

Vars <- readRDS(paste("Vars_", tx[runID], "_", n_tx[runID], "_", rep2[1], 
                      ".rds", sep = ""))
Vars2 <- readRDS(paste("Vars_", tx[runID], "_", n_tx[runID], "_", 
                       rep2[1] + 10000, ".rds", sep = ""))
if (min(rep) %in% c(1, 2501, 5001, 7501)) Vars <- Vars[1:500]
if (min(rep) %in% c(501, 3001, 5501, 8001)) Vars <- Vars[501:1000]
if (min(rep) %in% c(1001, 3501, 6001, 8501)) Vars <- Vars[1001:1500]
if (min(rep) %in% c(1501, 4001, 6501, 9001)) Vars <- Vars[1501:2000]
if (min(rep) %in% c(2001, 4501, 7001, 9501)) Vars <- Vars[2001:2500]
if (min(rep) %in% c(1, 2501, 5001, 7501)) Vars2 <- Vars2[1:500]
if (min(rep) %in% c(501, 3001, 5501, 8001)) Vars2 <- Vars2[501:1000]
if (min(rep) %in% c(1001, 3501, 6001, 8501)) Vars2 <- Vars2[1001:1500]
if (min(rep) %in% c(1501, 4001, 6501, 9001)) Vars2 <- Vars2[1501:2000]
if (min(rep) %in% c(2001, 4501, 7001, 9501)) Vars2 <- Vars2[2001:2500]

if (n_tx[runID] == 100) k_folds <- 5 else k_folds <- 10

seed <- rep
ncv <- 5

# find best models in training sets

Earth1_cv <- list() 
Earth2 <- list()
data0 <- list()
data2 <- list()

start <- Sys.time()
for (r in 1:length(rep)) {
  earth1_cv <- NULL
  Xy <- as.data.frame(Vars[[r]])
  for (i in 2:51) {
    Xy[, i] <- as.factor(Xy[, i])
  }
  Xy2 <- as.data.frame(Vars2[[r]])
  for (i in 2:51) {
    Xy2[, i] <- as.factor(Xy2[, i])
  }
  y <- Xy[, 1]
  Xy <- model.matrix(y ~ . - 1, data = Xy)
  Xy2 <- model.matrix(y ~ . - 1, data = Xy2)
  both <- intersect(colnames(Xy), colnames(Xy2))
  Xy <- Xy[, which(colnames(Xy) %in% both)]
  data0[[r]] <- Xy
  data2[[r]] <- Xy2[, which(colnames(Xy2) %in% both)]
  set.seed(seed[r])
  earth1_cv <- tryCatch({
    earth(x = Xy, y = y, keepxy = TRUE, degree = 1, ncross = ncv, 
      nfold = k_folds, Scale.y = FALSE, pmethod = "cv")
  }, error=function(e){})
  if (length(earth1_cv) == 0) {
    Earth1_cv <- append(Earth1_cv, list(NULL))
  } else {
    Earth1_cv[[r]] <- earth1_cv
  }
}
stop <- Sys.time()
Time <- round(difftime(stop, start, units = "secs"), 8)

# use best models to predict responses in training sets

for (r in 1:length(rep)) {
  if (length(Earth1_cv[[r]]) == 0) {
    Earth2 <- append(Earth2, list(NULL))
  } else {
    set.seed(seed[r])
    Earth2[[r]] <- predict(Earth1_cv[[r]], newdata = data0[[r]])
  }
}

saveRDS(Earth2, paste("train_Earth2_", tx[runID], "_", n_tx[runID], "_", 
                      rep[1], ".rds", sep = ""))

# bootstrap models in training sets and use best models to predict responses 
# in validation sets

Earth2 <- list() 

for (r in 1:length(rep)) {
  if (length(Earth1_cv[[r]]) == 0) {
    Earth2 <- append(Earth2, list(NULL))
  } else {
    set.seed(seed[r])
    Earth2[[r]] <- predict(Earth1_cv[[r]], newdata = data2[[r]])
  }
}

Earth1_cv <- lapply(Earth1_cv, function(x) list(coefficients = x$coefficients))
saveRDS(Earth1_cv, paste("Earth1_cv_", tx[runID], "_", n_tx[runID], "_", 
  rep2[1], "_", rep[1], ".rds", sep = ""))
saveRDS(Earth2, paste("Earth2_", tx[runID], "_", n_tx[runID], "_", 
  rep2[1], "_", rep[1], ".rds", sep = ""))
saveRDS(Time, paste("Time_Earth1_cv_", tx[runID], "_", n_tx[runID], "_", 
  rep2[1], "_", rep[1], ".rds", sep = ""))

rm(list = ls())