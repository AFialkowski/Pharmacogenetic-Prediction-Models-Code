
runID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
alpha = c(1, 0.25, 0.5, 0.75)

library("BhGLM")
library("glmnet")
source("Main/bootstrap.bh2.R")
source("Main/glmNet2.R")

# treatment group:
# number of treatment groups
tx <- c(rep(2, 12), rep(3, 12), rep(5, 12))
# sample size per tx group
n_tx <- c(rep(100, 4), rep(250, 4), rep(500, 4), 
          rep(100, 4), rep(250, 4), rep(500, 4), 
          rep(100, 4), rep(250, 4), rep(500, 4))

if (runID %in% seq(1, 33, 4)) rep <- 1:2500
if (runID %in% seq(2, 34, 4)) rep <- 2501:5000
if (runID %in% seq(3, 35, 4)) rep <- 5001:7500
if (runID %in% seq(4, 36, 4)) rep <- 7501:10000

# training set
Vars <- readRDS(paste("Vars_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  ".rds", sep = ""))
# testing set
Vars2 <- readRDS(paste("Vars_", tx[runID], "_", n_tx[runID], "_", 
  rep[1] + 10000, ".rds", sep = ""))

EN1_cv <- list() # CV model with in-sample statistics
EN1_boot <- list() # bootstrap
EN2 <- list() # predictions

if (n_tx[runID] == 100) k_folds <- 5 else k_folds <- 10

seed <- rep
ncv <- 10
nbs <- 1000

# find best models in training sets

start <- Sys.time()
for (r in 1:length(rep)) {
  Xy <- Vars[[r]]
  set.seed(seed[r])
  EN1_cv[[r]] <- glmNet2(x = Xy[, -1], y = Xy[, 1], alpha = alpha[4], 
    nfolds = k_folds, ncv = ncv, verbose = FALSE, scale.y = FALSE)
}
stop <- Sys.time()
saveRDS(EN1_cv, paste("EN1_cv_", tx[runID], "_", n_tx[runID], "_", 
  "a75_", rep[1], ".rds", sep = ""))
Time <- round(difftime(stop, start, units = "secs"), 8)
saveRDS(Time, paste("Time_EN1_cv_", tx[runID], "_", n_tx[runID], "_", 
  "a75_", rep[1], ".rds", sep = ""))

# use best models to predict responses in training sets

for (r in 1:length(rep)) {
  Xy2 <- Vars[[r]]
  set.seed(seed[r])
  EN2[[r]] <- predict.bh(EN1_cv[[r]], new.x = Xy2[, -1], new.y = Xy2[, 1])
}
saveRDS(EN2, paste("train_EN2_", tx[runID], "_", n_tx[runID], "_", 
  "a75_", rep[1], ".rds", sep = ""))

# bootstrap models in training sets and use best models to predict responses 
# in validation sets

EN2 <- list() # predictions

for (r in 1:length(rep)) {
  Xy2 <- Vars2[[r]]
  EN1_boot[[r]] <- bootstrap.bh2(EN1_cv[[r]], nbs = nbs, verbose = FALSE, 
    alpha = 0.05, boot.type = "data", bca.ci = FALSE, seed = seed[r], 
    scale.y = FALSE)
  set.seed(seed[r])
  EN2[[r]] <- predict.bh(EN1_cv[[r]], new.x = Xy2[, -1], new.y = Xy2[, 1])
}
saveRDS(EN1_boot, paste("EN1_boot_", tx[runID], "_", n_tx[runID], "_", 
  "a75_", rep[1], ".rds", sep = ""))
saveRDS(EN2, paste("EN2_", tx[runID], "_", n_tx[runID], "_", 
  "a75_", rep[1], ".rds", sep = ""))

rm(list = ls())