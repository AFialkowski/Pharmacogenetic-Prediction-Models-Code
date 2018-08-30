
runID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
alpha = c(1, 0.25, 0.5, 0.75)

library("BhGLM")
library("glmnet")

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

seed <- rep

# find best models in training sets

Blasso1 <- list()
Blasso2 <- list()
Sum_Blasso1 <- list()

start <- Sys.time()
for (r in 1:length(rep)) {
  Xy <- Vars[[r]]
  set.seed(seed[r])
  Blasso1[[r]] <- bglm(y ~ ., data = Xy, prior = "de", 
   prior.scale = EN1_cv[[r]]$prior.scale)
}
stop <- Sys.time()
saveRDS(Blasso1, paste("Blasso1_", tx[runID], "_", n_tx[runID], "_", 
  rep[1], ".rds", sep = ""))
Time <- round(difftime(stop, start, units = "secs"), 8)
saveRDS(Time, paste("Time_Blasso1_", tx[runID], "_", n_tx[runID], "_", 
  rep[1], ".rds", sep = ""))

# use best models to predict responses in training sets

for (r in 1:length(rep)) {
  Xy2 <- Vars[[r]]
  set.seed(seed[r])
  Blasso2[[r]] <- predict.bh(Blasso1[[r]], new.x = Xy2[, -1], new.y = Xy2[, 1])
}
saveRDS(Blasso2, paste("train_Blasso2_", tx[runID], "_", n_tx[runID], "_", 
                       rep[1], ".rds", sep = ""))

# summarize models in training sets and use best models to predict responses 
# in validation sets

Blasso2 <- list() # predictions

for (r in 1:length(rep)) {
  Xy2 <- Vars2[[r]]
  Sum_Blasso1[[r]] <- summary.bh(Blasso1[[r]])
  set.seed(seed[r])
  Blasso2[[r]] <- predict.bh(Blasso1[[r]], new.x = Xy2[, -1], new.y = Xy2[, 1])
}
saveRDS(Sum_Blasso1, paste("Sum_Blasso1_", tx[runID], "_", n_tx[runID], "_", rep[1], 
                           ".rds", sep = ""))
saveRDS(Blasso2, paste("Blasso2_", tx[runID], "_", n_tx[runID], "_", 
                       rep[1], ".rds", sep = ""))


rm(list = ls())