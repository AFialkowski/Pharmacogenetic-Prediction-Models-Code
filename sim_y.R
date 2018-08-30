
###############################################################################
## Code to simulate outcomes
###############################################################################

# 10000 * 2 repetitions run in groups of 2500 * 2 for 9 numbers of treatment 
# groups (tx) * sample size per treatment group (n_tx) combinations >> 
# 36 jobs total

runID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library("SimCorrMix")
f11 <- "Beta1/" # small tx effects
f12 <- "Beta2/" # medium tx effects
f13 <- "Beta3/" # large tx effects

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

# number of chr/genes
chr <- 5
# number of SNP's per chr/gene
snp <- 10

# SNP effects
qtl <- 1:10 # 10 non-zero SNP effects in C1
neg_beta <- c(3, 7) # which SNP effects are negative
# SNP effects
beta_snps <- c(0.7, 0.8, 0.9, 1, 1.25, 1, 0.9, 0.8, 0.7, 0.6)
beta_snps[neg_beta] <- -beta_snps[neg_beta]

# covariate effects
beta_sex <- 0.5
beta_age <- 0.5
beta_race <- matrix(c(0.5, -0.5), 2, 1)

# training sets

Sim <- readRDS(paste("Sim_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  sep = ""))

###############################################################################
## Small treatment effects
###############################################################################

beta_tx <- list(NULL, matrix(0.5, 1, 1), matrix(c(0.5, -0.25), 2, 1), NULL, 
  matrix(c(0.5, -0.5, 0.75, -0.75), 4, 1))

Vars <- list()

for (r in 1:length(rep)) {
  X <- cbind(Sim[[r]]$Y_cat, Sim[[r]]$Y_cont[, 1])
  colnames(X) <- c(paste("C1M", 1:10, sep = ""), 
    paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""), 
    paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""), 
    "Sex", "Race", "Tx", "Age")
  y <- X[, qtl] %*% matrix(beta_snps, length(beta_snps), 1) + 
    Sim[[r]]$Y_cont[, 2]
  X <- as.data.frame(X)
  X$Race_1 <- ifelse(X$Race == 1, 1, 0)
  X$Race_2 <- ifelse(X$Race == 2, 1, 0)
  X <- X[, -which(colnames(X) == "Race")]
  if (tx[runID] == 3) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  if (tx[runID] == 5) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X$Tx_3 <- ifelse(X$Tx == 3, 1, 0)
    X$Tx_4 <- ifelse(X$Tx == 4, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  X <- as.matrix(X)
  y <- y + X[, "Age"] * beta_age + X[, "Sex"] * beta_sex + 
    X[, grep("Race", colnames(X))] %*% beta_race + 
    X[, grep("Tx", colnames(X))] %*% beta_tx[[tx[runID]]]
  Xy <- as.data.frame(cbind(y, X))
  colnames(Xy)[1] <- "y"
  Vars[[r]] <- Xy
}
saveRDS(Vars, paste(f11, "Vars_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  ".rds", sep = ""))

###############################################################################
## Medium treatment effects
###############################################################################

beta_tx <- list(NULL, matrix(1, 1, 1), matrix(c(1, -0.5), 2, 1), NULL, 
  matrix(c(0.75, -0.75, 1.25, -1.25), 4, 1))

Vars <- list()

for (r in 1:length(rep)) {
  X <- cbind(Sim[[r]]$Y_cat, Sim[[r]]$Y_cont[, 1])
  colnames(X) <- c(paste("C1M", 1:10, sep = ""), 
    paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""), 
    paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""), 
    "Sex", "Race", "Tx", "Age")
  y <- X[, qtl] %*% matrix(beta_snps, length(beta_snps), 1) + 
    Sim[[r]]$Y_cont[, 2]
  X <- as.data.frame(X)
  X$Race_1 <- ifelse(X$Race == 1, 1, 0)
  X$Race_2 <- ifelse(X$Race == 2, 1, 0)
  X <- X[, -which(colnames(X) == "Race")]
  if (tx[runID] == 3) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  if (tx[runID] == 5) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X$Tx_3 <- ifelse(X$Tx == 3, 1, 0)
    X$Tx_4 <- ifelse(X$Tx == 4, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  X <- as.matrix(X)
  y <- y + X[, "Age"] * beta_age + X[, "Sex"] * beta_sex + 
    X[, grep("Race", colnames(X))] %*% beta_race + 
    X[, grep("Tx", colnames(X))] %*% beta_tx[[tx[runID]]]
  Xy <- as.data.frame(cbind(y, X))
  colnames(Xy)[1] <- "y"
  Vars[[r]] <- Xy
}
saveRDS(Vars, paste(f12, "Vars_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  ".rds", sep = ""))

###############################################################################
## Large treatment effects
###############################################################################

beta_tx <- list(NULL, matrix(1.5, 1, 1), matrix(c(1.5, -0.75), 2, 1), NULL, 
  matrix(c(1, -1, 1.5, -1.5), 4, 1))

Vars <- list()

for (r in 1:length(rep)) {
  X <- cbind(Sim[[r]]$Y_cat, Sim[[r]]$Y_cont[, 1])
  colnames(X) <- c(paste("C1M", 1:10, sep = ""), 
    paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""), 
    paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""), 
    "Sex", "Race", "Tx", "Age")
  y <- X[, qtl] %*% matrix(beta_snps, length(beta_snps), 1) + 
    Sim[[r]]$Y_cont[, 2]
  X <- as.data.frame(X)
  X$Race_1 <- ifelse(X$Race == 1, 1, 0)
  X$Race_2 <- ifelse(X$Race == 2, 1, 0)
  X <- X[, -which(colnames(X) == "Race")]
  if (tx[runID] == 3) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  if (tx[runID] == 5) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X$Tx_3 <- ifelse(X$Tx == 3, 1, 0)
    X$Tx_4 <- ifelse(X$Tx == 4, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  X <- as.matrix(X)
  y <- y + X[, "Age"] * beta_age + X[, "Sex"] * beta_sex + 
    X[, grep("Race", colnames(X))] %*% beta_race + 
    X[, grep("Tx", colnames(X))] %*% beta_tx[[tx[runID]]]
  Xy <- as.data.frame(cbind(y, X))
  colnames(Xy)[1] <- "y"
  Vars[[r]] <- Xy
}
saveRDS(Vars, paste(f13, "Vars_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  ".rds", sep = ""))

# validation sets

rep <- rep + 10000

Sim <- readRDS(paste("Sim_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  ".rds", sep = ""))

###############################################################################
## Small treatment effects
###############################################################################

beta_tx <- list(NULL, matrix(0.5, 1, 1), matrix(c(0.5, -0.25), 2, 1), NULL, 
  matrix(c(0.5, -0.5, 0.75, -0.75), 4, 1))

Vars <- list()

for (r in 1:length(rep)) {
  X <- cbind(Sim[[r]]$Y_cat, Sim[[r]]$Y_cont[, 1])
  colnames(X) <- c(paste("C1M", 1:10, sep = ""), 
    paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""), 
    paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""), 
    "Sex", "Race", "Tx", "Age")
  y <- X[, qtl] %*% matrix(beta_snps, length(beta_snps), 1) + 
    Sim[[r]]$Y_cont[, 2]
  X <- as.data.frame(X)
  X$Race_1 <- ifelse(X$Race == 1, 1, 0)
  X$Race_2 <- ifelse(X$Race == 2, 1, 0)
  X <- X[, -which(colnames(X) == "Race")]
  if (tx[runID] == 3) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  if (tx[runID] == 5) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X$Tx_3 <- ifelse(X$Tx == 3, 1, 0)
    X$Tx_4 <- ifelse(X$Tx == 4, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  X <- as.matrix(X)
  y <- y + X[, "Age"] * beta_age + X[, "Sex"] * beta_sex + 
    X[, grep("Race", colnames(X))] %*% beta_race + 
    X[, grep("Tx", colnames(X))] %*% beta_tx[[tx[runID]]]
  Xy <- as.data.frame(cbind(y, X))
  colnames(Xy)[1] <- "y"
  Vars[[r]] <- Xy
}
saveRDS(Vars, paste(f11, "Vars_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  ".rds", sep = ""))

###############################################################################
## Medium treatment effects
###############################################################################

beta_tx <- list(NULL, matrix(1, 1, 1), matrix(c(1, -0.5), 2, 1), NULL, 
  matrix(c(0.75, -0.75, 1.25, -1.25), 4, 1))

Vars <- list()

for (r in 1:length(rep)) {
  X <- cbind(Sim[[r]]$Y_cat, Sim[[r]]$Y_cont[, 1])
  colnames(X) <- c(paste("C1M", 1:10, sep = ""), 
    paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""), 
    paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""), 
    "Sex", "Race", "Tx", "Age")
  y <- X[, qtl] %*% matrix(beta_snps, length(beta_snps), 1) + 
    Sim[[r]]$Y_cont[, 2]
  X <- as.data.frame(X)
  X$Race_1 <- ifelse(X$Race == 1, 1, 0)
  X$Race_2 <- ifelse(X$Race == 2, 1, 0)
  X <- X[, -which(colnames(X) == "Race")]
  if (tx[runID] == 3) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  if (tx[runID] == 5) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X$Tx_3 <- ifelse(X$Tx == 3, 1, 0)
    X$Tx_4 <- ifelse(X$Tx == 4, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  X <- as.matrix(X)
  y <- y + X[, "Age"] * beta_age + X[, "Sex"] * beta_sex + 
    X[, grep("Race", colnames(X))] %*% beta_race + 
    X[, grep("Tx", colnames(X))] %*% beta_tx[[tx[runID]]]
  Xy <- as.data.frame(cbind(y, X))
  colnames(Xy)[1] <- "y"
  Vars[[r]] <- Xy
}
saveRDS(Vars, paste(f12, "Vars_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  ".rds", sep = ""))

###############################################################################
## Large treatment effects
###############################################################################

beta_tx <- list(NULL, matrix(1.5, 1, 1), matrix(c(1.5, -0.75), 2, 1), NULL, 
  matrix(c(1, -1, 1.5, -1.5), 4, 1))

Vars <- list()

for (r in 1:length(rep)) {
  X <- cbind(Sim[[r]]$Y_cat, Sim[[r]]$Y_cont[, 1])
  colnames(X) <- c(paste("C1M", 1:10, sep = ""), 
    paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""), 
    paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""), 
    "Sex", "Race", "Tx", "Age")
  y <- X[, qtl] %*% matrix(beta_snps, length(beta_snps), 1) + 
    Sim[[r]]$Y_cont[, 2]
  X <- as.data.frame(X)
  X$Race_1 <- ifelse(X$Race == 1, 1, 0)
  X$Race_2 <- ifelse(X$Race == 2, 1, 0)
  X <- X[, -which(colnames(X) == "Race")]
  if (tx[runID] == 3) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  if (tx[runID] == 5) {
    X$Tx_1 <- ifelse(X$Tx == 1, 1, 0)
    X$Tx_2 <- ifelse(X$Tx == 2, 1, 0)
    X$Tx_3 <- ifelse(X$Tx == 3, 1, 0)
    X$Tx_4 <- ifelse(X$Tx == 4, 1, 0)
    X <- X[, -which(colnames(X) == "Tx")]
  }
  X <- as.matrix(X)
  y <- y + X[, "Age"] * beta_age + X[, "Sex"] * beta_sex + 
    X[, grep("Race", colnames(X))] %*% beta_race + 
    X[, grep("Tx", colnames(X))] %*% beta_tx[[tx[runID]]]
  Xy <- as.data.frame(cbind(y, X))
  colnames(Xy)[1] <- "y"
  Vars[[r]] <- Xy
}
saveRDS(Vars, paste(f13, "Vars_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  ".rds", sep = ""))


rm(list = ls())