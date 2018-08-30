
###############################################################################
## Code to simulate GWAS SNP variants + covariates: Logistic error terms
###############################################################################

# 10000 * 2 repetitions run in groups of 2500 * 2 for 9 numbers of treatment 
# groups (tx) * sample size per treatment group (n_tx) combinations >> 
# 36 jobs total

# fixed:
# 1) number of chr/genes + number of SNP's per chr/gene
# 2) correlation parameters
# 3) minor allele frequencies
# 4) distributional parameters
# 5) SNP effects + locations
# 6) beta coefficients

# varying:
# 1) sample size per tx group = n_tx
# 2) number of tx groups = tx
#
# so marginal, support, total sample size changes


runID <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

library("SimCorrMix")

# folder to hold datasets with logistic errors
f2 <- "Logistic_Errors/"

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

# correlation parameters
rho_dchr <- 0.05 # between SNP's of different chr/gene
rho_scov <- 0.2 # between SNP's and covariates
rho_cov <- 0.2 # between covariates

# simulate minor allele frequencies from U(0.1, 0.5) distribution
set.seed(1)
p <- matrix(runif(chr * snp, 0.05, 0.5), 10, 5)
p <- apply(p, 2, function(x) sort(x, decreasing = TRUE))
q <- apply(p, 2, function(x) rep(1, snp) - x)

# generate marginal distributions under Hardy-Weinberg Equilibrium
marginal <- sapply(as.vector(q), function(x) c(x^2, x^2 + 2 * (1 - x) * x),
  simplify = FALSE)
# code minor allele counts as 0, 1, 2
support <- rep(list(c(0, 1, 2)), chr * snp)

# Covariates: Sex, age, ethnicity
# Sex: female = 0; male = 1
# ethnicity: Caucasian = 0; African-American = 1; Other = 2
marginal <- append(marginal, list(0.4, c(0.6, 0.9)))
support <- append(support, list(0:1, 0:2))

# age: Normal(age_mu, age_sd)
age_mu <- 21
age_sd <- 2

# error term: Logistic(0, 1)
Log <- calc_theory("Logistic", c(0, 1))

means <- c(age_mu, Log[1])
vars <- c(age_sd^2, Log[2]^2)
skews <- c(0, Log[3])
skurts <- c(0, Log[4])
fifths <- c(0, Log[5])
sixths <- c(0, Log[6])
Six <- list(NULL, 1.75)

k_cat <- snp * chr + 3
k_cont <- 2
k_mix <- 0
k_comp <- 0
k_pois <- 0
k_nb <- 0
k_total <- k_cat + k_cont + k_comp + k_pois + k_nb

# correlation matrix
rho <- matrix(rho_dchr, k_total, k_total)
colnames(rho) <- rownames(rho) <- c(paste("C1M", 1:10, sep = ""),
  paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""),
  paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""), 
  "Sex", "Race", "Tx", "Age", "E")

# correlation between SNP's of same chr/gene set based on distance

max_rho <- 0.8
rho_chr <- diag(snp)
for (i in 1:(snp - 1)) {
  for (j in (i + 1):snp) {
    if ((j - i) == 1) rho_chr[i, j] <- max_rho
    if ((j - i) == 2) rho_chr[i, j] <- max_rho - 0.1
    if ((j - i) == 3) rho_chr[i, j] <- max_rho - 0.2
    if ((j - i) == 4) rho_chr[i, j] <- max_rho - 0.3
    if ((j - i) == 5) rho_chr[i, j] <- max_rho - 0.4
    if ((j - i) == 6) rho_chr[i, j] <- max_rho - 0.5
    if ((j - i) == 7) rho_chr[i, j] <- max_rho - 0.6
    if ((j - i) == 8) rho_chr[i, j] <- max_rho - 0.65
    if ((j - i) == 9) rho_chr[i, j] <- max_rho - 0.7
    rho_chr[j, i] <- rho_chr[i, j]
  }
}

rho[1:snp, 1:snp] <- rho[(snp + 1):(2 * snp), (snp + 1):(2 * snp)] <- 
  rho[(2 * snp + 1):(3 * snp), (2 * snp + 1):(3 * snp)] <- 
  rho[(3 * snp + 1):(4 * snp), (3 * snp + 1):(4 * snp)] <- 
  rho[(4 * snp + 1):(5 * snp), (4 * snp + 1):(5 * snp)] <- rho_chr

# correlation between SNP's and covariates set at rho_scov
rho[1:(snp * chr),
  which(colnames(rho) == "Sex"):which(colnames(rho) == "Age")] <-
  matrix(rho_scov, snp * chr, 4)
rho[which(rownames(rho) == "Sex"):which(rownames(rho) == "Age"),
  1:(snp * chr)] <- matrix(rho_scov, 4, snp * chr)

# correlation between covariates set at rho_cov
rho[which(rownames(rho) == "Sex"):which(rownames(rho) == "Age"),
  which(colnames(rho) == "Sex"):which(colnames(rho) == "Age")] <-
  matrix(rho_cov, 4, 4)

# error term uncorrelated with other variables
rho[which(rownames(rho) == "E"), ] <- rep(0, ncol(rho))
rho[, which(colnames(rho) == "E")] <- rep(0, nrow(rho))

diag(rho) <- 1

Sim <- list()

n <- tx[runID] * n_tx[runID]
marginal <- append(marginal, list(seq(1/tx[runID], 1 - 1/tx[runID], 
  1/tx[runID])))
support <- append(support, list(0:(tx[runID] - 1)))

# 10,000 training sets

for (r in 1:length(rep)) {
  Com_sim <- corrvar(n, k_cat, k_cont, k_mix, k_pois, k_nb,
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six,
    marginal = marginal, support = support, rho = rho, seed = rep[r], 
    use.nearPD = FALSE, quiet = TRUE)
  for (k in 1:(chr * snp)) {
    if (sum(Com_sim$Y_cat[, k] == 2) < 5) {
      Com_sim$Y_cat[, k] <- replace(Com_sim$Y_cat[, k], 
        which(Com_sim$Y_cat[, k] == 2), 1)
    }
  }
  Sim[[r]] <- Com_sim
}
saveRDS(Sim, paste(f2, "Sim_", tx[runID], "_", n_tx[runID], "_", rep[1], 
  sep = ""))

# 10,000 validation sets

seed <- rep + 10000

for (r in 1:length(rep)) {
  Com_sim <- corrvar(n, k_cat, k_cont, k_mix, k_pois, k_nb, 
    "Polynomial", means, vars, skews, skurts, fifths, sixths, Six, 
    marginal = marginal, support = support, rho = rho, seed = seed[r], 
    use.nearPD = FALSE, quiet = TRUE)
  for (k in 1:(chr * snp)) {
    if (sum(Com_sim$Y_cat[, k] == 2) < 5) {
      Com_sim$Y_cat[, k] <- replace(Com_sim$Y_cat[, k], 
        which(Com_sim$Y_cat[, k] == 2), 1)
    }
  }
  Sim[[r]] <- Com_sim
}
saveRDS(Sim, paste(f2, "Sim_", tx[runID], "_", n_tx[runID], "_", seed[1], 
                   ".rds", sep = ""))

rm(list = ls())