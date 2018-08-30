###############################################################################
## Code to summarize model results with medium tx effects: Logistic error terms
###############################################################################

library("psych")
library("dplyr")

f2 <- "Logistic_Errors/"

# code to combine earth files
ind1 <- seq(1, 2001, 500)
ind2501 <- ind1 + 2500
ind5001 <- ind1 + 5000
ind7501 <- ind1 + 7500

ind <- list(ind1, ind2501, ind5001, ind7501)

# treatment group:
# number of treatment groups
tx <- c(2, 3, 5)
# sample size per tx group
n_tx <- c(100, 250, 500)

for (i in 1:length(tx)) {
  for (j in 1:length(n_tx)) {
    for (k in 1:length(ind)) {
      mars <- mars2 <- mars3 <- list()
      Time <- 0
      for (l in 1:length(ind[[k]])) {
        mars <- append(mars, readRDS(paste("Earth1_cv_", tx[i], "_", n_tx[j], 
          "_", ind[[k]][1], "_", ind[[k]][l], ".rds", sep = "")))
        mars2 <- append(mars2, readRDS(paste("Earth2_", tx[i], "_", n_tx[j], 
          "_", ind[[k]][1], "_", ind[[k]][l], ".rds", sep = "")))
        mars3 <- append(mars3, readRDS(paste("train_Earth2_", tx[i], "_", 
          n_tx[j], "_", ind[[k]][1], "_", ind[[k]][l], ".rds", sep = "")))
        Time <- sum(Time, as.numeric(readRDS(paste("Time_Earth1_cv_", tx[i], 
          "_", n_tx[j], "_", ind[[k]][1], "_", ind[[k]][l], ".rds", 
          sep = ""))))
      }
      saveRDS(Time, paste("Time_Earth1_cv_", tx[i], "_", n_tx[j], "_", 
        ind[[k]][1], ".rds", sep = ""))
      saveRDS(mars, paste("Earth1_cv_", tx[i], "_", n_tx[j], "_", ind[[k]][1], 
        ".rds", sep = ""))
      saveRDS(mars2, paste("Earth2_", tx[i], "_", n_tx[j], "_", ind[[k]][1], 
        ".rds", sep = ""))
      saveRDS(mars2, paste("train_Earth2_", tx[i], "_", n_tx[j], "_", 
        ind[[k]][1], ".rds", sep = ""))
    }
  }
}


# 1st element of file indices
rep <- c(1, 2501, 5001, 7501)

# treatment group:
# number of treatment groups
tx <- c(rep(2, 3), rep(3, 3), rep(5, 3))
# sample size per tx group
n_tx <- rep(c(100, 250, 500), 3)

# computation time
for (i in 1:length(tx)) {
  Time0 <- Time1 <- Time2 <- Time3 <- Time4 <- Time5 <- 0
  for (r in 1:length(rep)) {
    Time0 <- sum(Time0, as.numeric(readRDS(paste("Time_EN1_cv_", tx[i], "_", 
      n_tx[i], "_", rep[r], ".rds", sep = ""))))
    Time1 <- sum(Time1, as.numeric(readRDS(paste("Time_EN1_cv_", tx[i], "_", 
      n_tx[i], "_a25_", rep[r], ".rds", sep = ""))))
    Time2 <- sum(Time2, as.numeric(readRDS(paste("Time_EN1_cv_", tx[i], "_", 
      n_tx[i], "_a50_", rep[r], ".rds", sep = ""))))
    Time3 <- sum(Time3, as.numeric(readRDS(paste("Time_EN1_cv_", tx[i], "_", 
      n_tx[i], "_a75_", rep[r], ".rds", sep = ""))))
    Time4 <- sum(Time4, as.numeric(readRDS(paste("Time_Blasso1_", tx[i], "_", 
      n_tx[i], "_", rep[r], ".rds", sep = ""))))
    Time5 <- sum(Time5, as.numeric(readRDS(paste("Time_Earth1_cv_", tx[i], "_", 
      n_tx[i], "_", rep[r], ".rds", sep = ""))))
  }
  saveRDS(Time0, paste(f2, "Time_EN1_cv_", tx[i], "_", n_tx[i], "_betas2.rds", 
    sep = ""))
  saveRDS(Time1, paste(f2, "Time_EN1_cv_", tx[i], "_a25_", n_tx[i], 
    "_betas2.rds", sep = ""))
  saveRDS(Time2, paste(f2, "Time_EN1_cv_", tx[i], "_a50_", n_tx[i], 
    "_betas2.rds", sep = ""))
  saveRDS(Time2, paste(f2, "Time_EN1_cv_", tx[i], "_a75_", n_tx[i], 
    "_betas2.rds", sep = ""))
  saveRDS(Time4, paste(f2, "Time_Blasso1_", tx[i], "_", n_tx[i], 
    "_betas2.rds", sep = ""))
  saveRDS(Time5, paste(f2, "Time_Earth1_cv_", tx[i], "_", n_tx[i], 
    "_betas2.rds", sep = ""))
}
    
# summarize outcomes
n_tx <- c(100, 250, 500)
Sum_y <- matrix(0, 9, 9)
Sum_y[, 1] <- c(2, 2, 2, 3, 3, 3, 5, 5, 5)
Sum_y[, 2] <- rep(n_tx, 3)
for (j in 1:length(n_tx)) {
  y <- NULL
  for (r in 1:length(rep)) {
    Vars <- readRDS(paste("Vars_", 2, "_", n_tx[j], "_", rep[r], ".rds", sep = ""))
    for (k in 1:length(Vars)) {
      y <- cbind(y, Vars[[k]][, 1])
    }
  }
  y <- as.data.frame(y)
  colnames(y) <- paste("y", 1:10000, sep = "")
  y_sum <- describe(y, quant = c(0.25, 0.75))
  write.table(y_sum, paste(f2, "Sum_y_", 2, "_", n_tx[j], "_betas2.txt", sep = ""), row.names = T, col.names = T, 
              sep = "\t", quote = F)
  Sum_y[j, 3:9] <- colMeans(y_sum)[c("n", "mean", "sd", "min", "Q0.25", "Q0.75", "max")]
  cat(2, "tx groups,", n_tx[j], "per tx group")
}
for (j in 1:length(n_tx)) {
  y <- NULL
  for (r in 1:length(rep)) {
    Vars <- readRDS(paste("Vars_", 3, "_", n_tx[j], "_", rep[r], ".rds", sep = ""))
    for (k in 1:length(Vars)) {
      y <- cbind(y, Vars[[k]][, 1])
    }
  }
  y <- as.data.frame(y)
  colnames(y) <- paste("y", 1:10000, sep = "")
  y_sum <- describe(y, quant = c(0.25, 0.75))
  write.table(y_sum, paste(f2, "Sum_y_", 3, "_", n_tx[j], "_betas2.txt", sep = ""), row.names = T, col.names = T, 
              sep = "\t", quote = F)
  Sum_y[j + 3, 3:9] <- colMeans(y_sum)[c("n", "mean", "sd", "min", "Q0.25", "Q0.75", "max")]
  cat(3, "tx groups,", n_tx[j], "per tx group")
}
for (j in 1:length(n_tx)) {
  y <- NULL
  for (r in 1:length(rep)) {
    Vars <- readRDS(paste("Vars_", 5, "_", n_tx[j], "_", rep[r], ".rds", sep = ""))
    for (k in 1:length(Vars)) {
      y <- cbind(y, Vars[[k]][, 1])
    }
  }
  y <- as.data.frame(y)
  colnames(y) <- paste("y", 1:10000, sep = "")
  y_sum <- describe(y, quant = c(0.25, 0.75))
  write.table(y_sum, paste(f2, "Sum_y_", 5, "_", n_tx[j], "_betas2.txt", sep = ""), row.names = T, col.names = T, 
              sep = "\t", quote = F)
  Sum_y[j + 6, 3:9] <- colMeans(y_sum)[c("n", "mean", "sd", "min", "Q0.25", "Q0.75", "max")]
  cat(5, "tx groups,", n_tx[j], "per tx group")
}
Sum_y <- as.data.frame(Sum_y)
colnames(Sum_y) <- c("Tx", "N_Tx", "N", "Mean", "SD", "Min", "Q0.25", "Q0.75", "Max")
write.table(Sum_y, paste(f2, "Sum_y_betas2.txt", sep = ""), row.names = T, col.names = T, 
            sep = "\t", quote = F)


# Find errors for estimated treatment effects

terr <- NULL
rep <- c(1, 2501, 5001, 7501)

tx <- rep(2, 3)
beta_tx <- 1

n_tx <- c(100, 250, 500)
for (i in 1:length(n_tx)) {
  # Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR <- sapply(pred, function(x) (x$coefs.est["Tx", 1] - beta_tx[1]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(0, length(pred)), TERR)
  terr <- rbind(terr, terr2)
  # EN: alpha = 0.25
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_a25_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR <- sapply(pred, function(x) (x$coefs.est["Tx", 1] - beta_tx[1]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(1, length(pred)), TERR)
  terr <- rbind(terr, terr2)
  # EN: alpha = 0.50
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_a50_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR <- sapply(pred, function(x) (x$coefs.est["Tx", 1] - beta_tx[1]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(2, length(pred)), TERR)
  terr <- rbind(terr, terr2)
  # EN: alpha = 0.75
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_a75_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR <- sapply(pred, function(x) (x$coefs.est["Tx", 1] - beta_tx[1]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(3, length(pred)), TERR)
  terr <- rbind(terr, terr2)
    # Bayesian Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("Sum_Blasso1_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR <- sapply(pred, function(x) (x["Tx", 1] - beta_tx[1]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(4, length(pred)), TERR)
  terr <- rbind(terr, terr2)
  # MARS
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("Earth1_cv_", tx[i], "_", n_tx[i], "_", 
      rep[r], ".rds", sep = "")))
  }
  TERR <- sapply(pred, function(x) 
    if (length(x) == 0) 
      NA 
    else if ("Tx" %in% rownames(x$coefficients)) 
      (x$coefficients["Tx", 1] - beta_tx[1]) else (0 - beta_tx[1]))
  TERR2 <- sapply(TERR, function(x) if (is.na(x)) NA else x)
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(5, length(pred)), TERR2)
  terr <- rbind(terr, terr2)
  print(i)
}

colnames(terr) <- c("Tx", "N_Tx", "Model", "TERR")
terr <- as.data.frame(terr)
write.table(terr, paste(f2, "tx2_terr_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

terr <- NULL
tx <- rep(3, 3)
beta_tx <- c(1, -0.5)

n_tx <- c(100, 250, 500)
for (i in 1:length(n_tx)) {
  # Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x$coefs.est["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x$coefs.est["Tx_2", 1] - beta_tx[2]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(0, length(pred)), TERR1, TERR2)
  terr <- rbind(terr, terr2)
  # EN: alpha = 0.25
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_a25_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x$coefs.est["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x$coefs.est["Tx_2", 1] - beta_tx[2]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(1, length(pred)), TERR1, TERR2)
  terr <- rbind(terr, terr2)
  # EN: alpha = 0.50
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_a50_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x$coefs.est["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x$coefs.est["Tx_2", 1] - beta_tx[2]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(2, length(pred)), TERR1, TERR2)
  terr <- rbind(terr, terr2)
  # EN: alpha = 0.75
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_a75_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x$coefs.est["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x$coefs.est["Tx_2", 1] - beta_tx[2]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(3, length(pred)), TERR1, TERR2)
  terr <- rbind(terr, terr2)
    # Bayesian Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("Sum_Blasso1_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x["Tx_2", 1] - beta_tx[2]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(4, length(pred)), TERR1, TERR2)
  terr <- rbind(terr, terr2)
  # MARS
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("Earth1_cv_", tx[i], "_", n_tx[i], "_", 
      rep[r], ".rds", sep = "")))
  }
  TERR <- sapply(pred, function(x) 
    if (length(x) == 0) 
      NA 
    else if ("Tx_1" %in% rownames(x$coefficients)) 
      (x$coefficients["Tx_1", 1] - beta_tx[1]) else (0 - beta_tx[1]))
  TERR1 <- sapply(TERR, function(x) if (is.na(x)) NA else x)
  TERR <- sapply(pred, function(x) 
    if (length(x) == 0) 
      NA 
    else if ("Tx_2" %in% rownames(x$coefficients)) 
      (x$coefficients["Tx_2", 1] - beta_tx[2]) else (0 - beta_tx[2]))
  TERR2 <- sapply(TERR, function(x) if (is.na(x)) NA else x)
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(5, length(pred)), TERR1, TERR2)
  terr <- rbind(terr, terr2)
  print(i)
}

colnames(terr) <- c("Tx", "N_Tx", "Model", "TERR1", "TERR2")
terr <- as.data.frame(terr)
write.table(terr, paste(f2, "tx3_terr_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)


terr <- NULL
rep <- c(1, 2501, 5001, 7501)
tx <- rep(5, 3)
beta_tx <- c(0.75, -0.75, 1.25, -1.25)
n_tx <- c(100, 250, 500)
for (i in 1:length(n_tx)) {
  # Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x$coefs.est["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x$coefs.est["Tx_2", 1] - beta_tx[2]))
  TERR3 <- sapply(pred, function(x) (x$coefs.est["Tx_3", 1] - beta_tx[3]))
  TERR4 <- sapply(pred, function(x) (x$coefs.est["Tx_4", 1] - beta_tx[4]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(0, length(pred)), TERR1, TERR2, TERR3, TERR4)
  terr <- rbind(terr, terr2)
  # EN: alpha = 0.25
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_a25_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x$coefs.est["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x$coefs.est["Tx_2", 1] - beta_tx[2]))
  TERR3 <- sapply(pred, function(x) (x$coefs.est["Tx_3", 1] - beta_tx[3]))
  TERR4 <- sapply(pred, function(x) (x$coefs.est["Tx_4", 1] - beta_tx[4]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(1, length(pred)), TERR1, TERR2, TERR3, TERR4)
  terr <- rbind(terr, terr2)
  # EN: alpha = 0.50
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_a50_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x$coefs.est["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x$coefs.est["Tx_2", 1] - beta_tx[2]))
  TERR3 <- sapply(pred, function(x) (x$coefs.est["Tx_3", 1] - beta_tx[3]))
  TERR4 <- sapply(pred, function(x) (x$coefs.est["Tx_4", 1] - beta_tx[4]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(2, length(pred)), TERR1, TERR2, TERR3, TERR4)
  terr <- rbind(terr, terr2)
  # EN: alpha = 0.75
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN1_boot_", tx[i], "_", n_tx[i], "_a75_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x$coefs.est["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x$coefs.est["Tx_2", 1] - beta_tx[2]))
  TERR3 <- sapply(pred, function(x) (x$coefs.est["Tx_3", 1] - beta_tx[3]))
  TERR4 <- sapply(pred, function(x) (x$coefs.est["Tx_4", 1] - beta_tx[4]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(3, length(pred)), TERR1, TERR2, TERR3, TERR4)
  terr <- rbind(terr, terr2)
    # Bayesian Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("Sum_Blasso1_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  TERR1 <- sapply(pred, function(x) (x["Tx_1", 1] - beta_tx[1]))
  TERR2 <- sapply(pred, function(x) (x["Tx_2", 1] - beta_tx[2]))
  TERR3 <- sapply(pred, function(x) (x["Tx_3", 1] - beta_tx[3]))
  TERR4 <- sapply(pred, function(x) (x["Tx_4", 1] - beta_tx[4]))
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(4, length(pred)), TERR1, TERR2, TERR3, TERR4)
  terr <- rbind(terr, terr2)
  # MARS
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("Earth1_cv_", tx[i], "_", n_tx[i], "_", 
      rep[r], ".rds", sep = "")))
  }
  TERR <- sapply(pred, function(x) {
    if ("Tx_1" %in% rownames(x$coefficients)) 
      (x$coefficients["Tx_1", 1] - beta_tx[1]) else 
        (0 - beta_tx[1])
  })
  TERR1 <- sapply(TERR, function(x) if (is.na(x)) NA else x)
  TERR <- sapply(pred, function(x) {
    if ("Tx_2" %in% rownames(x$coefficients)) 
      (x$coefficients["Tx_2", 1] - beta_tx[2]) else 
        (0 - beta_tx[2])
  })
  TERR2 <- sapply(TERR, function(x) x)
  TERR <- sapply(pred, function(x) {
    if ("Tx_3" %in% rownames(x$coefficients)) 
      (x$coefficients["Tx_3", 1] - beta_tx[3]) else 
        (0 - beta_tx[3])
  })
  TERR3 <- sapply(TERR, function(x) x)
  TERR <- sapply(pred, function(x) {
    if ("Tx_4" %in% rownames(x$coefficients)) 
      (x$coefficients["Tx_4", 1] - beta_tx[4]) else
        (0 - beta_tx[4])
  })
  TERR4 <- sapply(TERR, function(x) x)
  terr2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(5, length(pred)), TERR1, TERR2, TERR3, TERR4)
  terr <- rbind(terr, terr2)
  print(i)
}

colnames(terr) <- c("Tx", "N_Tx", "Model", "TERR1", "TERR2", "TERR3", "TERR4")
terr <- as.data.frame(terr)
write.table(terr, paste(f2, "tx5_terr_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

# Find prediction errors and R^2 values

# 1st element of file indices
rep <- c(1, 2501, 5001, 7501)

# treatment group:
# number of treatment groups
tx <- c(rep(2, 3), rep(3, 3), rep(5, 3))
# sample size per tx group
n_tx <- rep(c(100, 250, 500), 3)


# RMSE yhat and R^2: testing sets
rmse <- NULL

for (i in 1:length(tx)) {
  Vars2 <- list()
  for (r in 1:length(rep)) {
    Vars2 <- append(Vars2, readRDS(paste("Vars_", tx[i], "_", n_tx[i], "_", 
                                       rep[r] + 10000, ".rds", sep = "")))
  }
  # Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN2_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(0, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
  # EN: alpha = 0.25
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN2_", tx[i], "_", n_tx[i], "_a25_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(1, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
  # EN: alpha = 0.50
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN2_", tx[i], "_", n_tx[i], "_a50_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(2, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
  # EN: alpha = 0.75
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("EN2_", tx[i], "_", n_tx[i], "_a75_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(3, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
    # Bayesian Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("Blasso2_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(4, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
  # MARS
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("Earth2_", tx[i], "_", n_tx[i], "_", 
      rep[r], ".rds", sep = "")))
  }
  pred <- lapply(pred, function(x) if (length(x) == 0) NA else x)
  mse <- mapply(function(x, y) 
    if (is.null(dim(y))) NA else mean((x[, 1] - y[, 1])^2, na.rm = TRUE)/var(x[, 1]), 
    Vars2, pred)
  mse2 <- sapply(mse, function(x) if (is.na(x)) NA else sqrt(x))
  r2 <- mapply(function(x, y) 
    if (is.null(dim(y))) NA else cor(x[, 1], y[, 1])^2, 
    Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(5, length(pred)), mse2, r2)
  rmse <- rbind(rmse, rmse2)
  print(i)
}

colnames(rmse) <- c("Tx", "N_Tx", "Model", "RMSE", "R2")
rmse <- as.data.frame(rmse)
write.table(rmse, paste(f2, "rmse_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

# RMSE yhat and R^2: training sets
rmse <- NULL

for (i in 1:length(tx)) {
  Vars2 <- list()
  for (r in 1:length(rep)) {
    Vars2 <- append(Vars2, readRDS(paste("Vars_", tx[i], "_", n_tx[i], "_", 
                                         rep[r], ".rds", sep = "")))
  }
  # Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("train_EN2_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(0, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
  # EN: alpha = 0.25
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("train_EN2_", tx[i], "_", n_tx[i], "_a25_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(1, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
  # EN: alpha = 0.50
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("train_EN2_", tx[i], "_", n_tx[i], "_a50_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(2, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
  # EN: alpha = 0.75
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("train_EN2_", tx[i], "_", n_tx[i], "_a75_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(3, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
    # Bayesian Lasso
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("train_Blasso2_", tx[i], "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  mse <- mapply(function(x, y) mean((x[, 1] - y$y.fitted)^2, na.rm = TRUE)/var(x[, 1]), Vars2, pred)
  r2 <- mapply(function(x, y) cor(x[, 1], y$y.fitted)^2, Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(4, length(pred)), sqrt(mse), r2)
  rmse <- rbind(rmse, rmse2)
  # MARS
  pred <- list()
  for (r in 1:length(rep)) {
    pred <- append(pred, readRDS(paste("train_Earth2_", tx[i], "_", n_tx[i], "_", 
      rep[r], ".rds", sep = "")))
  }
  pred <- lapply(pred, function(x) if (length(x) == 0) NA else x)
  mse <- mapply(function(x, y) 
    if (is.null(dim(y))) NA else mean((x[, 1] - y[, 1])^2, na.rm = TRUE)/var(x[, 1]), 
    Vars2, pred)
  mse2 <- sapply(mse, function(x) if (is.na(x)) NA else sqrt(x))
  r2 <- mapply(function(x, y) 
    if (is.null(dim(y))) NA else cor(x[, 1], y[, 1])^2, 
    Vars2, pred)
  rmse2 <- cbind(rep(tx[i], length(pred)), rep(n_tx[i], length(pred)), 
    rep(5, length(pred)), mse2, r2)
  rmse <- rbind(rmse, rmse2)
  print(i)
}

colnames(rmse) <- c("Tx", "N_Tx", "Model", "RMSE", "R2")
rmse <- as.data.frame(rmse)
write.table(rmse, paste(f2, "train_rmse_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

# Find false and true positive rates

# 1st element of file indices
rep <- c(1, 2501, 5001, 7501)

alpha1 <- 0.001
alpha2 <- 0.05
alpha3 <- 0.01

## False positive rate = # false H_0 rejects/# zero effects = # false H_0 rejects/40
# fpr1 pvalue < alpha1
# fpr2 pvalue < alpha2
# fpr3 pvalue < alpha3

## True positive rate = # true H_0 rejects/# non-zero effects = # true H_0 rejects/14
# tpr1 pvalue < alpha1
# tpr2 pvalue < alpha2
# tpr3 pvalue < alpha3

# indices of bootstrap CI in coes.est for EN models: 5, 6 = percentile; 8, 9 = normal

# Percentile CI:

ind1 <- 5; ind2 <- 6

# 2 tx groups:

# number of treatment groups
tx <- rep(2, 3)
# sample size per tx group
n_tx <- c(100, 250, 500)

False <- c(paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""),
  paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""))
True <- c(paste("C1M", 1:10, sep = ""), "Sex", "Age", "Tx")

pr_2 <- NULL

# for tx = 2
for (i in 1:length(n_tx)) {
  # Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 2, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(0, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  # EN: alpha = 0.25
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 2, "_", n_tx[i], "_a25_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(1, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  # EN: alpha = 0.5
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 2, "_", n_tx[i], "_a50_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(2, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  # EN: alpha = 0.75
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 2, "_", n_tx[i], "_a75_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(3, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  # Bayesian Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Sum_Blasso1_", 2, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr1 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha1)/40)
  fpr2 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha2)/40)
  fpr3 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha3)/40)
  tpr1 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha1)
      if ((x["Race_1", "pvalue"] < alpha1) | 
          (x["Race_2", "pvalue"] < alpha1)) m <- m + 1
      m/14
    })
  tpr2 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha2)
      if ((x["Race_1", "pvalue"] < alpha2) | 
          (x["Race_2", "pvalue"] < alpha2)) m <- m + 1
      m/14
    })
  tpr3 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha3)
      if ((x["Race_1", "pvalue"] < alpha3) | 
          (x["Race_2", "pvalue"] < alpha3)) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(4, length(Model)), fpr1, fpr2, fpr3, tpr1, tpr2, tpr3)
  pr_2 <- rbind(pr_2, pr)
  # MARS
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Earth1_cv_", 2, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA
    } else {
      m <- 0
      for (j in 2:5) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/40
    } 
  })
  tpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA 
    } else {
      m <- 0
      if ("Sex" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Sex", 1]), 3) > 0) m <- m + 1
      }
      if ("Age" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Age", 1]), 3) > 0) m <- m + 1
      }
      if ("Tx" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx", 1]), 3) > 0) m <- m + 1
      }
      p1 <- 1; p2 <- 1
      if ("Race_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Race_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_2", 1]), 3) > 0) p2 <- 0
      }
      if (p1 == 0 | p2 == 0) m <- m + 1
      for (j in 1) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/14
    }
  })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(5, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  print(i)
}

colnames(pr_2) <- c("Tx", "N_Tx", "Model", "fpr1", "fpr2", "fpr3", "tpr1", "tpr2", "tpr3")
pr_2 <- as.data.frame(pr_2)
write.table(pr_2, paste(f2, "perc_pr_2_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

# 3 tx groups:
# number of treatment groups
tx <- rep(3, 3)
# sample size per tx group
n_tx <- c(100, 250, 500)

False <- c(paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""),
  paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""))
True <- c(paste("C1M", 1:10, sep = ""), "Sex", "Age")

pr_3 <- NULL

for (i in 1:length(n_tx)) {
  # Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 3, "_", n_tx[i], "_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(0, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  # EN: alpha = 0.25
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 3, "_", n_tx[i], "_a25_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(1, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  # EN: alpha = 0.5
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 3, "_", n_tx[i], "_a50_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(2, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  # EN: alpha = 0.75
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 3, "_", n_tx[i], "_a75_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(3, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  # Bayesian Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Sum_Blasso1_", 3, "_", n_tx[i], "_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr1 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha1)/40)
  fpr2 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha2)/40)
  fpr3 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha3)/40)
  tpr1 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha1)
      if ((x["Race_1", "pvalue"] < alpha1) | 
          (x["Race_2", "pvalue"] < alpha1)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha1) | 
          (x["Tx_2", "pvalue"] < alpha1)) m <- m + 1
      m/14
    })
  tpr2 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha2)
      if ((x["Race_1", "pvalue"] < alpha2) | 
          (x["Race_2", "pvalue"] < alpha2)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha2) | 
          (x["Tx_2", "pvalue"] < alpha2)) m <- m + 1
      m/14
    })
  tpr3 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha3)
      if ((x["Race_1", "pvalue"] < alpha3) | 
          (x["Race_2", "pvalue"] < alpha3)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha3) | 
          (x["Tx_2", "pvalue"] < alpha3)) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(4, length(Model)), fpr1, fpr2, fpr3, tpr1, tpr2, tpr3)
  pr_3 <- rbind(pr_3, pr)
  # MARS
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Earth1_cv_", 3, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA
    } else {
      m <- 0
      for (j in 2:5) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/40
    } 
  })
  tpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA 
    } else {
      m <- 0
      if ("Sex" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Sex", 1]), 3) > 0) m <- m + 1
      }
      if ("Age" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Age", 1]), 3) > 0) m <- m + 1
      }
      p1 <- 1; p2 <- 1
      if ("Tx_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Tx_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_2", 1]), 3) > 0) p2 <- 0
      }
      if (p1 == 0 | p2 == 0) m <- m + 1
      p1 <- 1; p2 <- 1
      if ("Race_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Race_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_2", 1]), 3) > 0) p2 <- 0
      }
      if (p1 == 0 | p2 == 0) m <- m + 1
      for (j in 1) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/14
    }
  })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(5, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  print(i)
}

colnames(pr_3) <- c("Tx", "N_Tx", "Model", "fpr1", "fpr2", "fpr3", "tpr1", "tpr2", "tpr3")
pr_3 <- as.data.frame(pr_3)
write.table(pr_3, paste(f2, "perc_pr_3_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

# 5 tx groups:
# number of treatment groups
tx <- rep(5, 3)
# sample size per tx group
n_tx <- c(100, 250, 500)

False <- c(paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""),
           paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""))
True <- c(paste("C1M", 1:10, sep = ""), "Sex", "Age")

pr_5 <- NULL

for (i in 1:length(n_tx)) {
  # Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 5, "_", n_tx[i], "_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2])) | 
          !(between(0, x$coefs.est["Tx_3", ind1], x$coefs.est["Tx_3", ind2])) | 
          !(between(0, x$coefs.est["Tx_4", ind1], x$coefs.est["Tx_4", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(0, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  # EN: alpha = 0.25
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 5, "_", n_tx[i], "_a25_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2])) | 
          !(between(0, x$coefs.est["Tx_3", ind1], x$coefs.est["Tx_3", ind2])) | 
          !(between(0, x$coefs.est["Tx_4", ind1], x$coefs.est["Tx_4", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(1, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  # EN: alpha = 0.5
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 5, "_", n_tx[i], "_a50_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2])) | 
          !(between(0, x$coefs.est["Tx_3", ind1], x$coefs.est["Tx_3", ind2])) | 
          !(between(0, x$coefs.est["Tx_4", ind1], x$coefs.est["Tx_4", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(2, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  # EN: alpha = 0.75
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 5, "_", n_tx[i], "_a75_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2])) | 
          !(between(0, x$coefs.est["Tx_3", ind1], x$coefs.est["Tx_3", ind2])) | 
          !(between(0, x$coefs.est["Tx_4", ind1], x$coefs.est["Tx_4", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(3, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  # Bayesian Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Sum_Blasso1_", 5, "_", n_tx[i], "_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr1 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha1)/40)
  fpr2 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha2)/40)
  fpr3 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha3)/40)
  tpr1 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha1)
      if ((x["Race_1", "pvalue"] < alpha1) | 
          (x["Race_2", "pvalue"] < alpha1)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha1) | 
          (x["Tx_2", "pvalue"] < alpha1) | 
          (x["Tx_3", "pvalue"] < alpha1) | 
          (x["Tx_4", "pvalue"] < alpha1)) m <- m + 1
      m/14
    })
  tpr2 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha2)
      if ((x["Race_1", "pvalue"] < alpha2) | 
          (x["Race_2", "pvalue"] < alpha2)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha2) | 
          (x["Tx_2", "pvalue"] < alpha2) | 
          (x["Tx_3", "pvalue"] < alpha2) | 
          (x["Tx_4", "pvalue"] < alpha2)) m <- m + 1
      m/14
    })
  tpr3 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha3)
      if ((x["Race_1", "pvalue"] < alpha3) | 
          (x["Race_2", "pvalue"] < alpha3)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha3) | 
          (x["Tx_2", "pvalue"] < alpha3) | 
          (x["Tx_3", "pvalue"] < alpha3) | 
          (x["Tx_4", "pvalue"] < alpha3)) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(4, length(Model)), fpr1, fpr2, fpr3, tpr1, tpr2, tpr3)
  pr_5 <- rbind(pr_5, pr)
  # MARS
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Earth1_cv_", 5, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA
    } else {
      m <- 0
      for (j in 2:5) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/40
    } 
  })
  tpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA 
    } else {
      m <- 0
      if ("Sex" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Sex", 1]), 3) > 0) m <- m + 1
      }
      if ("Age" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Age", 1]), 3) > 0) m <- m + 1
      }
      p1 <- 1; p2 <- 1; p3 <- 1; p4 <- 1
      if ("Tx_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Tx_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_2", 1]), 3) > 0) p2 <- 0
      }
      if ("Tx_3" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_3", 1]), 3) > 0) p3 <- 0
      }
      if ("Tx_4" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_4", 1]), 3) > 0) p4 <- 0
      }
      if (p1 == 0 | p2 == 0 | p3 == 0 | p4 == 0) m <- m + 1
      p1 <- 1; p2 <- 1
      if ("Race_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Race_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_2", 1]), 3) > 0) p2 <- 0
      }
      if (p1 == 0 | p2 == 0) m <- m + 1
      for (j in 1) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/14
    }
  })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(5, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  print(i)
}

colnames(pr_5) <- c("Tx", "N_Tx", "Model", "fpr1", "fpr2", "fpr3", "tpr1", "tpr2", "tpr3")
pr_5 <- as.data.frame(pr_5)
write.table(pr_5, paste(f2, "perc_pr_5_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

# Normal-theory CI:

ind1 <- 8; ind2 <- 9

# 2 tx groups:

# number of treatment groups
tx <- rep(2, 3)
# sample size per tx group
n_tx <- c(100, 250, 500)

False <- c(paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""),
  paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""))
True <- c(paste("C1M", 1:10, sep = ""), "Sex", "Age", "Tx")

pr_2 <- NULL

# for tx = 2
for (i in 1:length(n_tx)) {
  # Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 2, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(0, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  # EN: alpha = 0.25
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 2, "_", n_tx[i], "_a25_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(1, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  # EN: alpha = 0.5
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 2, "_", n_tx[i], "_a50_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(2, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  # EN: alpha = 0.75
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 2, "_", n_tx[i], "_a75_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(3, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  # Bayesian Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Sum_Blasso1_", 2, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr1 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha1)/40)
  fpr2 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha2)/40)
  fpr3 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha3)/40)
  tpr1 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha1)
      if ((x["Race_1", "pvalue"] < alpha1) | 
          (x["Race_2", "pvalue"] < alpha1)) m <- m + 1
      m/14
    })
  tpr2 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha2)
      if ((x["Race_1", "pvalue"] < alpha2) | 
          (x["Race_2", "pvalue"] < alpha2)) m <- m + 1
      m/14
    })
  tpr3 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha3)
      if ((x["Race_1", "pvalue"] < alpha3) | 
          (x["Race_2", "pvalue"] < alpha3)) m <- m + 1
      m/14
    })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(4, length(Model)), fpr1, fpr2, fpr3, tpr1, tpr2, tpr3)
  pr_2 <- rbind(pr_2, pr)
  # MARS
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Earth1_cv_", 2, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA
    } else {
      m <- 0
      for (j in 2:5) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/40
    } 
  })
  tpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA 
    } else {
      m <- 0
      if ("Sex" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Sex", 1]), 3) > 0) m <- m + 1
      }
      if ("Age" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Age", 1]), 3) > 0) m <- m + 1
      }
      if ("Tx" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx", 1]), 3) > 0) m <- m + 1
      }
      p1 <- 1; p2 <- 1
      if ("Race_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Race_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_2", 1]), 3) > 0) p2 <- 0
      }
      if (p1 == 0 | p2 == 0) m <- m + 1
      for (j in 1) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/14
    }
  })
  pr <- cbind(rep(2, length(Model)), rep(n_tx[i], length(Model)), 
              rep(5, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_2 <- rbind(pr_2, pr)
  print(i)
}

colnames(pr_2) <- c("Tx", "N_Tx", "Model", "fpr1", "fpr2", "fpr3", "tpr1", "tpr2", "tpr3")
pr_2 <- as.data.frame(pr_2)
write.table(pr_2, paste(f2, "pr_2_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

# 3 tx groups:
# number of treatment groups
tx <- rep(3, 3)
# sample size per tx group
n_tx <- c(100, 250, 500)

False <- c(paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""),
  paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""))
True <- c(paste("C1M", 1:10, sep = ""), "Sex", "Age")

pr_3 <- NULL

for (i in 1:length(n_tx)) {
  # Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 3, "_", n_tx[i], "_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(0, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  # EN: alpha = 0.25
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 3, "_", n_tx[i], "_a25_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(1, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  # EN: alpha = 0.5
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 3, "_", n_tx[i], "_a50_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(2, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  # EN: alpha = 0.75
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 3, "_", n_tx[i], "_a75_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(3, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  # Bayesian Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Sum_Blasso1_", 3, "_", n_tx[i], "_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr1 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha1)/40)
  fpr2 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha2)/40)
  fpr3 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha3)/40)
  tpr1 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha1)
      if ((x["Race_1", "pvalue"] < alpha1) | 
          (x["Race_2", "pvalue"] < alpha1)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha1) | 
          (x["Tx_2", "pvalue"] < alpha1)) m <- m + 1
      m/14
    })
  tpr2 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha2)
      if ((x["Race_1", "pvalue"] < alpha2) | 
          (x["Race_2", "pvalue"] < alpha2)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha2) | 
          (x["Tx_2", "pvalue"] < alpha2)) m <- m + 1
      m/14
    })
  tpr3 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha3)
      if ((x["Race_1", "pvalue"] < alpha3) | 
          (x["Race_2", "pvalue"] < alpha3)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha3) | 
          (x["Tx_2", "pvalue"] < alpha3)) m <- m + 1
      m/14
    })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(4, length(Model)), fpr1, fpr2, fpr3, tpr1, tpr2, tpr3)
  pr_3 <- rbind(pr_3, pr)
  # MARS
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Earth1_cv_", 3, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA
    } else {
      m <- 0
      for (j in 2:5) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/40
    } 
  })
  tpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA 
    } else {
      m <- 0
      if ("Sex" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Sex", 1]), 3) > 0) m <- m + 1
      }
      if ("Age" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Age", 1]), 3) > 0) m <- m + 1
      }
      p1 <- 1; p2 <- 1
      if ("Tx_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Tx_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_2", 1]), 3) > 0) p2 <- 0
      }
      if (p1 == 0 | p2 == 0) m <- m + 1
      p1 <- 1; p2 <- 1
      if ("Race_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Race_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_2", 1]), 3) > 0) p2 <- 0
      }
      if (p1 == 0 | p2 == 0) m <- m + 1
      for (j in 1) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/14
    }
  })
  pr <- cbind(rep(3, length(Model)), rep(n_tx[i], length(Model)), 
              rep(5, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_3 <- rbind(pr_3, pr)
  print(i)
}

colnames(pr_3) <- c("Tx", "N_Tx", "Model", "fpr1", "fpr2", "fpr3", "tpr1", "tpr2", "tpr3")
pr_3 <- as.data.frame(pr_3)
write.table(pr_3, paste(f2, "pr_3_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

# 5 tx groups:
# number of treatment groups
tx <- rep(5, 3)
# sample size per tx group
n_tx <- c(100, 250, 500)

False <- c(paste("C2M", 1:10, sep = ""), paste("C3M", 1:10, sep = ""),
           paste("C4M", 1:10, sep = ""), paste("C5M", 1:10, sep = ""))
True <- c(paste("C1M", 1:10, sep = ""), "Sex", "Age")

pr_5 <- NULL

for (i in 1:length(n_tx)) {
  # Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 5, "_", n_tx[i], "_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2])) | 
          !(between(0, x$coefs.est["Tx_3", ind1], x$coefs.est["Tx_3", ind2])) | 
          !(between(0, x$coefs.est["Tx_4", ind1], x$coefs.est["Tx_4", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(0, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  # EN: alpha = 0.25
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 5, "_", n_tx[i], "_a25_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2])) | 
          !(between(0, x$coefs.est["Tx_3", ind1], x$coefs.est["Tx_3", ind2])) | 
          !(between(0, x$coefs.est["Tx_4", ind1], x$coefs.est["Tx_4", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(1, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  # EN: alpha = 0.5
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 5, "_", n_tx[i], "_a50_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2])) | 
          !(between(0, x$coefs.est["Tx_3", ind1], x$coefs.est["Tx_3", ind2])) | 
          !(between(0, x$coefs.est["Tx_4", ind1], x$coefs.est["Tx_4", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(2, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  # EN: alpha = 0.75
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("EN1_boot_", 5, "_", n_tx[i], "_a75_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% False, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      m/40
    })
  tpr <- sapply(Model, 
    function(x) {
      x2 <- x$coefs.est[rownames(x$coefs.est) %in% True, ]
      m <- 0
      for (j in 1:nrow(x2)) {
        if (!(between(0, x2[j, ind1], x2[j, ind2]))) m <- m + 1
      }
      if (!(between(0, x$coefs.est["Race_1", ind1], x$coefs.est["Race_1", ind2])) | 
          !(between(0, x$coefs.est["Race_2", ind1], x$coefs.est["Race_2", ind2]))) m <- m + 1
      if (!(between(0, x$coefs.est["Tx_1", ind1], x$coefs.est["Tx_1", ind2])) | 
          !(between(0, x$coefs.est["Tx_2", ind1], x$coefs.est["Tx_2", ind2])) | 
          !(between(0, x$coefs.est["Tx_3", ind1], x$coefs.est["Tx_3", ind2])) | 
          !(between(0, x$coefs.est["Tx_4", ind1], x$coefs.est["Tx_4", ind2]))) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(3, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  # Bayesian Lasso
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Sum_Blasso1_", 5, "_", n_tx[i], "_", 
                                         rep[r], ".rds", sep = "")))
  }
  fpr1 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha1)/40)
  fpr2 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha2)/40)
  fpr3 <- sapply(Model, 
    function(x) sum(x[rownames(x) %in% False, "pvalue"] < alpha3)/40)
  tpr1 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha1)
      if ((x["Race_1", "pvalue"] < alpha1) | 
          (x["Race_2", "pvalue"] < alpha1)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha1) | 
          (x["Tx_2", "pvalue"] < alpha1) | 
          (x["Tx_3", "pvalue"] < alpha1) | 
          (x["Tx_4", "pvalue"] < alpha1)) m <- m + 1
      m/14
    })
  tpr2 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha2)
      if ((x["Race_1", "pvalue"] < alpha2) | 
          (x["Race_2", "pvalue"] < alpha2)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha2) | 
          (x["Tx_2", "pvalue"] < alpha2) | 
          (x["Tx_3", "pvalue"] < alpha2) | 
          (x["Tx_4", "pvalue"] < alpha2)) m <- m + 1
      m/14
    })
  tpr3 <- sapply(Model, 
    function(x) {
      m <- sum(x[rownames(x) %in% True, "pvalue"] < alpha3)
      if ((x["Race_1", "pvalue"] < alpha3) | 
          (x["Race_2", "pvalue"] < alpha3)) m <- m + 1
      if ((x["Tx_1", "pvalue"] < alpha3) | 
          (x["Tx_2", "pvalue"] < alpha3) | 
          (x["Tx_3", "pvalue"] < alpha3) | 
          (x["Tx_4", "pvalue"] < alpha3)) m <- m + 1
      m/14
    })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(4, length(Model)), fpr1, fpr2, fpr3, tpr1, tpr2, tpr3)
  pr_5 <- rbind(pr_5, pr)
  # MARS
  Model <- list()
  for (r in 1:length(rep)) {
    Model <- append(Model, readRDS(paste("Earth1_cv_", 5, "_", n_tx[i], "_", 
                                       rep[r], ".rds", sep = "")))
  }
  fpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA
    } else {
      m <- 0
      for (j in 2:5) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/40
    } 
  })
  tpr <- sapply(Model, function(x) {
    if (length(x) == 0) {
      NA 
    } else {
      m <- 0
      if ("Sex" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Sex", 1]), 3) > 0) m <- m + 1
      }
      if ("Age" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Age", 1]), 3) > 0) m <- m + 1
      }
      p1 <- 1; p2 <- 1; p3 <- 1; p4 <- 1
      if ("Tx_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Tx_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_2", 1]), 3) > 0) p2 <- 0
      }
      if ("Tx_3" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_3", 1]), 3) > 0) p3 <- 0
      }
      if ("Tx_4" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Tx_4", 1]), 3) > 0) p4 <- 0
      }
      if (p1 == 0 | p2 == 0 | p3 == 0 | p4 == 0) m <- m + 1
      p1 <- 1; p2 <- 1
      if ("Race_1" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_1", 1]), 3) > 0) p1 <- 0
      }
      if ("Race_2" %in% rownames(x$coefficients)) {
        if (round(abs(x$coefficients["Race_2", 1]), 3) > 0) p2 <- 0
      }
      if (p1 == 0 | p2 == 0) m <- m + 1
      for (j in 1) {
        for (k in 1:10) {
          p1 <- 1; p2 <- 1
          if (paste("C", j, "M", k, "1", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "1", sep = ""), 1]), 3) > 0) p1 <- 0
          }
          if (paste("C", j, "M", k, "2", sep = "") %in% rownames(x$coefficients)) {
            if (round(abs(x$coefficients[paste("C", j, "M", k, "2", sep = ""), 1]), 3) > 0) p2 <- 0
          }
          if (p1 == 0 | p2 == 0) m <- m + 1
        }
      }
      m/14
    }
  })
  pr <- cbind(rep(5, length(Model)), rep(n_tx[i], length(Model)), 
              rep(5, length(Model)), fpr, fpr, fpr, tpr, tpr, tpr)
  pr_5 <- rbind(pr_5, pr)
  print(i)
}

colnames(pr_5) <- c("Tx", "N_Tx", "Model", "fpr1", "fpr2", "fpr3", "tpr1", "tpr2", "tpr3")
pr_5 <- as.data.frame(pr_5)
write.table(pr_5, paste(f2, "pr_5_betas_2.txt", sep = ""), col.names = T, row.names = F, 
            sep = "\t", quote = F)

# Five simulated and predicted response values
ind <- c(1, 500, 1000, 1500, 2000)

# treatment group:
# number of treatment groups
tx <- c(rep(2, 3), rep(3, 3), rep(5, 3))
# sample size per tx group
n_tx <- rep(c(100, 250, 500), 3)

for (i in 1:length(tx)) {
  Vars2 <- list()
  for (r in 1:length(ind)) {
    Vars2[[r]] <- matrix(0, tx[i] * n_tx[i], 7)
    colnames(Vars2[[r]]) <- c("y", "EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS")
    Vars2[[r]][, 1] <- readRDS(paste("Vars_", tx[i], "_", n_tx[i], "_", 10001, ".rds", sep = ""))[[ind[r]]][, 1]
    Vars2[[r]][, 2] <- readRDS(paste("EN2_", tx[i], "_", n_tx[i], "_a25_", 1, ".rds", sep = ""))[[ind[r]]]$y.fitted
    Vars2[[r]][, 3] <- readRDS(paste("EN2_", tx[i], "_", n_tx[i], "_a50_", 1, ".rds", sep = ""))[[ind[r]]]$y.fitted
    Vars2[[r]][, 4] <- readRDS(paste("EN2_", tx[i], "_", n_tx[i], "_a75_", 1, ".rds", sep = ""))[[ind[r]]]$y.fitted
    Vars2[[r]][, 5] <- readRDS(paste("EN2_", tx[i], "_", n_tx[i], "_", 1, ".rds", sep = ""))[[ind[r]]]$y.fitted
    Vars2[[r]][, 6] <- readRDS(paste("Blasso2_", tx[i], "_", n_tx[i], "_", 1, ".rds", sep = ""))[[ind[r]]]$y.fitted
    Vars2[[r]][, 7] <- readRDS(paste("Earth2_", tx[i], "_", n_tx[i], "_", 1, ".rds", sep = ""))[[ind[r]]][, 1]
  }
  saveRDS(Vars2, paste(f2, "Vars2_", tx[i], "_", n_tx[i], "_betas_2.rds", sep = ""))
}
