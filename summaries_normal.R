###############################################################################
## Code to summarize model results with Normal error terms
###############################################################################

library("ggplot2")
library("plyr")
library("xtable")
library("Rmisc")
library("gridExtra")
library("grid")


f1 <- "Normal_Errors/"

# function to extract legend from multiple ggplot2 graphs and place at bottom
# using grid.arrange;
# code adapted from https://github.com/tidyverse/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs

grid_arrange_shared_legend <- function(...) {
  plots <- list(...)
  g <- ggplotGrob(plots[[1]] + 
                    theme(legend.position="bottom", legend.direction = "horizontal", 
                          legend.key.size = unit(1, "cm"), 
                          legend.margin=margin(t=-0.15, r=0, b=0, l=0, unit="cm")) + 
                    guides(colour = guide_legend(nrow = 1)))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  plots2 <- arrangeGrob(..., ncol = 1)
  grid.arrange(plots2, legend, ncol = 1,
               heights = unit.c(unit(1, "npc") - lheight, lheight))
}

# summary of outcomes (Suppl. Info. Table 1)
y1 <- read.table(paste(f1, "Sum_y_betas1.txt", sep = ""), header = T, sep = "\t", quote = "")
y2 <- read.table(paste(f1, "Sum_y_betas2.txt", sep = ""), header = T, sep = "\t", quote = "")
y3 <- read.table(paste(f1, "Sum_y_betas3.txt", sep = ""), header = T, sep = "\t", quote = "")

y <- rbind(y1, y2, y3)
print(xtable(y, align = c("l", rep("c", 9)), digits = c(rep(0, 4), rep(3, 6))), 
      include.rownames = F, include.colnames = T)


# computation times: minutes per 1000 models (Table 2)
# treatment group:
# number of treatment groups
tx <- c(rep(2, 3), rep(3, 3), rep(5, 3))
# sample size per tx group
n_tx <- rep(c(100, 250, 500), 3)

tEN25 <- NULL
tEN50 <- NULL
tEN75 <- NULL
tL <- NULL
tBL <- NULL
tM <- NULL

for (i in 1:length(tx)) {
  for (k in 1:3) {
    tEN25 <- c(tEN25, readRDS(paste(f1, "Time_EN1_cv_", tx[i], "_a25_", n_tx[i], "_betas", k, ".rds", sep = "")))
    tEN50 <- c(tEN50, readRDS(paste(f1, "Time_EN1_cv_", tx[i], "_a50_", n_tx[i], "_betas", k, ".rds", sep = "")))
    tEN75 <- c(tEN75, readRDS(paste(f1, "Time_EN1_cv_", tx[i], "_a75_", n_tx[i], "_betas", k, ".rds", sep = "")))
    tL <- c(tL, readRDS(paste(f1, "Time_EN1_cv_", tx[i], "_", n_tx[i], "_betas", k, ".rds", sep = "")))
    tBL <- c(tBL, readRDS(paste(f1, "Time_Blasso1_", tx[i], "_", n_tx[i], "_betas", k, ".rds", sep = "")))
    tM <- c(tM, readRDS(paste(f1, "Time_Earth1_cv_", tx[i], "_", n_tx[i], "_betas", k, ".rds", sep = "")))
    names(tEN25)[length(tEN25)] <- names(tEN50)[length(tEN25)] <- names(tEN75)[length(tEN25)] <- 
      names(tL)[length(tEN25)] <- names(tBL)[length(tEN25)] <- names(tM)[length(tEN25)] <- 
      paste(tx[i], "_", n_tx[i], "_", k, sep = "")
  }
}
times <- rbind(rbind(tEN25[1:9], tEN50[1:9], tEN75[1:9], tL[1:9], tBL[1:9], tM[1:9]), 
               rbind(tEN25[10:18], tEN50[10:18], tEN75[10:18], tL[10:18], tBL[10:18], tM[10:18]), 
               rbind(tEN25[19:27], tEN50[19:27], tEN75[19:27], tL[19:27], tBL[19:27], tM[19:27]))
times <- as.data.frame(apply(times, 2, function(x) x/(10*60)))
times$Model <- rep(c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"), 3)
print(xtable(times[, c(10, 1:9)], align = c("l", "l", rep("c", 9)), digits = c(0, 0, rep(2, 9))), 
             include.rownames = F, include.colnames = F)

# Mean absolute error beta_tx
# 2 tx groups
terr1 <- read.table(paste(f1, "tx2_terr_betas_1.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")
terr2 <- read.table(paste(f1, "tx2_terr_betas_2.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")
terr3 <- read.table(paste(f1, "tx2_terr_betas_3.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")

terr <- rbind(terr1, terr2, terr3)
terr$beta <- c(rep(1, nrow(terr)/3), rep(2, nrow(terr)/3), rep(3, nrow(terr)/3))

terr$beta <- factor(terr$beta)
terr$n <- factor(terr$Tx * terr$N_Tx)
terr$Model <- ifelse(terr$Model == 0, 4, ifelse(terr$Model == 4, 5, 
  ifelse(terr$Model == 5, 6, terr$Model)))

terr$Model <- factor(terr$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
terr$Tx <- factor(terr$Tx)
terr$N_Tx <- factor(terr$N_Tx)
terr$TERR <- ifelse(terr$beta == 1, terr$TERR/0.5, ifelse(terr$beta == 2, terr$TERR/1, terr$TERR/1.5))

# graph by beta and N_Tx (Figure 1)
sum_terr <- summarySE(data=terr, measurevar = "TERR", groupvars = c("beta", "N_Tx", "Model"), na.rm = TRUE, 
                      conf.interval = 0.95, .drop = TRUE)

pd <- position_dodge(0.1)
color <- c("#009900", "#CC0000", "#FF9933", "#990099", "#0000FF", "#666666")

beta_names <- list(
  '1'=expression(paste("Small Tx Effect: ", beta[Tx[1]], " = 0.5")),
  '2'=expression(paste("Medium Tx Effect: ", beta[Tx[1]], " = 1")),
  '3'=expression(paste("Large Tx Effect: ", beta[Tx[1]], " = 1.5"))
)
beta_labeller <- function(variable, value){
  return(beta_names[value])
}

terr_plot1 <- ggplot(sum_terr, 
                     aes(x = N_Tx, y = TERR, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = TERR-ci, ymax = TERR+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~beta, labeller = beta_labeller) + theme_bw() + 
  xlab("Sample Size per Tx Group") + 
  scale_y_continuous(expression(paste("MAE Proportion for ", beta[Tx[1]]))) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + 
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.position = "bottom", legend.direction = "horizontal", 
        legend.key.size = unit(1, "cm"), 
        legend.margin=margin(t=-0.4, r=0, b=0, l=0, unit="cm")) + guides(colour = guide_legend(nrow = 1))
terr_plot1
ggsave(paste(f4, "Fialkowski_Figure1.tiff", sep = ""), dpi = 800, terr_plot1, 
       width = 6.7, height = 3, units = "in")

# table by total sample size n (labeled by sample size per group in tables)
sum_terr <- summarySE(data=terr, measurevar = "TERR", groupvars = c("beta", "n", "Model"), na.rm = TRUE, 
                      conf.interval = 0.95, .drop = TRUE)
sum_terr <- sum_terr[order(sum_terr$beta, sum_terr$n, sum_terr$Model), ]

sum_terr$CI <- paste(formatC(round(sum_terr$TERR, 3), 3, format = "f"), " (", 
  formatC(round(sum_terr$TERR, 3) - round(sum_terr$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_terr$TERR, 3) + round(sum_terr$ci, 3), 3, format = "f"), ")", sep = "")
sum_terr$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)
sum_terr$n <- rep(c("200", rep(NA, 5), "500", rep(NA, 5), "1000", rep(NA, 5)), 3)

sum_terr2 <- cbind(sum_terr[sum_terr$beta == 1, ], sum_terr[sum_terr$beta == 2, ], 
  sum_terr[sum_terr$beta == 3, ])
print(xtable(sum_terr2[, c(3, 9, 19, 29)], align = c("l", "l", rep("c", 3))), 
      caption.placement = "top", include.rownames = FALSE, sanitize.text.function = function(x){x})

# 3 tx groups
terr1 <- read.table(paste(f1, "tx3_terr_betas_1.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")
terr2 <- read.table(paste(f1, "tx3_terr_betas_2.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")
terr3 <- read.table(paste(f1, "tx3_terr_betas_3.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")

terr <- rbind(terr1, terr2, terr3)
terr$beta <- c(rep(1, nrow(terr)/3), rep(2, nrow(terr)/3), rep(3, nrow(terr)/3))

terr$beta <- factor(terr$beta)
terr$n <- factor(terr$Tx * terr$N_Tx)
terr$Model <- ifelse(terr$Model == 0, 4, ifelse(terr$Model == 4, 5, 
  ifelse(terr$Model == 5, 6, terr$Model)))

terr$Model <- factor(terr$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
terr$Tx <- factor(terr$Tx)
terr$N_Tx <- factor(terr$N_Tx)
terr$TERR1 <- ifelse(terr$beta == 1, terr$TERR1/0.5, ifelse(terr$beta == 2, terr$TERR1/1, terr$TERR1/1.5))
terr$TERR2 <- ifelse(terr$beta == 1, terr$TERR2/0.25, ifelse(terr$beta == 2, terr$TERR2/0.5, terr$TERR2/0.75))

# graph by beta and N_Tx (Figure 2)
sum_terr1 <- summarySE(terr, measurevar = "TERR1", groupvars = c("beta", "N_Tx", "Model"))
sum_terr2 <- summarySE(terr, measurevar = "TERR2", groupvars = c("beta", "N_Tx", "Model"))

pd <- position_dodge(0.1)
color <- c("#009900", "#CC0000", "#FF9933", "#990099", "#0000FF", "#666666")

beta_names <- list(
  '1'=expression(paste("Small Tx Effect: ", beta[Tx[1]], " = 0.5")),
  '2'=expression(paste("Medium Tx Effect: ", beta[Tx[1]], " = 1")),
  '3'=expression(paste("Large Tx Effect: ", beta[Tx[1]], " = 1.5"))
)
beta_labeller <- function(variable, value){
  return(beta_names[value])
}

terr_plot1 <- ggplot(sum_terr1, 
                     aes(x = N_Tx, y = TERR1, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = TERR1-ci, ymax = TERR1+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~beta, labeller = beta_labeller) + theme_bw() + 
  xlab("") + 
  scale_y_continuous(expression(paste("MAE Proportion for ", beta[Tx[1]]))) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
terr_plot1

beta_names2 <- list(
  '1'=expression(paste("Small Tx Effect: ", beta[Tx[2]], " = -0.25")),
  '2'=expression(paste("Medium Tx Effect: ", beta[Tx[2]], " = -0.5")),
  '3'=expression(paste("Large Tx Effect: ", beta[Tx[2]], " = -0.75"))
)
beta_labeller2 <- function(variable, value){
  return(beta_names2[value])
}

terr_plot2 <- ggplot(sum_terr2, 
                     aes(x = N_Tx, y = TERR2, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = TERR2-ci, ymax = TERR2+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~beta, labeller = beta_labeller2) + theme_bw() + 
  xlab("Sample Size per Tx Group") + 
  scale_y_continuous(expression(paste("MAE Proportion for ", beta[Tx[2]]))) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
terr_plot2

Fig2 <- grid_arrange_shared_legend(terr_plot1, terr_plot2)
ggsave(paste(f4, "Fialkowski_Figure2.png", sep = ""), dpi = 800, Fig2, 
       width = 6.7, height = 5.1, units = "in")

# table by total sample size n (labeled by sample size per group in tables)
sum_terr <- summarySE(terr, measurevar = "TERR1", groupvars = c("beta", "n", "Model"), na.rm = TRUE)
sum_terr <- sum_terr[order(sum_terr$beta, sum_terr$n, sum_terr$Model), ]

sum_terr$CI_1 <- paste(formatC(round(sum_terr$TERR1, 3), 3, format = "f"), " (", 
  formatC(round(sum_terr$TERR1, 3) - round(sum_terr$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_terr$TERR1, 3) + round(sum_terr$ci, 3), 3, format = "f"), ")", sep = "")
sum_terr$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)
sum_terr$n <- rep(c("300", rep(NA, 5), "750", rep(NA, 5), "1500", rep(NA, 5)), 3)

sum_terr2 <- summarySE(terr, measurevar = "TERR2", groupvars = c("beta", "n", "Model"), na.rm = TRUE)
sum_terr2 <- sum_terr2[order(sum_terr2$beta, sum_terr2$n, sum_terr2$Model), ]

sum_terr2$CI_2 <- paste(formatC(round(sum_terr2$TERR2, 3), 3, format = "f"), " (", 
  formatC(round(sum_terr2$TERR2, 3) - round(sum_terr2$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_terr2$TERR2, 3) + round(sum_terr2$ci, 3), 3, format = "f"), ")", sep = "")

sum_terr3 <- cbind(sum_terr[sum_terr$beta == 1, c("Model", "n", "CI_1")], 
  sum_terr2[sum_terr2$beta == 1, "CI_2"], sum_terr[sum_terr$beta == 2, c("CI_1")], 
  sum_terr2[sum_terr2$beta == 2, "CI_2"], sum_terr[sum_terr$beta == 3, "CI_1"], 
  sum_terr2[sum_terr2$beta == 3, "CI_2"])

print(xtable(sum_terr3[, -2], align = c("l", "l", rep("c", 6))), 
      caption.placement = "top", include.rownames = FALSE, 
      include.colnames = FALSE, sanitize.text.function = function(x){x})

# 5 tx groups
terr1 <- read.table(paste(f1, "tx5_terr_betas_1.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")
terr2 <- read.table(paste(f1, "tx5_terr_betas_2.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")
terr3 <- read.table(paste(f1, "tx5_terr_betas_3.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")

terr <- rbind(terr1, terr2, terr3)
terr$beta <- c(rep(1, nrow(terr)/3), rep(2, nrow(terr)/3), rep(3, nrow(terr)/3))

terr$beta <- factor(terr$beta)
terr$n <- factor(terr$Tx * terr$N_Tx)
terr$Model <- ifelse(terr$Model == 0, 4, ifelse(terr$Model == 4, 5, 
  ifelse(terr$Model == 5, 6, terr$Model)))

terr$Model <- factor(terr$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
terr$Tx <- factor(terr$Tx)
terr$N_Tx <- factor(terr$N_Tx)
terr$TERR1 <- ifelse(terr$beta == 1, terr$TERR1/0.5, ifelse(terr$beta == 2, terr$TERR1/0.75, terr$TERR1/1))
terr$TERR2 <- ifelse(terr$beta == 1, terr$TERR2/0.5, ifelse(terr$beta == 2, terr$TERR2/0.75, terr$TERR2/1))
terr$TERR3 <- ifelse(terr$beta == 1, terr$TERR3/0.75, ifelse(terr$beta == 2, terr$TERR3/1.25, terr$TERR3/1.5))
terr$TERR4 <- ifelse(terr$beta == 1, terr$TERR4/0.75, ifelse(terr$beta == 2, terr$TERR4/1.25, terr$TERR4/1.5))

# graph by beta and N_Tx (Figure 3)
sum_terr1 <- summarySE(terr, measurevar = "TERR1", groupvars = c("beta", "N_Tx", "Model"))
sum_terr2 <- summarySE(terr, measurevar = "TERR2", groupvars = c("beta", "N_Tx", "Model"))
sum_terr3 <- summarySE(terr, measurevar = "TERR3", groupvars = c("beta", "N_Tx", "Model"))
sum_terr4 <- summarySE(terr, measurevar = "TERR4", groupvars = c("beta", "N_Tx", "Model"))

pd <- position_dodge(0.1)
color <- c("#009900", "#CC0000", "#FF9933", "#990099", "#0000FF", "#666666")

beta_names <- list(
  '1'=expression(paste("Small Tx Effect: ", beta[Tx[1]], " = 0.5")),
  '2'=expression(paste("Medium Tx Effect: ", beta[Tx[1]], " = 0.75")),
  '3'=expression(paste("Large Tx Effect: ", beta[Tx[1]], " = 1"))
)
beta_labeller <- function(variable, value){
  return(beta_names[value])
}

terr_plot1 <- ggplot(sum_terr1, 
                     aes(x = N_Tx, y = TERR1, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = TERR1-ci, ymax = TERR1+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~beta, labeller = beta_labeller) + theme_bw() + 
  xlab("") + 
  scale_y_continuous(expression(paste("MAE Proportion for ", beta[Tx[1]]))) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
terr_plot1

beta_names2 <- list(
  '1'=expression(paste("Small Tx Effect: ", beta[Tx[2]], " = -0.5")),
  '2'=expression(paste("Medium Tx Effect: ", beta[Tx[2]], " = -0.75")),
  '3'=expression(paste("Large Tx Effect: ", beta[Tx[2]], " = -1"))
)
beta_labeller2 <- function(variable, value){
  return(beta_names2[value])
}

terr_plot2 <- ggplot(sum_terr2, 
                     aes(x = N_Tx, y = TERR2, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = TERR2-ci, ymax = TERR2+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~beta, labeller = beta_labeller2) + theme_bw() + 
  xlab("") + 
  scale_y_continuous(expression(paste("MAE Proportion for ", beta[Tx[2]]))) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
terr_plot2

beta_names3 <- list(
  '1'=expression(paste("Small Tx Effect: ", beta[Tx[3]], " = 0.75")),
  '2'=expression(paste("Medium Tx Effect: ", beta[Tx[3]], " = 1.25")),
  '3'=expression(paste("Large Tx Effect: ", beta[Tx[3]], " = 1.5"))
)
beta_labeller3 <- function(variable, value){
  return(beta_names3[value])
}

terr_plot3 <- ggplot(sum_terr3, 
                     aes(x = N_Tx, y = TERR3, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = TERR3-ci, ymax = TERR3+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~beta, labeller = beta_labeller3) + theme_bw() + 
  xlab("") + 
  scale_y_continuous(expression(paste("MAE Proportion for ", beta[Tx[3]]))) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
terr_plot3

beta_names4 <- list(
  '1'=expression(paste("Small Tx Effect: ", beta[Tx[4]], " = -0.75")),
  '2'=expression(paste("Medium Tx Effect: ", beta[Tx[4]], " = -1.25")),
  '3'=expression(paste("Large Tx Effect: ", beta[Tx[4]], " = -1.5"))
)
beta_labeller4 <- function(variable, value){
  return(beta_names4[value])
}

terr_plot4 <- ggplot(sum_terr4, 
                     aes(x = N_Tx, y = TERR4, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = TERR4-ci, ymax = TERR4+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~beta, labeller = beta_labeller4) + theme_bw() + 
  xlab("Sample Size per Tx Group") + 
  scale_y_continuous(expression(paste("MAE Proportion for ", beta[Tx[4]]))) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 10), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
terr_plot4

Fig3 <- grid_arrange_shared_legend(rmse_plot1, rmse_plot2, rmse_plot3, rmse_plot4)

ggsave(paste(f4, "Fialkowski_Figure3.png", sep = ""), dpi = 800, Fig3, 
       width = 6.7, height = 8.7, units = "in")


# table by total sample size n (labeled by sample size per group in tables)
sum_terr <- summarySE(terr, measurevar = "TERR1", groupvars = c("beta", "n", "Model"), na.rm = TRUE)
sum_terr <- sum_terr[order(sum_terr$beta, sum_terr$n, sum_terr$Model), ]

sum_terr$CI_1 <- paste(formatC(round(sum_terr$TERR1, 3), 3, format = "f"), " (", 
  formatC(round(sum_terr$TERR1, 3) - round(sum_terr$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_terr$TERR1, 3) + round(sum_terr$ci, 3), 3, format = "f"), ")", sep = "")
sum_terr$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)
sum_terr$n <- rep(c("500", rep(NA, 5), "1250", rep(NA, 5), "2500", rep(NA, 5)), 3)

sum_terr2 <- summarySE(terr, measurevar = "TERR2", groupvars = c("beta", "n", "Model"), na.rm = TRUE)
sum_terr2 <- sum_terr2[order(sum_terr2$beta, sum_terr2$n, sum_terr2$Model), ]

sum_terr2$CI_2 <- paste(formatC(round(sum_terr2$TERR2, 3), 3, format = "f"), " (", 
  formatC(round(sum_terr2$TERR2, 3) - round(sum_terr2$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_terr2$TERR2, 3) + round(sum_terr2$ci, 3), 3, format = "f"), ")", sep = "")

sum_terr4 <- summarySE(terr, measurevar = "TERR3", groupvars = c("beta", "n", "Model"), na.rm = TRUE)
sum_terr4 <- sum_terr4[order(sum_terr4$beta, sum_terr4$n, sum_terr4$Model), ]

sum_terr4$CI_3 <- paste(formatC(round(sum_terr4$TERR3, 3), 3, format = "f"), " (", 
  formatC(round(sum_terr4$TERR3, 3) - round(sum_terr4$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_terr4$TERR3, 3) + round(sum_terr4$ci, 3), 3, format = "f"), ")", sep = "")
sum_terr4$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)
sum_terr4$n <- rep(c("500", rep(NA, 5), "1250", rep(NA, 5), "2500", rep(NA, 5)), 3)

sum_terr5 <- summarySE(terr, measurevar = "TERR4", groupvars = c("beta", "n", "Model"), na.rm = TRUE)
sum_terr5 <- sum_terr5[order(sum_terr5$beta, sum_terr5$n, sum_terr5$Model), ]

sum_terr5$CI_4 <- paste(formatC(round(sum_terr5$TERR4, 3), 3, format = "f"), " (", 
  formatC(round(sum_terr5$TERR4, 3) - round(sum_terr5$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_terr5$TERR4, 3) + round(sum_terr5$ci, 3), 3, format = "f"), ")", sep = "")

sum_terr3 <- cbind(sum_terr[sum_terr$beta == 1, c("Model", "n", "CI_1")], 
  sum_terr2[sum_terr2$beta == 1, "CI_2"], sum_terr4[sum_terr4$beta == 1, "CI_3"], 
  sum_terr5[sum_terr5$beta == 1, "CI_4"])

print(xtable(sum_terr3[, -2], align = c("l", "l", rep("c", 4))), 
      caption.placement = "top", include.rownames = FALSE, 
      include.colnames = FALSE, sanitize.text.function = function(x){x})

sum_terr6 <- cbind(sum_terr[sum_terr$beta == 2, c("Model", "n", "CI_1")], 
  sum_terr2[sum_terr2$beta == 2, "CI_2"], sum_terr4[sum_terr4$beta == 2, "CI_3"], 
  sum_terr5[sum_terr5$beta == 2, "CI_4"])

print(xtable(sum_terr6[, -2], align = c("l", "l", rep("c", 4))), 
      caption.placement = "top", include.rownames = FALSE, 
      include.colnames = FALSE, sanitize.text.function = function(x){x})

sum_terr7 <- cbind(sum_terr[sum_terr$beta == 3, c("Model", "n", "CI_1")], 
  sum_terr2[sum_terr2$beta == 3, "CI_2"], sum_terr4[sum_terr4$beta == 3, "CI_3"], 
  sum_terr5[sum_terr5$beta == 3, "CI_4"])

print(xtable(sum_terr7[, -2], align = c("l", "l", rep("c", 4))), 
      caption.placement = "top", include.rownames = FALSE, 
      include.colnames = FALSE, sanitize.text.function = function(x){x})

# predictive MSE: testing sets
rmse1 <- read.table(paste(f1, "rmse_betas_1.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")
rmse2 <- read.table(paste(f1, "rmse_betas_2.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")
rmse3 <- read.table(paste(f1, "rmse_betas_3.txt", sep = ""), header = TRUE, sep = "\t", 
                   quote = "")

rmse <- rbind(rmse1, rmse2, rmse3)
rmse$beta <- c(rep(1, nrow(rmse)/3), rep(2, nrow(rmse)/3), 
  rep(3, nrow(rmse)/3))

rmse$beta <- factor(rmse$beta)
rmse$n <- factor(rmse$Tx * rmse$N_Tx)
rmse$Model <- ifelse(rmse$Model == 0, 4, ifelse(rmse$Model == 4, 5, 
  ifelse(rmse$Model == 5, 6, rmse$Model)))

rmse$Model <- factor(rmse$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
rmse$Tx <- factor(rmse$Tx)
rmse$N_Tx <- factor(rmse$N_Tx)
rmse$MSE <- rmse$RMSE^2 # use MSE instead
rmse <- rmse[, -which(colnames(rmse) == "RMSE")]

# MSE table by Tx and N_Tx
sum_rmse <- summarySE(rmse, measurevar = "MSE", 
  groupvars = c("beta", "Tx", "N_Tx", "Model"), na.rm = TRUE)
sum_rmse <- sum_rmse[order(sum_rmse$beta, sum_rmse$Tx, sum_rmse$N_Tx, sum_rmse$Model), ]
sum_rmse$CI <- paste(formatC(round(sum_rmse$MSE, 4), 4, format = "f"), " (", 
  formatC(round(sum_rmse$MSE, 4) - round(sum_rmse$ci, 4), 4, format = "f"), 
  ", ", formatC(round(sum_rmse$MSE, 4) + round(sum_rmse$ci, 4), 4, format = "f"), ")", sep = "")

sum_rmse$Tx <- c("2", rep(NA, 17), "3", rep(NA, 17), "5", rep(NA, 17))
sum_rmse$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)

sum_rmse2 <- cbind(sum_rmse[sum_rmse$beta == 1, ], sum_rmse[sum_rmse$beta == 2, ], 
  sum_rmse[sum_rmse$beta == 3, ])
colnames(sum_rmse2)[c(10, 20, 30)] <- c("CI_1", "CI_2", "CI_3")

print(xtable(sum_rmse2[, c("Model", "CI_1", "CI_2", "CI_3")], 
             align = c("l", "l", rep("c", 3))), include.rownames = FALSE,
      sanitize.text.function = function(x){x})

# MSE graph by Tx and N_Tx (Figure 4)
sum_rmse <- summarySE(rmse, measurevar = "MSE", 
  groupvars = c("beta", "Tx", "N_Tx", "Model"), na.rm = TRUE)
sum_rmse <- sum_rmse[order(sum_rmse$beta, sum_rmse$Tx, sum_rmse$N_Tx, sum_rmse$Model), ]

pd <- position_dodge(0.1)
color <- c("#009900", "#CC0000", "#FF9933", "#990099", "#0000FF", "#666666")

tx_names <- list(
  '2'="2 Tx Groups",
  '3'="3 Tx Groups",
  '5'="5 Tx Groups"
)
tx_labeller <- function(variable, value){
  return(tx_names[value])
}

mse_plot1 <- ggplot(sum_rmse[sum_rmse$beta == 1, ], 
                     aes(x = N_Tx, y = MSE, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = MSE-ci, ymax = MSE+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~Tx, labeller = tx_labeller) + theme_bw() + 
  xlab("") + 
  scale_y_continuous("Small Tx Effect", breaks = seq(0.08, 0.15, 0.01), 
                     limits = c(0.08, 0.15)) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
mse_plot1

mse_plot2 <- ggplot(sum_rmse[sum_rmse$beta == 2, ], 
                     aes(x = N_Tx, y = MSE, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = MSE-ci, ymax = MSE+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~Tx, labeller = tx_labeller) + theme_bw() + 
  xlab("") + 
  scale_y_continuous("Medium Tx Effect", breaks = seq(0.08, 0.15, 0.01), 
                     limits = c(0.08, 0.15)) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
mse_plot2

mse_plot3 <- ggplot(sum_rmse[sum_rmse$beta == 3, ], 
                     aes(x = N_Tx, y = MSE, colour = Model, group = Model)) + 
  geom_errorbar(aes(ymin = MSE-ci, ymax = MSE+ci), colour="black", 
                width = .3, size = 1.5, position = pd) +
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~Tx, labeller = tx_labeller) + theme_bw() + 
  xlab("Sample Size per Tx Group") + 
  scale_y_continuous("Large Tx Effect", breaks = seq(0.08, 0.15, 0.01), 
                     limits = c(0.08, 0.15)) +
  scale_colour_manual(name = "", values = color, 
                      labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
                                 "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
mse_plot3

Fig4 <- grid_arrange_shared_legend(mse_plot1, mse_plot2, mse_plot3)

ggsave(paste(f4, "Fialkowski_Figure4.tiff", sep = ""), dpi = 800, Fig4, 
  width = 6.7, height = 8.9, units = "in")

# R^2 table by Tx and N_Tx
sum_rmse <- summarySE(rmse, measurevar = "R2", 
  groupvars = c("beta", "Tx", "N_Tx", "Model"), na.rm = TRUE)
sum_rmse <- sum_rmse[order(sum_rmse$beta, sum_rmse$Tx, sum_rmse$N_Tx, sum_rmse$Model), ]
sum_rmse$CI <- paste(formatC(round(sum_rmse$R2, 3), 3, format = "f"), " (", 
  formatC(round(sum_rmse$R2, 3) - round(sum_rmse$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_rmse$R2, 3) + round(sum_rmse$ci, 3), 3, format = "f"), ")", sep = "")

sum_rmse$Tx <- c("2", rep(NA, 17), "3", rep(NA, 17), "5", rep(NA, 17))
sum_rmse$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)

sum_rmse2 <- cbind(sum_rmse[sum_rmse$beta == 1, ], sum_rmse[sum_rmse$beta == 2, ], 
  sum_rmse[sum_rmse$beta == 3, ])
colnames(sum_rmse2)[c(10, 20, 30)] <- c("CI_1", "CI_2", "CI_3")

print(xtable(sum_rmse2[, c("Model", "CI_1", "CI_2", "CI_3")], 
             align = c("l", "l", rep("c", 3))), include.rownames = FALSE,
      sanitize.text.function = function(x){x})

# MSE: training sets
rmse1 <- read.table(paste(f1, "train_rmse_betas_1.txt", sep = ""), header = TRUE, sep = "\t", 
                    quote = "")
rmse2 <- read.table(paste(f1, "train_rmse_betas_2.txt", sep = ""), header = TRUE, sep = "\t", 
                    quote = "")
rmse3 <- read.table(paste(f1, "train_rmse_betas_3.txt", sep = ""), header = TRUE, sep = "\t", 
                    quote = "")

rmse <- rbind(rmse1, rmse2, rmse3)
rmse$beta <- c(rep(1, nrow(rmse)/3), rep(2, nrow(rmse)/3), 
               rep(3, nrow(rmse)/3))

rmse$beta <- factor(rmse$beta)
rmse$n <- factor(rmse$Tx * rmse$N_Tx)
rmse$Model <- ifelse(rmse$Model == 0, 4, ifelse(rmse$Model == 4, 5, 
                                                ifelse(rmse$Model == 5, 6, rmse$Model)))

rmse$Model <- factor(rmse$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
rmse$Tx <- factor(rmse$Tx)
rmse$N_Tx <- factor(rmse$N_Tx)
rmse$MSE <- rmse$RMSE^2 # use MSE instead
rmse <- rmse[, -which(colnames(rmse) == "RMSE")]

# table by Tx and N_Tx
sum_rmse <- summarySE(rmse, measurevar = "MSE", 
  groupvars = c("beta", "Tx", "N_Tx", "Model"), na.rm = TRUE)
sum_rmse <- sum_rmse[order(sum_rmse$beta, sum_rmse$Tx, sum_rmse$N_Tx, sum_rmse$Model), ]
sum_rmse$CI <- paste(formatC(round(sum_rmse$MSE, 3), 3, format = "f"), " (", 
  formatC(round(sum_rmse$MSE, 3) - round(sum_rmse$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_rmse$MSE, 3) + round(sum_rmse$ci, 3), 3, format = "f"), ")", sep = "")

sum_rmse$Tx <- c("2", rep(NA, 17), "3", rep(NA, 17), "5", rep(NA, 17))
sum_rmse$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)

sum_rmse2 <- cbind(sum_rmse[sum_rmse$beta == 1, ], sum_rmse[sum_rmse$beta == 2, ], 
  sum_rmse[sum_rmse$beta == 3, ])
colnames(sum_rmse2)[c(10, 20, 30)] <- c("CI_1", "CI_2", "CI_3")

print(xtable(sum_rmse2[, c("Model", "CI_1", "CI_2", "CI_3")], 
             align = c("l", "l", rep("c", 3))), include.rownames = FALSE,
      sanitize.text.function = function(x){x})

# R^2 table by Tx and N_Tx
sum_rmse <- summarySE(rmse, measurevar = "R2", 
  groupvars = c("beta", "Tx", "N_Tx", "Model"), na.rm = TRUE)
sum_rmse <- sum_rmse[order(sum_rmse$beta, sum_rmse$Tx, sum_rmse$N_Tx, sum_rmse$Model), ]
sum_rmse$CI <- paste(formatC(round(sum_rmse$R2, 3), 3, format = "f"), " (", 
  formatC(round(sum_rmse$R2, 3) - round(sum_rmse$ci, 3), 3, format = "f"), 
  ", ", formatC(round(sum_rmse$R2, 3) + round(sum_rmse$ci, 3), 3, format = "f"), ")", sep = "")

sum_rmse$Tx <- c("2", rep(NA, 17), "3", rep(NA, 17), "5", rep(NA, 17))
sum_rmse$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)

sum_rmse2 <- cbind(sum_rmse[sum_rmse$beta == 1, ], sum_rmse[sum_rmse$beta == 2, ], 
  sum_rmse[sum_rmse$beta == 3, ])
colnames(sum_rmse2)[c(10, 20, 30)] <- c("CI_1", "CI_2", "CI_3")

print(xtable(sum_rmse2[, c("Model", "CI_1", "CI_2", "CI_3")], 
             align = c("l", "l", rep("c", 3))), include.rownames = FALSE,
      sanitize.text.function = function(x){x})

# False and true positive rates

# Percentile CI
# beta_1
pr_2 <- read.table(paste(f1, "nperc_pr_2_betas_1.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_3 <- read.table(paste(f1, "nperc_pr_3_betas_1.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_5 <- read.table(paste(f1, "nperc_pr_5_betas_1.txt", sep = ""), header = TRUE, sep = "\t", quote = "")

pr1 <- as.data.frame(rbind(pr_2, pr_3, pr_5))
pr1$Model <- ifelse(pr1$Model == 0, 4, ifelse(pr1$Model == 4, 5, 
  ifelse(pr1$Model == 5, 6, pr1$Model)))

pr1$Model <- factor(pr1$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
pr1$Tx <- factor(pr1$Tx)
pr1$N_Tx <- factor(pr1$N_Tx)

f_pr1 <- summarySE(pr1, measurevar = "fpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
f_pr1 <- f_pr1[order(f_pr1$Tx, f_pr1$N_Tx, f_pr1$Model), ]

t_pr1 <- summarySE(pr1, measurevar = "tpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
t_pr1 <- t_pr1[order(t_pr1$Tx, t_pr1$N_Tx, t_pr1$Model), ]

# beta_2
pr_2 <- read.table(paste(f1, "nperc_pr_2_betas_2.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_3 <- read.table(paste(f1, "nperc_pr_3_betas_2.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_5 <- read.table(paste(f1, "nperc_pr_5_betas_2.txt", sep = ""), header = TRUE, sep = "\t", quote = "")

pr2 <- as.data.frame(rbind(pr_2, pr_3, pr_5))
pr2$Model <- ifelse(pr2$Model == 0, 4, ifelse(pr2$Model == 4, 5, 
  ifelse(pr2$Model == 5, 6, pr2$Model)))

pr2$Model <- factor(pr2$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
pr2$Tx <- factor(pr2$Tx)
pr2$N_Tx <- factor(pr2$N_Tx)

f_pr2 <- summarySE(pr2, measurevar = "fpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
f_pr2 <- f_pr2[order(f_pr2$Tx, f_pr2$N_Tx, f_pr2$Model), ]

t_pr2 <- summarySE(pr2, measurevar = "tpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
t_pr2 <- t_pr2[order(t_pr2$Tx, t_pr2$N_Tx, t_pr2$Model), ]

# beta_3
pr_2 <- read.table(paste(f1, "nperc_pr_2_betas_3.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_3 <- read.table(paste(f1, "nperc_pr_3_betas_3.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_5 <- read.table(paste(f1, "nperc_pr_5_betas_3.txt", sep = ""), header = TRUE, sep = "\t", quote = "")

pr3 <- as.data.frame(rbind(pr_2, pr_3, pr_5))
pr3$Model <- ifelse(pr3$Model == 0, 4, ifelse(pr3$Model == 4, 5, 
  ifelse(pr3$Model == 5, 6, pr3$Model)))

pr3$Model <- factor(pr3$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
pr3$Tx <- factor(pr3$Tx)
pr3$N_Tx <- factor(pr3$N_Tx)

f_pr3 <- summarySE(pr3, measurevar = "fpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
f_pr3 <- f_pr3[order(f_pr3$Tx, f_pr3$N_Tx, f_pr3$Model), ]

t_pr3 <- summarySE(pr3, measurevar = "tpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
t_pr3 <- t_pr3[order(t_pr3$Tx, t_pr3$N_Tx, t_pr3$Model), ]

f_pr1$CI_1 <- paste(formatC(round(f_pr1$fpr2, 3), 3, format = "e"), " (", 
                   formatC(round(f_pr1$fpr2, 3) - round(f_pr1$ci, 3), 3, format = "e"), 
                   ", ", formatC(round(f_pr1$fpr2, 3) + round(f_pr1$ci, 3), 3, format = "e"), ")", sep = "")
f_pr2$CI_2 <- paste(formatC(round(f_pr2$fpr2, 3), 3, format = "e"), " (", 
                    formatC(round(f_pr2$fpr2, 3) - round(f_pr2$ci, 3), 3, format = "e"), 
                    ", ", formatC(round(f_pr2$fpr2, 3) + round(f_pr2$ci, 3), 3, format = "e"), ")", sep = "")
f_pr3$CI_3 <- paste(formatC(round(f_pr3$fpr2, 3), 3, format = "e"), " (", 
                    formatC(round(f_pr3$fpr2, 3) - round(f_pr3$ci, 3), 3, format = "e"), 
                    ", ", formatC(round(f_pr3$fpr2, 3) + round(f_pr3$ci, 3), 3, format = "e"), ")", sep = "")

fpr <- cbind(f_pr1[, c("Model", "Tx", "N_Tx", "fpr2")], f_pr2$fpr2, f_pr3$fpr2)
fpr$Tx <- c("2", rep(NA, 17), "3", rep(NA, 17), "5", rep(NA, 17))
fpr$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)

print(xtable(cbind(fpr[1:18, -c(2, 3)], fpr[19:36, -c(1, 2, 3)], fpr[37:nrow(fpr), -c(1, 2, 3)]), 
  align = c("l", "l", "c", "c", "c", "c", "c", "c", "c", "c", "c"), 
  digits = c(0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5)), include.rownames = FALSE,
  sanitize.text.function = function(x){x})

tpr <- cbind(t_pr1[, c("Model", "Tx", "N_Tx", "tpr2")], t_pr2$tpr2, t_pr3$tpr2)
tpr$Tx <- c("2", rep(NA, 17), "3", rep(NA, 17), "5", rep(NA, 17))
tpr$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)

print(xtable(cbind(tpr[1:18, -c(2, 3)], tpr[19:36, -c(1, 2, 3)], tpr[37:nrow(tpr), -c(1, 2, 3)]), 
             align = c("l", "l", "c", "c", "c", "c", "c", "c", "c", "c", "c"), 
             digits = c(0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5)), include.rownames = FALSE,
      sanitize.text.function = function(x){x})

# Figure 5
pd <- position_dodge(0.1)
color <- c("#009900", "#CC0000", "#FF9933", "#990099", "#0000FF", "#666666")

tx_names <- list(
  '2'="2 Tx Groups",
  '3'="3 Tx Groups",
  '5'="5 Tx Groups"
)
tx_labeller <- function(variable, value){
  return(tx_names[value])
}

t_plot1 <- ggplot(t_pr1, 
  aes(x = N_Tx, y = tpr2, colour = Model, group = Model)) + 
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~Tx, labeller = tx_labeller) + theme_bw() + 
  scale_y_continuous("Small Tx Effect", breaks = seq(0, 1, 0.1)) +
  xlab("") +
  scale_colour_manual(name = "", values = color, 
    labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
               "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
t_plot1

t_plot2 <- ggplot(t_pr2, 
  aes(x = N_Tx, y = tpr2, colour = Model, group = Model)) + 
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~Tx, labeller = tx_labeller) + theme_bw() + 
  scale_y_continuous("Medium Tx Effect", breaks = seq(0, 1, 0.1)) +
  xlab("") +
  scale_colour_manual(name = "", values = color, 
    labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
               "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
t_plot2

t_plot3 <- ggplot(t_pr3, 
  aes(x = N_Tx, y = tpr2, colour = Model, group = Model)) + 
  geom_line(position = pd, lwd = 1) + 
  geom_point(position = pd, size = 2, shape = 21, fill = "white") + 
  facet_wrap(~Tx, labeller = tx_labeller) + theme_bw() + 
  scale_y_continuous("Large Tx Effect", breaks = seq(0, 1, 0.1)) +
  xlab("Sample Size per Tx Group") +
  scale_colour_manual(name = "", values = color, 
    labels = c("EN25  ", "EN50  ", "EN75  ", "Lasso  ", 
               "BLasso  ", "MARS  ")) + #l = 40)
  theme(axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10), 
        strip.text = element_text(size = 12), 
        legend.text = element_text(size = 10), 
        legend.position = "none")
t_plot3

Fig5 <- grid_arrange_shared_legend(t_plot1, t_plot2, t_plot3)

ggsave(paste(f4, "Fialkowski_Figure5.png", sep = ""), dpi = 800, Fig5, 
       width = 6.7, height = 8.7, units = "in")

# Normal CI
# beta_1
pr_2 <- read.table(paste(f1, "npr_2_betas_1.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_3 <- read.table(paste(f1, "npr_3_betas_1.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_5 <- read.table(paste(f1, "npr_5_betas_1.txt", sep = ""), header = TRUE, sep = "\t", quote = "")

pr1 <- as.data.frame(rbind(pr_2, pr_3, pr_5))
pr1$Model <- ifelse(pr1$Model == 0, 4, ifelse(pr1$Model == 4, 5, 
  ifelse(pr1$Model == 5, 6, pr1$Model)))

pr1$Model <- factor(pr1$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
pr1$Tx <- factor(pr1$Tx)
pr1$N_Tx <- factor(pr1$N_Tx)

f_pr1 <- summarySE(pr1, measurevar = "fpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
f_pr1 <- f_pr1[order(f_pr1$Tx, f_pr1$N_Tx, f_pr1$Model), ]

t_pr1 <- summarySE(pr1, measurevar = "tpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
t_pr1 <- t_pr1[order(t_pr1$Tx, t_pr1$N_Tx, t_pr1$Model), ]

# beta_2
pr_2 <- read.table(paste(f1, "npr_2_betas_2.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_3 <- read.table(paste(f1, "npr_3_betas_2.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_5 <- read.table(paste(f1, "npr_5_betas_2.txt", sep = ""), header = TRUE, sep = "\t", quote = "")

pr2 <- as.data.frame(rbind(pr_2, pr_3, pr_5))
pr2$Model <- ifelse(pr2$Model == 0, 4, ifelse(pr2$Model == 4, 5, 
  ifelse(pr2$Model == 5, 6, pr2$Model)))

pr2$Model <- factor(pr2$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
pr2$Tx <- factor(pr2$Tx)
pr2$N_Tx <- factor(pr2$N_Tx)

f_pr2 <- summarySE(pr2, measurevar = "fpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
f_pr2 <- f_pr2[order(f_pr2$Tx, f_pr2$N_Tx, f_pr2$Model), ]

t_pr2 <- summarySE(pr2, measurevar = "tpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
t_pr2 <- t_pr2[order(t_pr2$Tx, t_pr2$N_Tx, t_pr2$Model), ]

# beta_3
pr_2 <- read.table(paste(f1, "npr_2_betas_3.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_3 <- read.table(paste(f1, "npr_3_betas_3.txt", sep = ""), header = TRUE, sep = "\t", quote = "")
pr_5 <- read.table(paste(f1, "npr_5_betas_3.txt", sep = ""), header = TRUE, sep = "\t", quote = "")

pr3 <- as.data.frame(rbind(pr_2, pr_3, pr_5))
pr3$Model <- ifelse(pr3$Model == 0, 4, ifelse(pr3$Model == 4, 5, 
  ifelse(pr3$Model == 5, 6, pr3$Model)))

pr3$Model <- factor(pr3$Model, labels = c("EN25", "EN50", "EN75", "Lasso", "BLasso", "MARS"))
pr3$Tx <- factor(pr3$Tx)
pr3$N_Tx <- factor(pr3$N_Tx)

f_pr3 <- summarySE(pr3, measurevar = "fpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
f_pr3 <- f_pr3[order(f_pr3$Tx, f_pr3$N_Tx, f_pr3$Model), ]

t_pr3 <- summarySE(pr3, measurevar = "tpr2", groupvars = c("Model", "Tx", "N_Tx"), na.rm = TRUE)
t_pr3 <- t_pr3[order(t_pr3$Tx, t_pr3$N_Tx, t_pr3$Model), ]

f_pr1$CI_1 <- paste(formatC(round(f_pr1$fpr2, 3), 3, format = "e"), " (", 
                   formatC(round(f_pr1$fpr2, 3) - round(f_pr1$ci, 3), 3, format = "e"), 
                   ", ", formatC(round(f_pr1$fpr2, 3) + round(f_pr1$ci, 3), 3, format = "e"), ")", sep = "")
f_pr2$CI_2 <- paste(formatC(round(f_pr2$fpr2, 3), 3, format = "e"), " (", 
                    formatC(round(f_pr2$fpr2, 3) - round(f_pr2$ci, 3), 3, format = "e"), 
                    ", ", formatC(round(f_pr2$fpr2, 3) + round(f_pr2$ci, 3), 3, format = "e"), ")", sep = "")
f_pr3$CI_3 <- paste(formatC(round(f_pr3$fpr2, 3), 3, format = "e"), " (", 
                    formatC(round(f_pr3$fpr2, 3) - round(f_pr3$ci, 3), 3, format = "e"), 
                    ", ", formatC(round(f_pr3$fpr2, 3) + round(f_pr3$ci, 3), 3, format = "e"), ")", sep = "")

fpr <- cbind(f_pr1[, c("Model", "Tx", "N_Tx", "fpr2")], f_pr2$fpr2, f_pr3$fpr2)
fpr$Tx <- c("2", rep(NA, 17), "3", rep(NA, 17), "5", rep(NA, 17))
fpr$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)

print(xtable(cbind(fpr[1:18, -c(2, 3)], fpr[19:36, -c(1, 2, 3)], fpr[37:nrow(fpr), -c(1, 2, 3)]), 
  align = c("l", "l", "c", "c", "c", "c", "c", "c", "c", "c", "c"), 
  digits = c(0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5)), include.rownames = FALSE,
  sanitize.text.function = function(x){x})

tpr <- cbind(t_pr1[, c("Model", "Tx", "N_Tx", "tpr2")], t_pr2$tpr2, t_pr3$tpr2)
tpr$Tx <- c("2", rep(NA, 17), "3", rep(NA, 17), "5", rep(NA, 17))
tpr$N_Tx <- rep(c("100", rep(NA, 5), "250", rep(NA, 5), "500", rep(NA, 5)), 3)

print(xtable(cbind(tpr[1:18, -c(2, 3)], tpr[19:36, -c(1, 2, 3)], tpr[37:nrow(tpr), -c(1, 2, 3)]), 
             align = c("l", "l", "c", "c", "c", "c", "c", "c", "c", "c", "c"), 
             digits = c(0, 0, 5, 5, 5, 5, 5, 5, 5, 5, 5)), include.rownames = FALSE,
      sanitize.text.function = function(x){x})
