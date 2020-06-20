rm(list=ls())

path <- "/nfs/general/repos/direct-approach-to-fdr/"

set.seed(1)

library(xtable)
library(latex2exp)

get_cor <- function(ro) {
  return(ro^2 / (1 + ro^2))
}

m <- 1000 # number of hypotheses
lambda <- 0.5 # fix value of tuning parameter lambda
ro <- 0  # 0.6
gamma <- .01 # rejection region
pi0_range <- seq(0.1, 0.9, 0.1) # proportion of null hypotheses range
m0_range <- round(pi0_range * m)

figure_size = 5

###########################################################################
###########################################################################
### Average detection power comparison of BH procedure and direct approach

# BH procedure and improved BH procedure where we plug hat_pi_0 in BH procedure
for (use_estimate in c(FALSE, TRUE)) {
  
  # Table 2, Table 3
  table_1_names <- c("FDR_S", "FDR_BH", "Power_S", "Power_BH", "FDR_hat_S", "pi_0_hat_S", "gamma_hat_BH")
  table_1 <- matrix(0, length(pi0_range), length(table_1_names))
  colnames(table_1) <- table_1_names
  rownames(table_1) <- pi0_range
  for (m0 in m0_range) {
      # test statistic
      Z <- cbind(matrix(rnorm(m0 * m, 0, 1), ncol = m0), # pi_0 * m under H_0: N(0, 1)
                 matrix(rnorm((m - m0) * m, 2, 1), ncol = (m - m0))) # (1 - pi_0) * m under H_1: N(2, 1)
      # we simulate dependence by adding the same vector M to all the Z's
      M <- matrix(rep(rnorm(m, 0, ro), m), ncol = m)
      Z <- Z + M
      P <- pnorm(Z, lower.tail = FALSE) # p-values
  
      #############################################
      # direct approach
      # estimate pi_0 using Storey's method and Pr(P \leq \gamma)
      W <- apply(P > lambda, 1, sum)
      R <- apply(P <= gamma, 1, sum)
      Pi_0_hat <- W/( (1 - lambda) * m )
      pi_0_hat <- mean(Pi_0_hat) # average pi_0 estimate using Storey's method
  
      PrP_hat <- R/m
      PrP_hat[PrP_hat == 0] <- 1/m # because of R(\gamma) v 1 in estimation of \hat{Pr(P \leq \gamma)}
  
      # estimate of FDR_lambda(gamma)
      FDR_hat <- (Pi_0_hat * gamma) / PrP_hat
      mean_FDR_hat <- mean(FDR_hat) # average FDR esstimate using Storey's MTP
  
  
      V <- apply(P[, 1:m0] <= gamma, 1, sum) # number of type I errors
  
      PrR_S <- (1 - sum(R==0)/m)
      if (PrR_S != 1) {
        cat(m0)
        cat("\n")
        cat("PrR_S = ", PrR_S) # some R=0 when using Storey's estimators in direct MTP
        cat("\n")
      }
      FDR_S <- mean(V[R!=0]/R[R!=0]) * PrR_S
  
      S <- apply(P[, (m0 + 1):m] <= gamma, 1, sum) # number of true positives
      Power_S <- mean(S/(m - m0)) # average power
  
      table_1[as.character(m0/m), "FDR_S"] <- FDR_S
      table_1[as.character(m0/m), "Power_S"] <- Power_S
      table_1[as.character(m0/m), "FDR_hat_S"] <- mean_FDR_hat
      table_1[as.character(m0/m), "pi_0_hat_S"] <- pi_0_hat
  
      #############################################
      # B-H procedure
      alphas <- FDR_hat # we want FDR <= alpha
      RHS <- t(matrix(rep(1:m, m), m, m))
      if (use_estimate == TRUE) {
        RHS <- RHS * (alphas/(Pi_0_hat*m))
      } else {
        RHS <- RHS * (alphas/m)
      }
  
      P_sort <- t(apply(P, 1, sort))
      tmp <- P_sort <= RHS
      K_kappas <- apply(tmp, 1, sum)
  
      # rejecting p_(1), ..., p_(k_hat) provides FDR <= alpha
      # estimate rejection region [0, gamma]
      Gamma_hat <- sapply(1:m, function(i) ifelse(K_kappas[i] == 0, 0, P_sort[i, K_kappas[i]]))
  
      R <- apply(P <= Gamma_hat, 1, sum);
      V <- apply(P[, 1:m0] <= Gamma_hat, 1, sum) # V = #{FP}
      PrR_BH <- (1 - sum(R==0)/m)
      if (PrR_BH != 1) {
        cat(m0)
        cat("\n")
        cat("PrR_BH = ", PrR_BH) # some R = 0 when using BH MTP
        cat("\n")
      }
  
      FDR_BH <- mean(V[R!=0]/R[R!=0]) * PrR_BH # TODO: report also pFDR or the ratio between them?
  
      S <- apply(P[, (m0 + 1):m] <= Gamma_hat, 1, sum) # S = #{TP}
      Power_BH <- mean(S/(m - m0))
  
      table_1[as.character(m0/m), "FDR_BH"] <- FDR_BH
      table_1[as.character(m0/m), "Power_BH"] <- Power_BH
      table_1[as.character(m0/m), "gamma_hat_BH"] <- mean(unlist(Gamma_hat))
  }
  print(xtable(round(table_1, 3), digits = c(2, rep(3, 6), 4)))
  
  # Figure 1, Figure 2
  if (use_estimate == TRUE) {
    pdf(paste(path, "power_pi0_BH", ".pdf", sep=""), height=figure_size, width=figure_size)
  } else {
    pdf(paste(path, "power", ".pdf", sep=""), height=figure_size, width=figure_size)
  }
  PS <- as.numeric(table_1[, 3])
  PBH <- as.numeric(table_1[, 4])
  p_range <- range(0, PS, PBH+0.1)
  plot(pi0_range, PBH, type = "o", col="blue", ylim = p_range,
       axes = FALSE, xlab = TeX("$\\pi_{0}$"), ylab="Average power")
  axis(1, at=pi0_range)
  axis(2, at=c(0, 0.1, 0.2, 0.3, 0.4))
  box()
  lines(pi0_range, PS, type = "o", pch=22, lty=2, col="red")
  legend(0.7, 0.1, c("BH", "Storey"), cex=0.8,
         col = c("blue", "red"), pch=21:22, lty=1:2)
  dev.off()
}

##########################################################################################
##########################################################################################
## Properties of Storey’s estimator of \hat{pi_{0}} as a function of \lambda
# Figure 3

pdf(paste(path, "pi0_of_lambda", ".pdf", sep=""), height=figure_size, width=figure_size)
m0 <- 900
ro <- 0
lambda_range <- seq(0, 0.99999, 0.001)
pi0_hats <- c()
for (lambda in lambda_range) {
  Z <- cbind(matrix(rnorm(m0, 0, 1), ncol = m0), # pi_0 * m under H_0: N(0, 1)
             matrix(rnorm((m - m0), 2, 1), ncol = (m - m0))) # (1 - pi_0) * m under H_1: N(2, 1)
  P <- pnorm(Z, lower.tail = FALSE) # p-values

  # estimate pi_0 using Storey's method and Pr(P \leq \gamma)
  W <- apply(P > lambda, 1, sum)
  R <- apply(P <= gamma, 1, sum)
  Pi_0_hat <- W/( (1 - lambda) * m )
  pi0_hats <- c(pi0_hats, min(1, Pi_0_hat))
}

plot(lambda_range, pi0_hats, type = "l", xlab = TeX("$\\lambda$"),
     ylab = TeX("$\\widehat{\\pi}_{0}$"))
dev.off()


##########################################################################################
# estimates of pi_0 for various mu_1 and lambda over pi_0 range
# Table 4
table_3_names <- c("2;0.5", "2;0.9", "1;0.5", "1;0.9")
table_3 <- matrix(0, length(pi0_range), length(table_3_names))
colnames(table_3) <- table_3_names
rownames(table_3) <- pi0_range

for (m0 in m0_range) {
  for (mu2 in c(2, 1)) {
    for (lambda in c(.5, 0.9)) {
      # test statistics
      Z <- cbind(matrix(rnorm(m0 * m, 0, 1), ncol = m0), # pi_0 * m under H_0: N(0, 1)
                 matrix(rnorm((m - m0) * m, mu2, 1), ncol = (m - m0))) # (1 - pi_0) * m under H_1: N(2, 1)
      P <- pnorm(Z, lower.tail = FALSE) # p-values

      # estimate pi_0 using Storey's method
      W <- apply(P > lambda, 1, sum)
      Pi_0_hat <- W/( (1 - lambda) * m )
      table_3[as.character(m0/m), paste(mu2, lambda, sep=";")] <- mean(Pi_0_hat) # average pi_0 estimate using Storey's method
    }
  }
}
print(xtable(table_3, digits = 3))


##########################################################################################
# Estimates of pi_0 using Storey's bootstrap method
# Table 5
on_HPC = FALSE # TRUE if running on fast CPU
if (on_HPC) {
  table_5_names <- c("pi_0", "E(pi_0_hat)")
  table_5 <- matrix(0, length(pi0_range), length(table_5_names))
  colnames(table_5) <- table_5_names
  rownames(table_5) <- pi0_range
  
  Rl <- seq(0, 0.95, 0.05) # lambda range for bootrap method
  B <- 10 # 1000
  mu2 <- 2
  
  for (m0 in m0_range) {
    # test statistics
    Z <- cbind(matrix(rnorm(m0 * m, 0, 1), ncol = m0), # pi_0 * m under H_0: N(0, 1)
               matrix(rnorm((m - m0) * m, mu2, 1), ncol = (m - m0))) # (1 - pi_0) * m under H_1: N(2, 1)
    P <- pnorm(Z, lower.tail = FALSE) # p-values
    # in P[i, ] are observed p-values of the first i'th simulation
    
    # first calculate min{\hat{\pi_{0}}}
    Pi_0_hat_of_lambda <- matrix(NA, nrow = 1000, ncol = length(Rl))
    for (j in 1:length(Rl)) {
      # estimate pi_0 using Storey's method
      W <- apply(P > Rl[j], 1, sum)
      Pi_0_hat <- W/( (1 - Rl[j]) * m )
      Pi_0_hat_of_lambda[, j] <- Pi_0_hat
    }
    Pi_0_hat_min <- apply(Pi_0_hat_of_lambda, 1, min) # min estimates
    
    MSE_hats <- matrix(NA, nrow = 1000, ncol = length(Rl))
    for (j in 1:length(Rl)) {
      
      # estimate pi_0(lambda) on bootstrap samples
      Pi_0_hat_bootstrap <- matrix(NA, nrow = 1000, ncol = B)
      for (b in 1:B) {
        Pb <- apply(P, 1, function(x) { sample(x, 1000, TRUE) })
        W <- apply(Pb > Rl[j], 1, sum)
        Pi_0_hat <- W/( (1 - Rl[j]) * m )
        Pi_0_hat_bootstrap[, b] <- Pi_0_hat #[hatpi0_1, ..., hatpi0_B] in each line in 1...1000
      }
      
      # estimate MSE(lambda)
      MSE_hat <- sapply(1:1000, function(i) {
        sum((Pi_0_hat_bootstrap[i, ] - Pi_0_hat_min[i])^2)/B
      })
      MSE_hats[, j] <- MSE_hat
    }
    
    lambda_hat <- sapply(1:1000, function(i) { Rl[which.min(MSE_hats[i, ])] })
    # finally calculate pi_0_hat(lambda_hat)
    W <- apply(P > lambda_hat, 1, sum)
    Pi_0_hat <- W/( (1 - lambda_hat) * m )
    table_5[as.character(m0/m), table_5_names[2]] <- mean(Pi_0_hat)
  }
  saveRDS(table_5, paste(path, "table_5", ".rds", sep=""))
} else {
  table_5 = readRDS(paste(path, "table_5", ".rds", sep=""))  
}
print(xtable(table_5, digits = 3))


##########################################################################################
##########################################################################################
### Case of dependence

pdf(paste(path, "FDR_of_ro", ".pdf", sep=""), height=figure_size, width=figure_size)
m0 <- 900
ro_range <- seq(0, 0.9, 0.05)
errors <- c()
lambda <- 0.5
for (ro in ro_range) {
  Z <- cbind(matrix(rnorm(m0 * m, 0, 1), ncol = m0), # pi_0 * m under H_0: N(0, 1)
             matrix(rnorm((m - m0) * m, 2, 1), ncol = (m - m0))) # (1 - pi_0) * m under H_1: N(2, 1)
  # simulate dependence by adding the same vector to all the variables
  M <- matrix(rep(rnorm(m, 0, ro), m), ncol = m)
  Z <- Z + M
  P <- pnorm(Z, lower.tail = FALSE) # p-values

  #############################################
  # direct approach
  # estimate pi_0 using Storey's method and Pr(P \leq \gamma)
  W <- apply(P > lambda, 1, sum)
  R <- apply(P <= gamma, 1, sum)
  Pi_0_hat <- W/( (1 - lambda) * m )
  pi_0_hat <- mean(Pi_0_hat) # average pi_0 estimate using Storey's method

  PrP_hat <- R/m
  PrP_hat[PrP_hat == 0] <- 1/m # because of R(\gamma) v 1 in estimation of \hat{Pr(P \leq \gamma)}

  # estimate of FDR_lambda(gamma)
  FDR_hat <- (Pi_0_hat * gamma) / PrP_hat
  sum(FDR_hat > 1)
  mean_FDR_hat <- mean(FDR_hat) # average FDR esstimate using Storey's MTP


  V <- apply(P[, 1:m0] <= gamma, 1, sum) # number of type I errors

  PrR_S <- (1 - sum(R==0)/m)
  if (PrR_S != 1) {
    cat(m0)
    cat("\n")
    cat("PrR_S = ", PrR_S) # some R=0 when using Storey's estimators in direct MTP
    cat("\n")
  }
  FDR_S <- mean(V[R!=0]/R[R!=0]) * PrR_S

  errors <- c(errors, mean_FDR_hat - FDR_S)
}
par(mar=c(5,6,4,2)+0.1)
plot(ro_range, errors, type = "l", xlab = TeX("$\\rho$"),
     ylab = TeX("$E(\\widehat{FDR}_{ST}) - FDR_{ST}$"))
dev.off()
