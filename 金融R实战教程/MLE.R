#############################
##########- 实验3 -##########
#############################

# 学习如何使用R进行极大似然估计

# 背景介绍：
# 极大似然估计（Maximum Likelihood Estimation，MLE）是用来估计一个概率模型的参数的一种方法
# 其核心思想是利用已知的样本结果信息，反推最具有可能（最大概率）导致这些样本结果出现的模型参数值
# 参考Sheppard（2013），Greene（2002）
# 本代码使用MLE进行最小二乘估计，并做三大检验，具体而言是实现Sheppard（2013）在2.7.3提到的内容
# 即对r - r_f = beta0 + beta1 * (r_m - r_f) + beta2 * r_s + beta3 * r_v + epsilon使用MLE进行估计

# 实验过程：
# 1. 读取收益率数据，得到估计的原始数据
# 2. 进行极大似然估计
# 3. 对估计的参数进行三大检验

# 实验目标：
# 1. 理解极大似然估计的思想和方法，以及如何使用R实现
# 2. 掌握使用R进行统计检验
# 3. 进一步熟练如何定义函数

# 预期结果：
# 估计的极大似然参数值，以及检验的统计量

# 作业
# 将GLD和GOOG两只证券的MLE参数估计出来


##########- 路径设置 -##########
setwd("~/../OneDrive/Works/金融科技兴趣小组/项目/Lab3_极大似然估计/")
library(bbmle)

##########- 定义函数 -##########
LL <- function(beta0, beta1, beta2, beta3, mu, sigma) {
  # 本函数得到P125中的第二部分方程：l(r|X;theta)=···
  # 先计算出残差，即：r-beta'*x
  rsd <- stk_rtn - MKT * beta1 - SMB * beta2 - HML * beta3 - beta0
  # 对数似然函数，dnorm得到正态分布的概率密度函数
  R <- suppressWarnings(dnorm(rsd, mu, sigma))
  -sum(log(R))
}

LL_H0 <- function(beta0, beta1, mu, sigma) {
  rsd <- stk_rtn - MKT * beta1 - beta0
  R <- suppressWarnings(dnorm(rsd, mu, sigma))
  -sum(log(R))
}


wald_test = function(L, thetahat, Vn, h=0) {
  # H0: L theta = h
  # Note Vn is the asymptotic covariance matrix, so it's the
  # Consistent estimator divided by n. For true Wald tests
  # based on numerical MLEs, just use the inverse of the Hessian.
  
  WaldTest <- numeric(3)
  names(WaldTest) <- c("W", "df", "p-value")
  r <- dim(L)[1]
  W <- t(L %*% thetahat - h) %*% solve(L %*% Vn %*% t(L)) %*% (L %*% thetahat - h)
  W <- as.numeric(W)
  pval <- 1 - pchisq(W, r)
  
  WaldTest[1] <- W
  WaldTest[2] <- r
  WaldTest[3] <- pval
  return(WaldTest)
}


##########- 数据处理 -##########
raw_data <- read.csv("./data.csv", stringsAsFactors=F, na.strings="NA", check.names=F)
stock <- "XOM"
data_df <- raw_data[!is.na(raw_data[, stock]), c("Mkt-RF", "SMB", "HML", "RF", stock)]

stk_rtn <- data_df[, stock] * 100 - data_df[, "RF"]
MKT <- data_df[, "Mkt-RF"]
SMB <- data_df[, "SMB"]
HML <- data_df[, "HML"]

##########- 参数估计 -##########
# 下面提到的页数均指Sheppard（2013）中的页码
# 使用bbmle包里的mle2函数进行估计
# MLE估计的核心是写出极大似然函数，然后对极大似然函数进行最大化
# 因为假设回归的残差服从正态分布，因此可以对残差写极大似然函数
# 下面是#P125中下半部分的估计
fit <- mle2(LL, start = list(beta0 = 0, beta1 = 1, beta2=0, beta3=0, mu=0, sigma=1), fixed = list(mu = 0))
# 估计结果对应P131第一个表格，数值有一些差别，但基本一致
summary(fit)
print(coef(fit))

# 有约束的估计，为了后续的检验做准备
# 这里的H0: beta2 = beta3 = 0，因此假设beta2和beta3等于0后进行估计
# 对应P128倒数第二段和P129第一段
##### 测试：根据fit的估计，写出fit_H0的估计 #####
fit_H0 <- mle2(LL_H0, start = list(beta0=0, beta1 = 1, mu=0, sigma=1), fixed = list(mu = 0))
summary(fit_H0)

##########- LR test -##########
# 似然比检验：LR=有约束似然值/无约束似然值。LR越接近1，说明所得模型是满足约束条件的。
# 2.7.3.3
LR_stats <- anova(fit_H0, fit)
print(paste("  LR:", round(LR_stats[2, "Chisq"], 2), "  ", round(LR_stats[2, "Pr(>Chisq)"], 2)))

##########- Wald Test -##########
# Wald检验：如果约束是有效的，那么在没有约束情况下估计出来的估计量应该渐进地满足约束条件，因为MLE是一致的
# 2.7.3.2
R_mtrx <- matrix(c(0,0,0,0,1,0,0,1,0,0), 2, 5)
r_vct <- c(0, 0)
wald_stats <- wald_test(R_mtrx, coef(fit)[c(1:4, 6)], vcov(fit), r_vct)
print(paste("Wald:", round(wald_stats[1], 2), "  ", round(wald_stats[3], 2)))

##########- Lagrange Multiplier tests -##########
# 拉格朗日乘子检验（LM）：在约束条件下，可以用拉格朗日方法构造目标函数，
# 如果约束有效，则最大化拉格朗日函数所得估计量应位于最大化无约束所得参数估计值附近。
# 2.7.3.4

# 有约束的估计量
var_tilde <- coef(fit_H0)["sigma"] ^ 2  # var_tilde = sigma_tilde^2
epsilon_tilde <- stk_rtn - 
  MKT * (coef(fit_H0)["beta1"] / (coef(fit_H0)["sigma"] ^ 2)) - 
  (coef(fit_H0)["beta0"] / var_tilde)
s_tilde_col_1 <- 1 / var_tilde * MKT * epsilon_tilde
s_tilde_col_2 <- - 1 / (2 * var_tilde) + (epsilon_tilde ^ 2 / (2 * var_tilde ^ 2))
s_tilde <- cbind(s_tilde_col_1, s_tilde_col_2)
s_tilde_bar <- colMeans(s_tilde)
S_cap_tilde <- (t(s_tilde) %*% s_tilde) / length(stk_rtn)

# 没有约束的估计量
var_hat <- coef(fit)["sigma"] ^ 2  # var_gat = sigma_hat^2
epsilon_hat <- stk_rtn - 
  MKT * (coef(fit)["beta1"] / (coef(fit)["sigma"] ^ 2)) - 
  SMB * (coef(fit)["beta2"] / (coef(fit)["sigma"] ^ 2)) - 
  HML * (coef(fit)["beta3"] / (coef(fit)["sigma"] ^ 2)) - 
  (coef(fit)["beta0"] / var_tilde)
s_hat_col_1 <- cbind(1/var_hat * MKT * epsilon_hat)
s_hat_col_2 <- - 1 / (2 * var_hat) + (epsilon_hat ^ 2 / (2 * var_hat ^ 2))
s_hat <- cbind(s_hat_col_1, s_hat_col_2)
s_hat_bar <- colMeans(s_tilde)
S_cap_hat <- (t(s_hat) %*% s_hat) / length(stk_rtn)

# 统计量，P129第二个公式：LM=···
LM_stats_tilde <- length(stk_rtn) * s_tilde_bar %*% solve(S_cap_tilde) %*% s_tilde_bar
LM_pval_tilde <- 1 - pchisq(LM_stats_tilde, 2)
print(paste("  LM(S_tilde):", round(LM_stats_tilde, 2), "  ", round(LM_pval_tilde, 2)))

##### 测试：写出LM(S_hat)的估计式  #####
LM_stats_hat <- length(stk_rtn) * s_hat_bar %*% solve(S_cap_hat) %*% s_hat_bar
LM_pval_hat <- 1 - pchisq(LM_stats_hat, 2)
print(paste("  LM(S_hat):", round(LM_stats_hat, 2), "  ", round(LM_pval_hat, 2)))


##########- 结果输出 -##########
cat(c(names(coef(fit))[1], ": ", round(coef(fit)[1], 4)), "\n", 
    c(names(coef(fit))[2], ": ", round(coef(fit)[2], 4)), "\n", 
    c(names(coef(fit))[3], ": ", round(coef(fit)[3], 4)), "\n", 
    c(names(coef(fit))[4], ": ", round(coef(fit)[4], 4)), "\n", 
    c(names(coef(fit))[5], ": ", round(coef(fit)[5], 4)), "\n", 
    paste("LR:", round(LR_stats[2, "Chisq"], 2), ",", round(LR_stats[2, "Pr(>Chisq)"], 2)), "\n",
    paste("Wald:", round(wald_stats[1], 2), ",", round(wald_stats[3], 2)), "\n",
    paste("LM(S_tilde):", round(LM_stats_tilde, 2), ",", round(LM_pval_tilde, 2)), "\n",
    paste("LM(S_hat):", round(LM_stats_hat, 2), ",", round(LM_pval_hat, 2)), "\n", 
    sep="")




