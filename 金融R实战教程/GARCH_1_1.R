#############################
##########- 实验4 -##########
#############################

# GARCH(1, 1)的数值模拟

# 背景介绍：
# GARCH模型(Generalized AutoRegressive Conditional Heteroskedasticity)
# 自从Engle(1982)提出ARCH模型分析时间序列的异方差性以后，
# Bollerslev(1986)又提出了GARCH模型，其对误差的方差进行了进一步的建模，主要用于波动性的分析和预测
# 参考讲义，以及Hamilton(1994)第21章和Fan and Yao(2006)第三章
# GARCH的基本形式参考讲义的(34)(35)式，需要估计ω、α和β值

# 实验过程：
# 1. 先生成一个GARCH(1, 1)模型的参数
# 2. 由参数和随机项模拟出收益率序列
# 3. 估计模拟收益率序列的GARCH(1, 1)参数
# 4. 对比估计的参数和实际的参数，分析GARCH(1, 1)模型的估计精度

# 实验目标：
# 1. 理解GARCH模型的形式和思想，如何用R实现
# 2. 学习GARCH(1, 1)的估计方法，理解GARCH模型的估计精确度
# 3. 理解数值模拟的思想方法
# 4. 优化输出结果，便于理解

# 预期结果：
# 估计出来的GARCH参数以及和实际参数的对比，可以看到GARCH模型能够很好地对波动率建模和预测


#####- 1. 计算theta_hat -#####
library(fGarch)
setwd("~/../OneDrive/Works/金融科技兴趣小组/项目/Lab4_单变量GARCH/")
# 使用S&P500指数的收益率序列得到用于数值模拟的参数(omega_hat, alpha_hat, beta_hat)

# 读取数据，并计算S&P500指数的收益率
raw_data <- read.csv("./data.csv", stringsAsFactors=F)
SP500 <- raw_data[, "SP500"]
SP500_rtn <- SP500[2: length(SP500)] / SP500[1: (length(SP500) - 1)] - 1

# 对S&P500指数收益率进行GARCH(1, 1)估计，得到参数(omega_hat, alpha_hat, beta_hat)
# 使用“fGarch”包的garchFit函数进行估计
# 第一个参数为GARCH模型形式，使用GARCH(1, 1)记性估计
# 第二个参数为收益率序列
# trace=F表示估计过程不打印在屏幕上
raw_fit <- garchFit(~garch(1, 1), data=SP500_rtn, trace=F) 
# 将估计结果中的参数作为原始值
# 估计出来的raw_fit是一个S4对象，因此使用“@”来从raw_fit中取出需要的参数list
omega_hat <- raw_fit@fit$coef["omega"]
alpha_hat <- raw_fit@fit$coef["alpha1"]
beta_hat <- raw_fit@fit$coef["beta1"]


#####- 2. 数值模拟 -#####
# 使用得到的(omega_hat, alpha_hat, beta_hat)和随机项模拟出收益率序列

# 参数设定
total_num <- 5000  # 收益率的个数
set.seed(666)  # 设定随机数种子
simu_rtn <- rep(NA, total_num)  # 定义收益率序列
sigma <- rep(NA, total_num)  # 定义sigma序列
epsilon <- rnorm(total_num, 0 ,1)  # epsilon~N(0, 1)，生成服从标准正态分布的随机残差项

# t=1时的初始值
sigma[1] <- sqrt(omega_hat / (1 - alpha_hat - beta_hat))   # 使用sigma的长期均值作为其t=1时的值
simu_rtn[1] <- sigma[1] * epsilon[1]  # t=1时的收益率值

# 根据GARCH(1, 1)模型，使用上面的参数值和随机项循环得到模拟收益率
for (i in 2: total_num) {
  sigma[i] <- sqrt(omega_hat + alpha_hat * (simu_rtn[i - 1]) ^ 2 + beta_hat * (sigma[i - 1]) ^ 2)
  simu_rtn[i] <- sigma[i] * epsilon[i]
}


#####- 3. 计算模拟收益率下的参数  -#####
# 使用第二步得到的收益率序列，进行GARCH(1, 1)的估计

# 分别使用250, 500, 750, 1000个收益率进行估计，分析估计参数的收益率个数和参数精确度之间的关系
cal_num <- c(250, 500, 750, 1000)

# 定义输出结果的格式
result <- matrix(NA, length(cal_num) + 1, 4)
result[1, ] <- raw_fit@fit$coef  # 原始参数作为对比放入

# 使用循环的方式估计不同个数收益率序列的GARCH(1, 1)参数
for (i in 1: length(cal_num)) {
  cal_num_i <- cal_num[i]
  # 这里的估计方式和第一部分一致
  ##### 测试：根据前面garchFit的用法，写出下面的代码，要求每次使用cal_num中的数据数目 #####
  result[i + 1, ] <- garchFit(~garch(1, 1), data=simu_rtn[(5000 - cal_num_i + 1): 5000], trace=F)@fit$coef
}


#####- 4. 优化输出结果形式  -#####
# 我们关注omega, alpha1, beta1的估计差别
result <- result[, -1]
rownames(result) <- paste0("theta_", c("hat", as.character(cal_num)))
colnames(result) <- c("omega", "alpha1", "beta1")

# 计算数值模拟估计出来的参数值和原始参数值的差
result_diff <- matrix(NA, length(cal_num) + 1, 3)
colnames(result_diff) <- c("omega_diff", "alpha1_diff", "beta1_diff")
for(i in 1: length(cal_num)) {
  ##### 测试：得到估计值和实际值的差值  #####
  result_diff[i + 1, ] <- result[i + 1, ] - result[1, ]
}
result <- cbind(result, result_diff)  # 把结果合并在一起

# 输出打印
print(result)





