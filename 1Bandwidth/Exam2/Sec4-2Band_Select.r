rm(list = ls())
gc()
library(MASS) # 多元正态分布
library(Rcpp) # R链接c++语言
library(ggplot2) # 绘图
library(reshape2) # 数据融合
library(compiler)
# 并行计算
library(foreach)
library(iterators)
library(parallel)
library(doParallel) # 仅此一个也行
sourceCpp('TestVCM.cpp')
# 带Bootstrap的统计量仿真函数: testCPLRM()

version



# combine_vector_counts <- function(prev, curr) {
#   # prev：累计计数列表（含4个位置的统计）
#   # curr：当前迭代的长度为4的向量
#   for (pos in 1) {  # 遍历向量的4个位置
#     val <- curr[1, ][[pos]]  # 当前位置的值
#     # print(val)
#     # 根据值更新对应计数
#     if (is.na(val)) { # 报错复杂, 该项不会累加 #   if (val == error_flag) 
#       prev[[pos]]$count_e <- prev[[pos]]$count_e + 1
#     } else if (val == 1) {
#       prev[[pos]]$count_1 <- prev[[pos]]$count_1 + 1
#     } else if (val == 0) {
#       prev[[pos]]$count_0 <- prev[[pos]]$count_0 + 1
#     }
#   }
#   prev  # 返回更新后的计数列表
# }

# GenData2 <- function(n = 100, lam = 0.2, C1 = 0, C2 = 0, std = 1){
#   X <- matrix(0, n, 2)
#   X[, 1] <- rnorm(n, 1, 1)
#   X[, 2] <- runif(n, 0.5, 1.5)
#   V <- runif(n, min = 0, max = 2)
#   U <- runif(n, min = 0, max = 2)
#   beta12_U <- 0.5 * U * (3 - U) + 1
#   beta1_U <- cbind(4 - 3 * cos((U - 0.5) * pi / 2) / 4, beta12_U)
#   beta2_U <- 2 * U^2
#   Y <- apply(X * beta1_U, 1, sum) + V * beta2_U + C1 * X[, 1]^2 * X[, 2] * U / 2 + (0.5 * beta12_U * X[, 2])^C2 * rnorm(n, sd = std)
#   C <- (0.7 - 3/2 * lam) * exp(4 * X[, 1]) # 0.1: 40%; 0.4: 20#
#   delta <- 1 * (V <= C)
#   T1 <- pmin(V, C)
#   # cat(mean(delta))
#   return(list(Y = Y, X = X, T1 = T1, delta = delta, U = U))
# }

GenData2 <- function(n = 100, lam = 0.2, C1 = 0, C2 = 0, std = 1){
    p <- 7
    Sigma <- matrix(0.5, p, p)
    diag(Sigma) <- 1
    X <- mvrnorm(n, rep(1, p), Sigma)
    V <- runif(n, min = 0, max = 6)
    U <- rnorm(n) #runif(n, min = 0, max = 1)#
    beta1 <- 1 * U^2 # # beta1 <- 4 * U - 1
    beta2 <- U^3 / 3
    beta3 <- 1 * U * (U - 1)
    beta4 <- 4.5 * cos(U * pi / 2)
    beta5 <- 4.5 * sin(U * pi / 2)
    beta6 <- exp(U / 2)
    beta7 <- 4 * log(abs(U))
    beta8 <- 1 / (1 + exp(-U))
    betaU <- cbind(beta1, beta2, beta3, beta4, beta5, beta6, beta7, beta8) * 2
    XV <- cbind(X, V)
#     cat(sd(apply(XV * betaU, 1, sum)), sd(X[, 2] * beta2), sd(V * beta8), sd(0.4 * beta8 * X[, 1]^2))
    Y <- apply(XV * betaU, 1, sum) + C1 * apply(X^2, 1, mean) * V + (beta8 * X[, 1]^2)^C2 * rnorm(n, sd = std)
    C <- rgamma(n, shape = 1, rate = lam / 2) / (1 + 2 * lam) * exp(0.8)
    delta <- 1 * (V <= C)
    T1 <- pmin(V, C)
#     cat(mean(delta))
    return(list(Y = Y, X = X, T1 = T1, delta = delta, U = U))
}

print2 <- function(reDGP, diff_con, diff_dn) {
    d1 <- method_num * length(diff_con)
    d2 <- length(diff_dn)
    ind1 <- seq(1, (2 * d1 * d2), 2)
    ind2 <- seq(2, (2 * d1 * d2), 2)
    ind3 <- c(gl(d2, 1, d1 * d2))
    reG1 <- reDGP[ind1,]
    reG2 <- reDGP[ind2,]
    Sim_real <- nrow(reG2)
    for (k in 1:(d1 * d2)) {
        outindex <- which(reG1[k, ] %in% boxplot.stats(reG1[k, ], coef = 2)$out) # out = 1.5 / diff_dn[ind3[k]]
        if (length(outindex) != 0) reG2[k, outindex] <- 2
    }

    re <- round(matrix(apply(reG2, 1, FUN = function(x) {sum(x == 1) / (Sim - sum(x == 2))}), d1, d2, byrow = TRUE), 4)
    rownames(re) <- paste0('h_n', ':', diff_con)
    colnames(re) <- paste0('dn', ':', diff_dn)
    rownames(re) <- paste0('h_n', ':', rep(diff_con, each = method_num))
    colnames(re) <- paste0('d_n', ':', diff_dn)
    result <- list('GenData' = GenData2, '样本量, 删失率, 仿真次数, 自举次数' = c(n, lam, Sim, B),
                   'mode, std, level' = c(mode, std, diff_alpha), 'C, C2' = c(diff_con, '\t', C2, '\t'),
                   '耗时' = t2 - t1, #'reDGP' = reDGP,
                   'SizeAndPower' = re)
    # 保存
    SaveDR <- paste0("Band_Test_", mode, '/')#c("Band_Test_Heter/")
    # sink(paste0(SaveDR, "exam_result", ".txt"), append = TRUE, split = FALSE)
    # result
    # sink()
    capture.output(result, file = paste0(SaveDR, mode, '_', n, '_', "exam1_result", ".txt"), append = TRUE)
    # 生成要提取的行索引：从1开始，步长为4
    write.csv(re, file = paste0(SaveDR, n, "+0.", lam * 10, ".csv"))
    
    MinIndex <- which.min(abs(as.matrix(re[, 1]) - diff_alpha)) # 给出得是首次最小的!!!!!
    cat('最佳带宽(不一定):')
    cat(rownames(re)[MinIndex]) # 原假设下最优带宽, 根据图示结果判断更好.
    return(result)
}

# C1 = 0
# C2 = 0
# n = 100
# lam = 0.4
# diff_alpha = 0.5
# method = 'Boot_VcM' # 'Boot_NonPara_LocLin1' # 
# method_num = 1
# std = 1
# n_con = 3
# B = 10
# con = 2.34
# Data <- GenData2(n = n, lam = lam, C1 = C1, C2 = C2, std = std) 
# df <- append(Data, 
#              list(alpha = diff_alpha, h_con = con, n_con = n_con, B = B))
# tryCatch({
#   do.call(method, df)
# }, error = function(e) {
#   return(rep(NA, method_num)) })

#########################################################################################

#########################################################################################
method = 'Boot_VcM' # 'Boot_NonPara_LocLin1' # 
stat_test <- c('CVM')
method_num <- length(stat_test)
### 样本量:  n = 4 00; 删失率:  lambda= 0.2, 0.4
n <- 200 # 样本量
lam <- 0.4 # diff.lam <- c(0.2, 0.4) # 删失率
diff_con <- c(seq(0, 2, by = 0.2) + 2.34)#[1:4] # 带宽
diff_dn <- c(0, 0) # c(0.9, 1.2, 1.5) #c(0, 0.3, 0.6)#, 0.9, 1.2, 1.5) ## 偏离
diff_alpha <- c(0.05) # 名义水平
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
C2 <- 0 # 异方差是否1
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std = 1
n_con <- 3

t1 = Sys.time()
cl <- makeCluster(detectCores()) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <- 
  foreach(con = diff_con, .combine = rbind) %:%
  foreach(dn = diff_dn, .combine = rbind) %:%
  foreach(i = 1:Sim, .combine = cbind,
          .packages = c('Rcpp', 'MASS'), .noexport = c(method), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            Data <- GenData2(n = n, lam = lam, C1 = dn, C2 = C2, std = std) 
            df <- append(Data, 
                         list(alpha = diff_alpha, h_con = con, n_con = n_con, B = B))
            tryCatch({
              do.call(method, df)
            }, error = function(e) {
              return(rep(NA, method_num)) })
          }
stopCluster(cl) # 关闭

t2 <- Sys.time()

print2(reDGP, diff_con, diff_dn)

#########################################################################################
method = 'Boot_VcM' # 'Boot_NonPara_LocLin1' # 
stat_test <- c('CVM')
method_num <- length(stat_test)
### 样本量:  n = 4 00; 删失率:  lambda= 0.2, 0.4
n <- 200 # 样本量
lam <- 0.4 # diff.lam <- c(0.2, 0.4) # 删失率
diff_con <- c(seq(0, 2, by = 0.2) + 2.34)#[1:4] # 带宽
diff_dn <- c(0, 0) # c(0.9, 1.2, 1.5) #c(0, 0.3, 0.6)#, 0.9, 1.2, 1.5) ## 偏离
diff_alpha <- c(0.05) # 名义水平
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
C2 <- 1 # 异方差是否1
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std = 1
n_con <- 3

t1 = Sys.time()
cl <- makeCluster(detectCores()) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <- 
  foreach(con = diff_con, .combine = rbind) %:%
  foreach(dn = diff_dn, .combine = rbind) %:%
  foreach(i = 1:Sim, .combine = cbind,
          .packages = c('Rcpp', 'MASS'), .noexport = c(method), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            Data <- GenData2(n = n, lam = lam, C1 = dn, C2 = C2, std = std) 
            df <- append(Data, 
                         list(alpha = diff_alpha, h_con = con, n_con = n_con, B = B))
            tryCatch({
              do.call(method, df)
            }, error = function(e) {
              return(rep(NA, method_num)) })
          }
stopCluster(cl) # 关闭

t2 <- Sys.time()

print2(reDGP, diff_con, diff_dn)

#########################################################################################
method = 'Boot_VcM' # 'Boot_NonPara_LocLin1' # 
stat_test <- c('CVM')
method_num <- length(stat_test)
### 样本量:  n = 4 00; 删失率:  lambda= 0.2, 0.4
n <- 100 # 样本量
lam <- 0.4 # diff.lam <- c(0.2, 0.4) # 删失率
diff_con <- c(seq(0, 2, by = 0.2) + 2.34)#[1:4] # 带宽
diff_dn <- c(0, 0, 0) # c(0.9, 1.2, 1.5) #c(0, 0.3, 0.6)#, 0.9, 1.2, 1.5) ## 偏离
diff_alpha <- c(0.05) # 名义水平
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
C2 <- 0 # 异方差是否1
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std = 1
n_con <- 3

t1 = Sys.time()
cl <- makeCluster(detectCores()) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <- 
  foreach(con = diff_con, .combine = rbind) %:%
  foreach(dn = diff_dn, .combine = rbind) %:%
  foreach(i = 1:Sim, .combine = cbind,
          .packages = c('Rcpp', 'MASS'), .noexport = c(method), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            Data <- GenData2(n = n, lam = lam, C1 = dn, C2 = C2, std = std) 
            df <- append(Data, 
                         list(alpha = diff_alpha, h_con = con, n_con = n_con, B = B))
            tryCatch({
              do.call(method, df)
            }, error = function(e) {
              return(rep(NA, method_num)) })
          }
stopCluster(cl) # 关闭

t2 <- Sys.time()

print2(reDGP, diff_con, diff_dn)

#########################################################################################
method = 'Boot_VcM' # 'Boot_NonPara_LocLin1' # 
stat_test <- c('CVM')
method_num <- length(stat_test)
### 样本量:  n = 4 00; 删失率:  lambda= 0.2, 0.4
n <- 100 # 样本量
lam <- 0.2 # diff.lam <- c(0.2, 0.4) # 删失率
diff_con <- c(seq(0, 2, by = 0.2) + 2.34)#[1:4] # 带宽
diff_dn <- c(0, 0, 0) # c(0.9, 1.2, 1.5) #c(0, 0.3, 0.6)#, 0.9, 1.2, 1.5) ## 偏离
diff_alpha <- c(0.05) # 名义水平
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
C2 <- 0 # 异方差是否1
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std = 1
n_con <- 3

t1 = Sys.time()
cl <- makeCluster(detectCores()) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <- 
  foreach(con = diff_con, .combine = rbind) %:%
  foreach(dn = diff_dn, .combine = rbind) %:%
  foreach(i = 1:Sim, .combine = cbind,
          .packages = c('Rcpp', 'MASS'), .noexport = c(method), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            Data <- GenData2(n = n, lam = lam, C1 = dn, C2 = C2, std = std) 
            df <- append(Data, 
                         list(alpha = diff_alpha, h_con = con, n_con = n_con, B = B))
            tryCatch({
              do.call(method, df)
            }, error = function(e) {
              return(rep(NA, method_num)) })
          }
stopCluster(cl) # 关闭

t2 <- Sys.time()

print2(reDGP, diff_con, diff_dn)



#########################################################################################
method = 'Boot_VcM' # 'Boot_NonPara_LocLin1' # 
stat_test <- c('CVM')
method_num <- length(stat_test)
### 样本量:  n = 4 00; 删失率:  lambda= 0.2, 0.4
n <- 100 # 样本量
lam <- 0.4 # diff.lam <- c(0.2, 0.4) # 删失率
diff_con <- c(seq(0, 2, by = 0.2) + 2.34)#[1:4] # 带宽
diff_dn <- c(0, 0, 0) # c(0.9, 1.2, 1.5) #c(0, 0.3, 0.6)#, 0.9, 1.2, 1.5) ## 偏离
diff_alpha <- c(0.05) # 名义水平
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
C2 <- 1 # 异方差是否1
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std = 1
n_con <- 3

t1 = Sys.time()
cl <- makeCluster(detectCores()) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <- 
  foreach(con = diff_con, .combine = rbind) %:%
  foreach(dn = diff_dn, .combine = rbind) %:%
  foreach(i = 1:Sim, .combine = cbind,
          .packages = c('Rcpp', 'MASS'), .noexport = c(method), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            Data <- GenData2(n = n, lam = lam, C1 = dn, C2 = C2, std = std) 
            df <- append(Data, 
                         list(alpha = diff_alpha, h_con = con, n_con = n_con, B = B))
            tryCatch({
              do.call(method, df)
            }, error = function(e) {
              return(rep(NA, method_num)) })
          }
stopCluster(cl) # 关闭

t2 <- Sys.time()

print2(reDGP, diff_con, diff_dn)







