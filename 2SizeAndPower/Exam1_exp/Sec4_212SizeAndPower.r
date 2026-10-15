# rm(list = ls())
# gc()
library(MASS) # 多元正态分布
library(Rcpp) # R链接c++语言
library(ggplot2) # 绘图
library(reshape2) # 数据融合
library(compiler)
library(survival)
# 并行计算
library(foreach)
library(iterators)
library(parallel)
library(doParallel) # 仅此一个也行
sourceCpp('TestVCM.cpp')
# 遍历所有 .R 文件并导入
r_files <- list.files('speedyCMI/')
for (file in r_files) {
    source(paste0('speedyCMI/', file))    
}

version

# combine_vector_counts <- function(prev, curr) {
#   # prev：累计计数列表（含4个位置的统计）
#   # curr：当前迭代的长度为4的向量
#   for (pos in 1:9) {  # 遍历向量的4个位置
#     val <- curr[pos]  # 当前位置的值
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
GenData1 <- function(n = 100, lam = 0.2, C1 = 0, C2 = 0, std = 1){
  X <- matrix(0, n, 2)
  X[, 1] <- rnorm(n, 1, 1)
  X[, 2] <- runif(n, 0.5, 1.5)
  V <- runif(n, min = 0, max = 2)
  U <- runif(n, min = 0, max = 2)
  beta12_U <- 0.5 * U * (3 - U) + 1
  beta1_U <- cbind(4 - 3 * cos((U - 0.5) * pi / 2) / 4, beta12_U)
  beta2_U <- 2 * U^2
  Y <- apply(X * beta1_U, 1, sum) + V * beta2_U + C1 * X[, 1]^2 * X[, 2] * U / 2 + (0.5 * beta12_U * X[, 2])^C2 * rnorm(n, sd = std)
  C <- (0.7 - 3/2 * lam) * exp(4 * X[, 1]) # 0.1: 40%; 0.4: 20#
  delta <- 1 * (V <= C)
  T1 <- pmin(V, C)
  # cat(mean(delta))
  return(list(Y = Y, X = X, T1 = T1, delta = delta, U = U, V = V))
}


# method <- 'Test_VCM'
# n <- 100
# lam = 0.4
# C1 = 0
# C2 = 0
# std = 1
# 
# m = 40
# diff_alpha = 0.05
# B = 10
# h_con = 2.34
# n_con = 3
# Data <- GenData1(n = n, lam = lam, C1 = C1, C2 = C2, std = std)
# h_seq = seq(0.01, 1, by = 0.1)
# df <-list(Y = Data$Y, X = Data$X, T1 = Data$T1, delta = Data$delta, U = Data$U, h_seq = h_seq)
# method1 <- 'bw_select'
# c(do.call(method1, df))
# 
# 
# dat = data.frame(z = Data$X[, 1], w = Data$T1, y = Data$Y, d = Data$delta)
# imp_dat1 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
#                              dist = "exponential",
#                              data = dat)
# imp_dat2 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
#                              dist = "weibull",
#                              data = dat)
# df <- append(Data, list(VE = imp_dat1$imputed_data$imp, VW = imp_dat2$imputed_data$imp, B = B, alpha = diff_alpha, n_con = n_con, m = m))
# c(do.call(method, df))






# 
###########################################################
method <- 'Test_VCM'
n <- 200 # 样本量
lam <- 0.2
m <- 50
diff_method = c('Ideal', 'Proposed', 'CC', 'Naive', 'Impute_CMI_exp', 'Impute_CMI_wei','Impute_Mean', 'KS_N', 'CvM_N')
C1 <- 0 # 偏离
method_num <- length(diff_method)
C2 <- 0 # 同方差否
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std <- 1
diff_alpha <- 0.05 # 名义水平
diff_h_con <- c(1.54, 1.74, 1.94, 3.74, 3.94, 4.14)
n_con <- 3
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
t1 = Sys.time()
ncores <- detectCores() - 1
cl <- makeCluster(ncores) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <-
  foreach(h_con = diff_h_con, .combine = cbind) %:%
  foreach(i = 1:Sim, .combine = '+', #.combine = combine_vector_counts,
          #.init = lapply(1:method_num, function(x) list(count_0 = 0, count_1 = 0, count_e = 0)),
          .noexport = c(method, 'cmi_fp_analytical'),
          .packages = c('Rcpp', 'MASS', 'survival'), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            for (file in r_files) {
              source(paste0('speedyCMI/', file))
            }
            Data <- GenData1(n = n, lam = lam, C1 = C1, C2 = C2, std = std)
            dat = data.frame(z = Data$X[, 1], w = Data$T1, y = Data$Y, d = Data$delta)
            imp_dat1 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "exponential",
                                         data = dat)
            imp_dat2 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "weibull",
                                         data = dat)
            df <- append(Data, list(VE = imp_dat1$imputed_data$imp, VW = imp_dat2$imputed_data$imp, B = B, alpha = diff_alpha, h_con = h_con, n_con = n_con, m = m))
            c(do.call(method, df))
          }
stopCluster(cl) # 关闭
t2 <- Sys.time()
m1 <- length(diff_h_con)
m2 <- length(diff_method)

re <- t(round(reDGP / Sim, 5))
rownames(re) <- paste0('dn', ':', diff_h_con)
colnames(re) <- diff_method
# re
re <- round(re, 5)
rownames(re) <- paste0('c_n', ':', diff_h_con)
colnames(re) <- diff_method
result <- list('GenData' = GenData1, '样本量, 仿真次数, 自举次数, 删失率' = c(n, Sim, B, lam),
               'mode, std, level' = c(mode, std, diff_alpha), 'C2, m' = c(C2, '\t', m),
               '耗时' = t2 - t1, 'SizeAndPower' = re)
result
# 保存
SaveDR <- c("SizeAndPower/")
capture.output(result, file = paste0(SaveDR, "exam1_result", '_', mode, ".txt"), append = TRUE)
###############################################################


# 
# ###########################################################
# method <- 'Test_VCM'
# n <- 300 # 样本量
# lam <- 0.4
# m <- 50
# diff_method = c('Ideal', 'Proposed', 'CC', 'Naive', 'Impute_CMI_exp', 'Impute_CMI_wei','Impute_Mean', 'KS_N', 'CvM_N')
# # C1 <- 0 # 偏离
# diff_dn <- c(0, 0, 0.2, 0.4, 0.6, 0.8, 1.0)
# method_num <- length(diff_method)
# C2 <- 0 # 同方差否
# if (C2 == 0) {
#   mode = 'Homo'
# } else {
#   mode = 'Heter'
# }
# std <- 1
# diff_alpha <- 0.05 # 名义水平
# h_con <- 2.34
# n_con <- 3
# Sim <- 500 # 仿真次数
# B <- 300 # Bootstrap次数
# t1 = Sys.time()
# ncores <- detectCores() - 1
# cl <- makeCluster(ncores) # 核心数, 少调用一个
# registerDoParallel(cl)
# reDGP <-
#   foreach(dn = diff_dn, .combine = cbind) %:%
#   foreach(i = 1:Sim, .combine = '+', #.combine = combine_vector_counts,
#           #.init = lapply(1:method_num, function(x) list(count_0 = 0, count_1 = 0, count_e = 0)),
#           .noexport = c(method, 'cmi_fp_analytical'),
#           .packages = c('Rcpp', 'MASS', 'survival'), .errorhandling = 'remove') %dopar%{
#             sourceCpp('TestVCM.cpp')
#             for (file in r_files) {
#               source(paste0('speedyCMI/', file))
#             }
#             Data <- GenData1(n = n, lam = lam, C1 = dn, C2 = C2, std = std)
#             dat = data.frame(z = Data$X[, 1], w = Data$T1, y = Data$Y, d = Data$delta)
#             imp_dat1 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
#                                          dist = "exponential",
#                                          data = dat)
#             imp_dat2 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
#                                          dist = "weibull",
#                                          data = dat)
#             df <- append(Data, list(VE = imp_dat1$imputed_data$imp, VW = imp_dat2$imputed_data$imp, B = B, alpha = diff_alpha, h_con = h_con, n_con = n_con, m = m))
#             c(do.call(method, df))
#           }
# stopCluster(cl) # 关闭
# t2 <- Sys.time()
# m1 <- length(diff_dn)
# m2 <- length(diff_method)
# 
# re <- t(round(reDGP / Sim, 5))
# rownames(re) <- paste0('dn', ':', diff_dn)
# colnames(re) <- diff_method
# # re
# re <- round(re, 5)
# rownames(re) <- paste0('c_n', ':', diff_dn)
# colnames(re) <- diff_method
# result <- list('GenData' = GenData1, '样本量, 仿真次数, 自举次数, 删失率' = c(n, Sim, B, lam),
#                'mode, std, level' = c(mode, std, diff_alpha), 'C2, m' = c(C2, '\t', m),
#                '耗时' = t2 - t1, 'SizeAndPower' = re)
# result
# # 保存
# SaveDR <- c("SizeAndPower/")
# capture.output(result, file = paste0(SaveDR, "exam1_result", '_', mode, ".txt"), append = TRUE)
# ##############################################################




###########################################################
method <- 'Test_VCM'
n <- 200 # 样本量
lam <- 0.2
m <- 50
diff_method = c('Ideal', 'Proposed', 'CC', 'Naive', 'Impute_CMI_exp', 'Impute_CMI_wei','Impute_Mean', 'KS_N', 'CvM_N')
# C1 <- 0 # 偏离
diff_dn <- c(0, 0, 0)
method_num <- length(diff_method)
C2 <- 1 # 同方差否
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std <- 1
diff_alpha <- 0.05 # 名义水平
h_con <- 2.34
n_con <- 3
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
t1 = Sys.time()
ncores <- detectCores() - 1
cl <- makeCluster(ncores) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <-
  foreach(dn = diff_dn, .combine = cbind) %:%
  foreach(i = 1:Sim, .combine = '+', #.combine = combine_vector_counts,
          #.init = lapply(1:method_num, function(x) list(count_0 = 0, count_1 = 0, count_e = 0)),
          .noexport = c(method, 'cmi_fp_analytical'),
          .packages = c('Rcpp', 'MASS', 'survival'), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            for (file in r_files) {
              source(paste0('speedyCMI/', file))
            }
            Data <- GenData1(n = n, lam = lam, C1 = dn, C2 = C2, std = std)
            dat = data.frame(z = Data$X[, 1], w = Data$T1, y = Data$Y, d = Data$delta)
            imp_dat1 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "exponential",
                                         data = dat)
            imp_dat2 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "weibull",
                                         data = dat)
            df <- append(Data, list(VE = imp_dat1$imputed_data$imp, VW = imp_dat2$imputed_data$imp, B = B, alpha = diff_alpha, h_con = h_con, n_con = n_con, m = m))
            c(do.call(method, df))
          }
stopCluster(cl) # 关闭
t2 <- Sys.time()
m1 <- length(diff_dn)
m2 <- length(diff_method)

re <- t(round(reDGP / Sim, 5))
rownames(re) <- paste0('dn', ':', diff_dn)
colnames(re) <- diff_method
# re
re <- round(re, 5)
rownames(re) <- paste0('c_n', ':', diff_dn)
colnames(re) <- diff_method
result <- list('GenData' = GenData1, '样本量, 仿真次数, 自举次数, 删失率' = c(n, Sim, B, lam),
               'mode, std, level' = c(mode, std, diff_alpha), 'C2, m' = c(C2, '\t', m),
               '耗时' = t2 - t1, 'SizeAndPower' = re)
result
# 保存
SaveDR <- c("SizeAndPower/")
capture.output(result, file = paste0(SaveDR, "exam1_result", '_', mode, ".txt"), append = TRUE)
##############################################################


###########################################################
method <- 'Test_VCM'
n <- 200 # 样本量
lam <- 0.4
m <- 50
diff_method = c('Ideal', 'Proposed', 'CC', 'Naive', 'Impute_CMI_exp', 'Impute_CMI_wei','Impute_Mean', 'KS_N', 'CvM_N')
# C1 <- 0 # 偏离
diff_dn <- c(0, 0, 0)
method_num <- length(diff_method)
C2 <- 1 # 同方差否
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std <- 1
diff_alpha <- 0.05 # 名义水平
h_con <- 2.34
n_con <- 3
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
t1 = Sys.time()
ncores <- detectCores() - 1
cl <- makeCluster(ncores) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <-
  foreach(dn = diff_dn, .combine = cbind) %:%
  foreach(i = 1:Sim, .combine = '+', #.combine = combine_vector_counts,
          #.init = lapply(1:method_num, function(x) list(count_0 = 0, count_1 = 0, count_e = 0)),
          .noexport = c(method, 'cmi_fp_analytical'),
          .packages = c('Rcpp', 'MASS', 'survival'), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            for (file in r_files) {
              source(paste0('speedyCMI/', file))
            }
            Data <- GenData1(n = n, lam = lam, C1 = dn, C2 = C2, std = std)
            dat = data.frame(z = Data$X[, 1], w = Data$T1, y = Data$Y, d = Data$delta)
            imp_dat1 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "exponential",
                                         data = dat)
            imp_dat2 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "weibull",
                                         data = dat)
            df <- append(Data, list(VE = imp_dat1$imputed_data$imp, VW = imp_dat2$imputed_data$imp, B = B, alpha = diff_alpha, h_con = h_con, n_con = n_con, m = m))
            c(do.call(method, df))
          }
stopCluster(cl) # 关闭
t2 <- Sys.time()
m1 <- length(diff_dn)
m2 <- length(diff_method)

re <- t(round(reDGP / Sim, 5))
rownames(re) <- paste0('dn', ':', diff_dn)
colnames(re) <- diff_method
# re
re <- round(re, 5)
rownames(re) <- paste0('c_n', ':', diff_dn)
colnames(re) <- diff_method
result <- list('GenData' = GenData1, '样本量, 仿真次数, 自举次数, 删失率' = c(n, Sim, B, lam),
               'mode, std, level' = c(mode, std, diff_alpha), 'C2, m' = c(C2, '\t', m),
               '耗时' = t2 - t1, 'SizeAndPower' = re)
result
# 保存
SaveDR <- c("SizeAndPower/")
capture.output(result, file = paste0(SaveDR, "exam1_result", '_', mode, ".txt"), append = TRUE)
##############################################################



###########################################################
method <- 'Test_VCM'
n <- 200 # 样本量
lam <- 0.4
m <- 50
diff_method = c('Ideal', 'Proposed', 'CC', 'Naive', 'Impute_CMI_exp', 'Impute_CMI_wei','Impute_Mean', 'KS_N', 'CvM_N')
# C1 <- 0 # 偏离
diff_dn <- c(0, 0, 0)
method_num <- length(diff_method)
C2 <- 0 # 同方差否
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std <- 1
diff_alpha <- 0.05 # 名义水平
h_con <- 2.34
n_con <- 3
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
t1 = Sys.time()
ncores <- detectCores() - 1
cl <- makeCluster(ncores) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <-
  foreach(dn = diff_dn, .combine = cbind) %:%
  foreach(i = 1:Sim, .combine = '+', #.combine = combine_vector_counts,
          #.init = lapply(1:method_num, function(x) list(count_0 = 0, count_1 = 0, count_e = 0)),
          .noexport = c(method, 'cmi_fp_analytical'),
          .packages = c('Rcpp', 'MASS', 'survival'), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            for (file in r_files) {
              source(paste0('speedyCMI/', file))
            }
            Data <- GenData1(n = n, lam = lam, C1 = dn, C2 = C2, std = std)
            dat = data.frame(z = Data$X[, 1], w = Data$T1, y = Data$Y, d = Data$delta)
            imp_dat1 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "exponential",
                                         data = dat)
            imp_dat2 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "weibull",
                                         data = dat)
            df <- append(Data, list(VE = imp_dat1$imputed_data$imp, VW = imp_dat2$imputed_data$imp, B = B, alpha = diff_alpha, h_con = h_con, n_con = n_con, m = m))
            c(do.call(method, df))
          }
stopCluster(cl) # 关闭
t2 <- Sys.time()
m1 <- length(diff_dn)
m2 <- length(diff_method)

re <- t(round(reDGP / Sim, 5))
rownames(re) <- paste0('dn', ':', diff_dn)
colnames(re) <- diff_method
# re
re <- round(re, 5)
rownames(re) <- paste0('c_n', ':', diff_dn)
colnames(re) <- diff_method
result <- list('GenData' = GenData1, '样本量, 仿真次数, 自举次数, 删失率' = c(n, Sim, B, lam),
               'mode, std, level' = c(mode, std, diff_alpha), 'C2, m' = c(C2, '\t', m),
               '耗时' = t2 - t1, 'SizeAndPower' = re)
result
# 保存
SaveDR <- c("SizeAndPower/")
capture.output(result, file = paste0(SaveDR, "exam1_result", '_', mode, ".txt"), append = TRUE)
##############################################################



###########################################################
method <- 'Test_VCM'
n <- 200 # 样本量
lam <- 0.2
m <- 50
diff_method = c('Ideal', 'Proposed', 'CC', 'Naive', 'Impute_CMI_exp', 'Impute_CMI_wei','Impute_Mean', 'KS_N', 'CvM_N')
# C1 <- 0 # 偏离
diff_dn <- c(0, 0, 0)
method_num <- length(diff_method)
C2 <- 0 # 同方差否
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std <- 1
diff_alpha <- 0.05 # 名义水平
h_con <- 2.34
n_con <- 3
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
t1 = Sys.time()
ncores <- detectCores() - 1
cl <- makeCluster(ncores) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <-
  foreach(dn = diff_dn, .combine = cbind) %:%
  foreach(i = 1:Sim, .combine = '+', #.combine = combine_vector_counts,
          #.init = lapply(1:method_num, function(x) list(count_0 = 0, count_1 = 0, count_e = 0)),
          .noexport = c(method, 'cmi_fp_analytical'),
          .packages = c('Rcpp', 'MASS', 'survival'), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            for (file in r_files) {
              source(paste0('speedyCMI/', file))
            }
            Data <- GenData1(n = n, lam = lam, C1 = dn, C2 = C2, std = std)
            dat = data.frame(z = Data$X[, 1], w = Data$T1, y = Data$Y, d = Data$delta)
            imp_dat1 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "exponential",
                                         data = dat)
            imp_dat2 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "weibull",
                                         data = dat)
            df <- append(Data, list(VE = imp_dat1$imputed_data$imp, VW = imp_dat2$imputed_data$imp, B = B, alpha = diff_alpha, h_con = h_con, n_con = n_con, m = m))
            c(do.call(method, df))
          }
stopCluster(cl) # 关闭
t2 <- Sys.time()
m1 <- length(diff_dn)
m2 <- length(diff_method)

re <- t(round(reDGP / Sim, 5))
rownames(re) <- paste0('dn', ':', diff_dn)
colnames(re) <- diff_method
# re
re <- round(re, 5)
rownames(re) <- paste0('c_n', ':', diff_dn)
colnames(re) <- diff_method
result <- list('GenData' = GenData1, '样本量, 仿真次数, 自举次数, 删失率' = c(n, Sim, B, lam),
               'mode, std, level' = c(mode, std, diff_alpha), 'C2, m' = c(C2, '\t', m),
               '耗时' = t2 - t1, 'SizeAndPower' = re)
result
# 保存
SaveDR <- c("SizeAndPower/")
capture.output(result, file = paste0(SaveDR, "exam1_result", '_', mode, ".txt"), append = TRUE)
##############################################################





###########################################################
method <- 'Test_VCM'
n <- 200 # 样本量
lam <- 0.4
m <- 50
diff_method = c('Ideal', 'Proposed', 'CC', 'Naive', 'Impute_CMI_exp', 'Impute_CMI_wei','Impute_Mean', 'KS_N', 'CvM_N')
# C1 <- 0 # 偏离
diff_dn <- c(0, 0, 0)
method_num <- length(diff_method)
C2 <- 0 # 同方差否
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std <- 1
diff_alpha <- 0.05 # 名义水平
h_con <- 2.34
n_con <- 3
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
t1 = Sys.time()
ncores <- detectCores() - 1
cl <- makeCluster(ncores) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <-
  foreach(dn = diff_dn, .combine = cbind) %:%
  foreach(i = 1:Sim, .combine = '+', #.combine = combine_vector_counts,
          #.init = lapply(1:method_num, function(x) list(count_0 = 0, count_1 = 0, count_e = 0)),
          .noexport = c(method, 'cmi_fp_analytical'),
          .packages = c('Rcpp', 'MASS', 'survival'), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            for (file in r_files) {
              source(paste0('speedyCMI/', file))
            }
            Data <- GenData1(n = n, lam = lam, C1 = dn, C2 = C2, std = std)
            dat = data.frame(z = Data$X[, 1], w = Data$T1, y = Data$Y, d = Data$delta)
            imp_dat1 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "exponential",
                                         data = dat)
            imp_dat2 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "weibull",
                                         data = dat)
            df <- append(Data, list(VE = imp_dat1$imputed_data$imp, VW = imp_dat2$imputed_data$imp, B = B, alpha = diff_alpha, h_con = h_con, n_con = n_con, m = m))
            c(do.call(method, df))
          }
stopCluster(cl) # 关闭
t2 <- Sys.time()
m1 <- length(diff_dn)
m2 <- length(diff_method)

re <- t(round(reDGP / Sim, 5))
rownames(re) <- paste0('dn', ':', diff_dn)
colnames(re) <- diff_method
# re
re <- round(re, 5)
rownames(re) <- paste0('c_n', ':', diff_dn)
colnames(re) <- diff_method
result <- list('GenData' = GenData1, '样本量, 仿真次数, 自举次数, 删失率' = c(n, Sim, B, lam),
               'mode, std, level' = c(mode, std, diff_alpha), 'C2, m' = c(C2, '\t', m),
               '耗时' = t2 - t1, 'SizeAndPower' = re)
result
# 保存
SaveDR <- c("SizeAndPower/")
capture.output(result, file = paste0(SaveDR, "exam1_result", '_', mode, ".txt"), append = TRUE)
##############################################################



###########################################################
method <- 'Test_VCM'
n <- 300 # 样本量
lam <- 0.4
m <- 50
diff_method = c('Ideal', 'Proposed', 'CC', 'Naive', 'Impute_CMI_exp', 'Impute_CMI_wei','Impute_Mean', 'KS_N', 'CvM_N')
# C1 <- 0 # 偏离
diff_dn <- c(0, 0)
method_num <- length(diff_method)
C2 <- 1 # 同方差否
if (C2 == 0) {
  mode = 'Homo'
} else {
  mode = 'Heter'
}
std <- 1
diff_alpha <- 0.05 # 名义水平
h_con <- 2.34
n_con <- 3
Sim <- 500 # 仿真次数
B <- 300 # Bootstrap次数
t1 = Sys.time()
ncores <- detectCores() - 1
cl <- makeCluster(ncores) # 核心数, 少调用一个
registerDoParallel(cl)
reDGP <-
  foreach(dn = diff_dn, .combine = cbind) %:%
  foreach(i = 1:Sim, .combine = '+', #.combine = combine_vector_counts,
          #.init = lapply(1:method_num, function(x) list(count_0 = 0, count_1 = 0, count_e = 0)),
          .noexport = c(method, 'cmi_fp_analytical'),
          .packages = c('Rcpp', 'MASS', 'survival'), .errorhandling = 'remove') %dopar%{
            sourceCpp('TestVCM.cpp')
            for (file in r_files) {
              source(paste0('speedyCMI/', file))
            }
            Data <- GenData1(n = n, lam = lam, C1 = dn, C2 = C2, std = std)
            dat = data.frame(z = Data$X[, 1], w = Data$T1, y = Data$Y, d = Data$delta)
            imp_dat1 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "exponential",
                                         data = dat)
            imp_dat2 = cmi_fp_analytical(imputation_model = Surv(time = w, event = d) ~ z,
                                         dist = "weibull",
                                         data = dat)
            df <- append(Data, list(VE = imp_dat1$imputed_data$imp, VW = imp_dat2$imputed_data$imp, B = B, alpha = diff_alpha, h_con = h_con, n_con = n_con, m = m))
            c(do.call(method, df))
          }
stopCluster(cl) # 关闭
t2 <- Sys.time()
m1 <- length(diff_dn)
m2 <- length(diff_method)

re <- t(round(reDGP / Sim, 5))
rownames(re) <- paste0('dn', ':', diff_dn)
colnames(re) <- diff_method
# re
re <- round(re, 5)
rownames(re) <- paste0('c_n', ':', diff_dn)
colnames(re) <- diff_method
result <- list('GenData' = GenData1, '样本量, 仿真次数, 自举次数, 删失率' = c(n, Sim, B, lam),
               'mode, std, level' = c(mode, std, diff_alpha), 'C2, m' = c(C2, '\t', m),
               '耗时' = t2 - t1, 'SizeAndPower' = re)
result
# 保存
SaveDR <- c("SizeAndPower/")
capture.output(result, file = paste0(SaveDR, "exam1_result", '_', mode, ".txt"), append = TRUE)
##############################################################



























