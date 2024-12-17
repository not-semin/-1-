# 응용통계학 특수문제 과제 -----------------------------------------------------
# 201902086 안세민 -------------------------------------------------------------

# 0. Packages and seed  --------------------------------------------------------
library(survival)
library(dplyr)
library(MASS)

set.seed(2024311134)

# 1 - (a) ----------------------------------------------------------------------

#parameter
n <- 30000
rho <- 0.5
gamma <- 1.5
beta1 <- log(2)
beta2 <- log(2)
beta3 <- -1
mu <- c(0, 0)
Sigma <- matrix(c(1, 0.75, 0.75, 1), nrow=2)

simulate_cox_model <- function(iter) {
  
  beta_estimates <- matrix(NA, nrow = iter, ncol = 3)
  se_estimates <- matrix(NA, nrow = iter, ncol = 3)
  
  for (i in 1:iter) {
    # Z1, Z2, W1, W2, u 생성
    Z <- mvrnorm(n, mu = mu, Sigma = Sigma)
    Z1 <- Z[, 1]
    Z2 <- Z[, 2]
    W1 <- rnorm(n, 0, 1)
    W2 <- rbinom(n, 1, 0.5)
    u <- runif(n)
    
    # Ti 계산
    lambda <- exp(beta1 * Z1 + beta2 * W1 + beta3 * W2) * rho
    Ti <- (-log(1 - u) / lambda)^(1 / gamma)
    
    data <- data.frame(Ti = Ti, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
    
    cox_model <- coxph(Surv(Ti) ~ Z1 + W1 + W2, data = data)
    
    beta_estimates[i, ] <- cox_model$coefficients
    se_estimates[i, ] <- sqrt(diag(cox_model$var))
  }
  return(list(beta_estimates = beta_estimates, se_estimates = se_estimates))
}

result_1_a <- as.data.frame(simulate_cox_model(iter = 1))

result_1_a <- data.frame(
  rbind(
    c(round(result_1_a[1], 4), round(result_1_a[2], 4), round(result_1_a[3], 4)),
    c(round(result_1_a[4], 4), round(result_1_a[5], 4), round(result_1_a[6], 4))
  )
)
colnames(result_1_a) <- c("beta1", "beta2", "beta3")
rownames(result_1_a) <- c("estimates", "se")

# 1 - (b) ----------------------------------------------------------------------

result_1_b <- simulate_cox_model(iter = 100)

beta_means <- colMeans(result_1_b$beta_estimates)
beta_sds <- apply(result_1_b$beta_estimates, 2, sd)
se_means <- colMeans(result_1_b$se_estimates)

result_1_b <- data.frame(
  True_Beta = c(round(beta1,4), round(beta2,4), round(beta3,4)),
  Mean_Estimated_Beta = round(beta_means,4),
  Mean_SE = round(se_means,4),
  SD_Estimated_Beta = round(beta_sds,4)
)
rownames(result_1_b) <- c("Beta1", "Beta2", "Beta3")


# 1 - (c) ----------------------------------------------------------------------

# censoring rate
censoring_rates <- c(0.10, 0.30, 0.90, 0.95, 0.99)
lambda_C_values <- numeric(length(censoring_rates))
names(lambda_C_values) <- paste0(censoring_rates * 100, "%")

N <- 1e6

# censoring rate에 대한 lambda_C 찾기
find_lambda_C <- function(desired_censoring_rate) {
  
  objective_function <- function(lambda_C) {
    Z <- mvrnorm(N, mu = mu, Sigma = Sigma)
    Z1 <- Z[,1]
    W1 <- rnorm(N, 0, 1)
    W2 <- rbinom(N, 1, 0.5)
    u <- runif(N)
    lambda <- exp(beta1 * Z1 + beta2 * W1 + beta3 * W2) * rho
    Ti <- (-log(1 - u) / lambda)^(1 / gamma)
    Ci <- rexp(N, rate = lambda_C)
    censoring_rate <- mean(Ti > Ci)
    return((censoring_rate - desired_censoring_rate)^2)
  }
  # 최적화를 사용하여 lambda_C 찾기
  result <- optimize(objective_function, interval = c(0.0001, 10))
  return(result$minimum)
}

# 원하는 검열 비율에 대해 반복
for (i in seq_along(censoring_rates)) {
  lambda_C_values[i] <- find_lambda_C(censoring_rates[i])
}

# lambda_C 값 출력
result_1_c <- round(lambda_C_values,4)

# 1 - (d) ----------------------------------------------------------------------

# 검열을 포함한 시뮬레이션 함수
simulation1 <- function(lambda_C) {
  beta_estimates <- matrix(NA, nrow=100, ncol=3)
  se_estimates <- matrix(NA, nrow=100, ncol=3)

  for (i in 1:100) {
    
    Z <- mvrnorm(n, mu = mu, Sigma = Sigma)
    Z1 <- Z[,1]
    Z2 <- Z[,2]
    W1 <- rnorm(n, 0, 1)
    W2 <- rbinom(n, 1, 0.5)
    u <- runif(n)
    lambda <- exp(beta1 * Z1 + beta2 * W1 + beta3 * W2) * rho
    Ti <- (-log(1 - u) / lambda)^(1 / gamma)
    
    Ci <- rexp(n, rate = lambda_C)
    Xi <- pmin(Ti, Ci)
    Delta <- as.numeric(Ti <= Ci)
    
    data <- data.frame(Xi = Xi, Delta = Delta, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
    
    cox_model <- coxph(Surv(Xi, Delta) ~ Z1 + W1 + W2, data = data)
    
    beta_estimates[i, ] <- cox_model$coefficients
    se_estimates[i, ] <- sqrt(diag(cox_model$var))
  }
  
  beta_means <- colMeans(beta_estimates)
  beta_sds <- apply(beta_estimates, 2, sd)
  se_means <- colMeans(se_estimates)
  
  results <- data.frame(
    True_Beta = c(round(beta1,4), round(beta2,4), round(beta3,4)),
    Mean_Estimated_Beta = round(beta_means,4),
    Bias = c(round(beta1 - beta_means[1],4),
             round(beta2 - beta_means[2],4),
             round(beta3 - beta_means[3],4)),
    Mean_SE = round(se_means,4),
    SD_Estimated_Beta = round(beta_sds,4)
  )
  rownames(results) <- c("Beta1", "Beta2", "Beta3")
  
  return(list(results = results))  
}

# 각 censoring rate 시뮬레이션 수행
censoring_results <- list()

for (i in seq_along(censoring_rates)) {
  lambda_C <- lambda_C_values[i]
  sim_result <- simulation1(lambda_C)
  
  censoring_results[[i]] <- sim_result$results
}

censoring <- c("10%", "30%", "90%", "95%", "99%")

result_1_d <- data.frame()

for (i in 1:length(censoring_results)) {
  df <- censoring_results[[i]]
  rownames(df) <- paste0(rownames(df), "(", censoring[i], ")")
  
  result_1_d <- rbind(result_1_d, df)
}

# 2 - (a) - (i) ----------------------------------------------------------------

lambda_C <- lambda_C_values[5]

simulation2 <- function() {
  simulated_data_99 <- vector("list", 500) 
  
  # 500번 반복
  for (i in 1:500) {
    
    Z <- mvrnorm(n, mu = mu, Sigma = Sigma)
    Z1 <- Z[,1]
    Z2 <- Z[,2]
    W1 <- rnorm(n, 0, 1)
    W2 <- rbinom(n, 1, 0.5)
    u <- runif(n)
    lambda <- exp(beta1 * Z1 + beta2 * W1 + beta3 * W2) * rho
    Ti <- (-log(1 - u) / lambda)^(1 / gamma)
    
    Ci <- rexp(n, rate = lambda_C) 
    Xi <- pmin(Ti, Ci)
    Delta <- as.numeric(Ti <= Ci)
    
    data <- data.frame(Xi = Xi, Delta = Delta, Z1 = Z1, Z2 = Z2, W1 = W1, W2 = W2)
    
    # 500개의 데이터를 저장
    simulated_data_99[[i]] <- data
  }
  
  # 99% 검열 데이터 반환
  return(simulated_data_99)  
}

# 검열 비율이 99%인 경우에 대한 시뮬레이션 수행
simulated_data_99 <- simulation2()

num_pop <- 30000

#함수
perform_cox_analysis <- function(num_subcohort, iter) {
  # 결과를 저장할 리스트 및 행렬 초기화
  sampling_data_all <- list()
  sampling_results <- list()
  n_c <- numeric(iter)
  n_c_tilde <- numeric(iter)
  sampling_count_failure <- numeric(iter)
  sampling_count_sub <- numeric(iter)
  sampling_count_sub_failure <- numeric(iter)
  
  beta_estimates <- matrix(NA, nrow = iter, ncol = 3)
  se_estimates <- matrix(NA, nrow = iter, ncol = 3)
  
  for (i in 1:iter) {
    # 데이터를 샘플링
    sampling_data <- as.data.frame(simulated_data_99[i])
    temp_sampling_data <- sampling_data[sample(num_pop),]
    subco_ind <- c(rep(1, num_subcohort), rep(0, num_pop - num_subcohort))
    sampling_sub_data <- data.frame(temp_sampling_data, I(subco_ind))
    sampling_sub_data <- sampling_sub_data[order(sampling_sub_data[, 1]), 
                                           1:ncol(sampling_sub_data)]
    
    # 통계 계산
    sampling_count_failure[i] <- nrow(sampling_sub_data %>% filter(Delta == 1))
    sampling_count_sub[i] <- nrow(sampling_sub_data %>% filter(subco_ind == 1)) +
      nrow(sampling_sub_data %>% filter(Delta == 1 & subco_ind == 0))
    
    sampling_count_sub_failure[i] <- nrow(sampling_sub_data 
                                          %>% filter(Delta == 1 & subco_ind == 1))
    
    n_c[i] <- nrow(sampling_sub_data %>% filter(Delta == 0))
    n_c_tilde[i] <- nrow(sampling_sub_data %>% filter(Delta == 0 & subco_ind == 1))
    
    # 가중치 계산
    sampling_sub_data <- sampling_sub_data %>%
      mutate(weight = case_when(
        subco_ind == 0 & Delta == 0 ~ 0,
        subco_ind == 0 & Delta == 1 ~ 1,
        subco_ind == 1 & Delta == 0 ~ n_c[i] / n_c_tilde[i],
        subco_ind == 1 & Delta == 1 ~ 1
      ))
    
    sampling_sub_data$Z1_weighted <- ifelse(sampling_sub_data$weight > 0, 
                                            sampling_sub_data$Z1, NA)
    
    # 데이터 저장 및 모델 적합
    sampling_data_all[[i]] <- sampling_sub_data
    sampling_result <- coxph(Surv(Xi, Delta) ~ Z1_weighted + W1 + W2, 
                             data = sampling_sub_data, weights = weight)
    sampling_results[[i]] <- sampling_result
    
    # 계수 및 표준 오차 저장
    beta_estimates[i, ] <- sampling_result$coefficients
    se_estimates[i, ] <- sqrt(diag(sampling_result$var))
  }
  
  # 결과 반환
  return(list(
    sampling_data_all = sampling_data_all,
    sampling_results = sampling_results,
    beta_estimates = beta_estimates,
    se_estimates = se_estimates,
    count_failure = sampling_count_failure,
    count_sub = sampling_count_sub,
    count_sub_failure = sampling_count_sub_failure
  ))
}

result_2_a_1 <- perform_cox_analysis(100, 1)

# 2 - (a) - (ii) ---------------------------------------------------------------

result_2_a_2 <- data.frame(result_2_a_1$count_failure,
                           result_2_a_1$count_sub_failure,
                           result_2_a_1$count_sub)
colnames(result_2_a_2) <- c("Count_Failure", "Failure in subcohort", 
                            "Count_Subcohort")

# 2 - (a) - (iii) --------------------------------------------------------------

result_2_a_3_1 <- result_2_a_1$sampling_results

result_2_a_3_2 <- result_2_a_1$beta_estimates
result_2_a_3_2 <- t(data.frame(t(round(result_2_a_3_2,4))))
colnames(result_2_a_3_2) <- c("beta1", "beta2", "beta3")
rownames(result_2_a_3_2) <- "estimator"

# 2 - (a) - (iv) ---------------------------------------------------------------

result_2_a_4 <- perform_cox_analysis(100, 500)

# 2 - (a) - (v) ----------------------------------------------------------------

result_2_a_5 <- data.frame(mean(result_2_a_4$count_failure),
                           mean(result_2_a_4$count_sub_failure),
                           mean(result_2_a_4$count_sub))
colnames(result_2_a_5) <- c("Count_Failure", "Failure in subcohort", 
                            "Count_Subcohort")


# 2 - (a) - (vi) ---------------------------------------------------------------

beta_means <- colMeans(result_2_a_4$beta_estimates, na.rm = TRUE)
beta_sds <- apply(result_2_a_4$beta_estimates, 2, sd, na.rm = TRUE)
se_means <- colMeans(result_2_a_4$se_estimates, na.rm = TRUE)

result_2_a_6 <- data.frame(
  True_Beta = c(round(beta1,4), round(beta2,4), round(beta3,4)),
  Mean_Estimated_Beta = round(beta_means,4),
  Mean_SE = round(se_means,4),
  SD_Estimated_Beta = round(beta_sds,4)
)
rownames(result_2_a_6) <- c("Beta1", "Beta2", "Beta3")

# 2 - (a) - (vii) --------------------------------------------------------------

#  1번 추출
result_2_a_7_1_0 <- perform_cox_analysis(300, 1)

result_2_a_7_1_1 <- data.frame(result_2_a_7_1_0$count_failure,
                             result_2_a_7_1_0$count_sub_failure,
                             result_2_a_7_1_0$count_sub)
colnames(result_2_a_7_1_1) <- c("Count_Failure", "Failure in subcohort", 
                                "Count_Subcohort")

result_2_a_7_1_2 <- result_2_a_7_1_0$sampling_results

result_2_a_7_1_3 <- result_2_a_7_1_0$beta_estimates
result_2_a_7_1_3 <- t(data.frame(t(round(result_2_a_7_1_3,4))))
colnames(result_2_a_7_1_3) <- c("beta1", "beta2", "beta3")
rownames(result_2_a_7_1_3) <- "estimator"

# 500번 추출
result_2_a_7_2_0 <- perform_cox_analysis(100, 500)

result_2_a_7_2_1 <- data.frame(mean(result_2_a_7_2_0$count_failure),
                           mean(result_2_a_7_2_0$count_sub_failure),
                           mean(result_2_a_7_2_0$count_sub))
colnames(result_2_a_7_2_1) <- c("Count_Failure", "Failure in subcohort", 
                                "Count_Subcohort")

beta_means <- colMeans(result_2_a_7_2_0$beta_estimates, na.rm = TRUE)
beta_sds <- apply(result_2_a_7_2_0$beta_estimates, 2, sd, na.rm = TRUE)
se_means <- colMeans(result_2_a_7_2_0$se_estimates, na.rm = TRUE)

result_2_a_7_2_2 <- data.frame(
  True_Beta = c(round(beta1,4), round(beta2,4), round(beta3,4)),
  Mean_Estimated_Beta = round(beta_means,4),
  Mean_SE = round(se_means,4),
  SD_Estimated_Beta = round(beta_sds,4)
)
rownames(result_2_a_7_2_2) <- c("Beta1", "Beta2", "Beta3")

# All results ------------------------------------------------------------------

result_1_a
result_1_a
result_1_b
result_1_c
result_1_d

result_2_a_2
result_2_a_3_1
result_2_a_3_2
result_2_a_5
result_2_a_6
result_2_a_7_1_1
result_2_a_7_1_2
result_2_a_7_1_3
result_2_a_7_2_1
result_2_a_7_2_2
