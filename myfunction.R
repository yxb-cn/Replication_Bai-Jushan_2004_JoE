##### Functions to replicate Bai(2004)JoE #####
####### Code by Xuanbin Yang; 2022/11/12; at Nankai University ########


# Given k, estimate F (factors) in standard (restricted) dynamic factor model
# by principal component method
estimate_F = function(X, k, normalization = 1){ 
  #---------------------------------------------------
  #
  # Input:
  #   X: Sample matrix
  #   k: Number of factors
  #   normalization: By default, normalization = 1. 
  #      if normalization = 1, F'F/T^2=I (suitable for T<N), 
  #      if normalization = 2, L'L/N=I (suitable for N<T), 
  #      if normalization = 3, chose normalization method automatically 
  # Output: 
  #   A list containing: 
  #     F_hat: Factors
  #     Lambda_hat: Factor loading
  #     X_hat: Estimated X
  #     V_NT: A diagonal matrix consisting of eigenvalues, for calculating 
  #           confidence intervals of factors in standard (restricted) dynamic   
  #           factor model 
  #---------------------------------------------------
  T = dim(X)[1]
  N = dim(X)[2]
  if(normalization == 1){
    svdresult = eigs_sym(tcrossprod(X), k, which = "LM")
    F_hat = svdresult$vectors[,1:k]*T
    Lambda_hat = t(X)%*%F_hat/T^2
  }else if(normalization == 2){
    svdresult = eigs_sym(crossprod(X), k, which = "LM")
    Lambda_hat = svdresult$vectors[,1:k]*sqrt(N)
    F_hat = X%*%Lambda_hat/N
  }else{
    if(T<N){
      svdresult = eigs_sym(tcrossprod(X), k, which = "LM")
      F_hat = svdresult$vectors[,1:k]*T
      Lambda_hat = t(X)%*%F_hat/T^2
    }else{
      svdresult = eigs_sym(crossprod(X), k, which = "LM")
      Lambda_hat = svdresult$vectors[,1:k]*sqrt(N)
      F_hat = X%*%Lambda_hat/N
    }
  }
  V_NT = diag(svdresult$values/(T^2*N))
  return(list(F_hat = F_hat, Lambda_hat = Lambda_hat, X_hat = F_hat%*%t(Lambda_hat),V_NT = V_NT))
}


# Given r and q, estimate F (factors) in generalized dynamic factor model
# by principal component method
estimate_F_GDFM = function(X, r, q){
  #---------------------------------------------------
  #
  # Input:
  #   X: Sample matrix
  #   r: Number of I(1) factors
  #   q: Number of I(0) factors
  # Output:
  #   A list containing:
  #     F_hat: Factors
  #     Lambda_hat: Factor loading
  #     X_hat: Estimated X
  #     V_NT: A diagonal matrix consisting of eigenvalues, for calculating
  #           confidence intervals of factors in generalized dynamic
  #           factor model
  #---------------------------------------------------
  options(warn = -1)
  T = dim(X)[1]
  N = dim(X)[2]
  k = r+q
  svdresult = eigs_sym(tcrossprod(X), k, which = "LM")
  D_1T = diag(c(rep(T,r),rep(sqrt(T),q)))
  F_hat = svdresult$vectors%*%D_1T
  Upsilon = diag(c(rep(T,r),rep(sqrt(T),q))^(-2))
  Lambda_hat = t(X)%*%F_hat%*%Upsilon
  V_NT = diag(c(svdresult$values[1:r]/(T^2*N),svdresult$values[(r+1):(r+q)]/(T*N)))
  return(list(F_hat = F_hat, Lambda_hat = Lambda_hat, X_hat = F_hat%*%t(Lambda_hat),V_NT = V_NT))
}

# The sum of squared residuals (divided by NT) when k trends are estimated
V = function(X,k,normalization = 1){
  #---------------------------------------------------
  #
  # Input:
  #   X: Sample matrix
  #   k: Number of factors
  #   normalization: By default, normalization = 1.
  #      if normalization = 1, F'F/T^2=I (suitable for T<N), 
  #      if normalization = 2, L'L/N=I (suitable for N<T), 
  #      if normalization = 3, chose normalization method automatically 
  # Output: 
  #   V_k: The sum of squared residuals (divided by NT)
  #---------------------------------------------------
  T = dim(X)[1]
  N = dim(X)[2]
  estimationlist = estimate_F(X, k, normalization = normalization)
  X_hat = estimationlist$X_hat
  V_k = sum((X-X_hat)^2)/(N*T)
  return(V_k)
}

# criteria for determining the number of trends (factors)
PC_IPC = function(X,k,kmax,normalization = 1,datatype){
  #---------------------------------------------------
  #
  # Input:
  #   k: Number of factors
  #   X: Sample matrix
  #   kmax: Upper limit of number of factors
  #   normalization: By default, normalization = 1.
  #      if normalization = 1, F'F/T^2=I (suitable for T<N), 
  #      if normalization = 2, L'L/N=I (suitable for N<T), 
  #      if normalization = 3, chose normalization method automatically 
  #   datatype: datatype = 1 means using data in differences,
  #             datatype = 2 means using data in levels
  # Output: 
  #   PC_IPC_value: A vector whose elements are the three PC or IPC values
  #---------------------------------------------------
  T = dim(X)[1]
  N = dim(X)[2]
  PC_IPC_value = rep(0,3)
  if(datatype == 1){
    alpha = 1
  }else{
    alpha = T/(4*log(log(T)))
  }
  V_k = V(X,k,normalization = normalization)
  sigma_square = V(X,kmax,normalization = normalization)
  PC_IPC_value[1] = V_k+k*sigma_square*alpha*(N+T)/(N*T)*log(N*T/(N+T))
  PC_IPC_value[2] = V_k+k*sigma_square*alpha*(N+T)/(N*T)*log(min(N,T))
  PC_IPC_value[3] = V_k+k*sigma_square*alpha*(N+T-k)/(N*T)*log(N*T)
  return(PC_IPC_value)
}

# Estimate the number of factors
estimate_factornum_byPC_IPC = function(X,kmax,normalization = 1,datatype){
  #---------------------------------------------------
  #
  # Input:
  #   X: Sample matrix
  #   k: Number of factors
  #   kmax: Upper limit of number of factors
  #   normalization: By default, normalization = 1.
  #      if normalization = 1, F'F/T^2=I (suitable for T<N),
  #      if normalization = 2, L'L/N=I (suitable for N<T),
  #      if normalization = 3, chose normalization method automatically
  #   datatype: datatype = 1 means using data in differences,
  #             datatype = 2 means using data in levels
  # Output:
  #   k_hat: A vector whose elements are the k_hat determined by the three PC or IPC
  #          criteria, respectively
  #---------------------------------------------------
  k_hat = rep(1,3)
  best = PC_IPC(X,1,kmax,normalization = normalization,datatype = datatype)

  for(k in 2:kmax){
    temp = PC_IPC(X,k,kmax,normalization = normalization,datatype = datatype)
    if(temp[1] < best[1]){
      k_hat[1] = k
      best[1] = temp[1]
    }
    if(temp[2] < best[2]){
      k_hat[2] = k
      best[2] = temp[2]
    }
    if(temp[3] < best[3]){
      k_hat[3] = k
      best[3] = temp[3]
    }
  }
  return(k_hat)
}


# Generate data from standard (restricted) dynamic factor model
generate_standard_DFM_data = function(N,T,r,rou,theta,seed = FALSE){
  #---------------------------------------------------
  #
  # Input:
  #   N: The second dimension of X
  #   T: The first dimension of X
  #   r: Parameter
  #   rou: Parameter
  #   theta: Parameter
  #   seed: Random seed, not used by default
  # Output:
  #   X_list: A list containing simulated X(T*N), differenced X(T-1*N) and the simulated F (named F0)
  #---------------------------------------------------
  Lambda = mvrnorm(N,rep(0,r),diag(r))
  if(seed){
    set.seed(seed)
  }
  burn_num = 500
  # Generate u,v, drop the first elements at last
  u = mvrnorm(T+burn_num,rep(0,r),diag(r))
  v = mvrnorm(T+burn_num,rep(0,N),diag(N))
  # Generate e, drop the first elements at last
  e = matrix(0,T+burn_num,N)
  for(t in 2:(T+burn_num)){
    e[t,] = rou*e[t-1,]+v[t,]+theta*v[t-1,]
  }
  # Generate F, drop the first elements at last
  F = matrix(0,T+burn_num,r)
  F[1,] = u[1,]
  for(t in 2:(T+burn_num)){
    F[t,] = F[t-1,]+u[t,]
  }
  # Generate X, and first difference of X
  F = F[(1+burn_num):(T+burn_num),]
  e = e[(1+burn_num):(T+burn_num),]
  X = F%*%t(Lambda)+e
  X_diff = diff(X)
  return(list(X = X,X_diff = X_diff,F0 = F))
}

# Generate data from generalized dynamic factor model 
generate_generalized_DFM_data = function(N,T,r,p,rou,theta,seed = FALSE){
  #---------------------------------------------------
  #
  # Input: 
  #   N: The second dimension of X
  #   T: The first dimension of X
  #   r: Parameter
  #   p: Parameter
  #   rou: Parameter
  #   theta: Parameter
  #   seed: Random seed, not used by default
  # Output: 
  #   X_list: A list containing simulated X(T*N), differenced X(T-1*N) and the simulated F (named F0)
  #---------------------------------------------------
  Lambda = array(0,c(N,r,p+1))
  for(k in 1:(p+1)){
    Lambda[,,k] = mvrnorm(N,rep(0,r),diag(r))
  }
  if(seed){
    set.seed(seed)
  }
  burn_num = 500
  # Generate u,v, drop the first elements at last
  u = mvrnorm(T+burn_num+p,rep(0,r),diag(r))
  v = mvrnorm(T+burn_num+1,rep(0,N),diag(N))
  # Generate e, drop the first elements at last
  e = matrix(0,T+burn_num,N)
  for(t in 2:(T+burn_num)){
    e[t,] = rou*e[t-1,]+v[t,]+theta*v[t-1,]
  }
  # Generate F, drop the first elements at last
  F = matrix(0,T+burn_num+p,r)
  F[1,] = u[1,]
  for(t in 2:(T+burn_num+p)){
    F[t,] = F[t-1,]+u[t,]
  }
  F_lag = array(0,c(T+burn_num,r,p+1))
  # F_lag[,,1] = F[(p+1):(T+burn_num+1),]
  for(k in 1:(p+1)){
    F_lag[,,k] = F[(p+2-k):(T+burn_num+p+1-k),]
  }
  # Generate X, and first difference of X
  F_lag = F_lag[(1+burn_num):(T+burn_num),,]
  e = e[(1+burn_num):(T+burn_num),]
  X = e
  for(k in 1:(p+1)){
    X = X + F_lag[,,k]%*%t(Lambda[,,k])
  }
  X_diff = diff(X)
  F = F[(2+burn_num):(T+burn_num+1),]
  return(list(X = X,X_diff = X_diff,F0 = F))
}

# Calculate the 95% confidence interval of standard dynamic factor model and the
# sample correlation coefficient between the rotated estimated factor and true factor
F_confidence_interval_SDFM <- function(X, F0, k, normalization = 1){
  #---------------------------------------------------
  #
  # Input: 
  #   X: Sample matrix
  #   F0: True factors,True factors, a matrix
  #   k: Number of factors to be estimated
  #   normalization: By default, normalization = 1.
  #      if normalization = 1, F'F/T^2=I (suitable for T<N), 
  #      if normalization = 2, L'L/N=I (suitable for N<T), 
  #      if normalization = 3, chose normalization method automatically 
  # Output: 
  #   A list containing confidence interval and sample correlation coefficient
  #---------------------------------------------------
  T = dim(X)[1]
  N = dim(X)[2]
  F0_num = dim(F0)[2]
  estimationlist = estimate_F(X, k, normalization = normalization)
  F_tilde = estimationlist$F_hat
  Lambda_tilde = estimationlist$Lambda_hat
  X_tilde = estimationlist$X_hat
  V_NT = estimationlist$V_NT
  e_tilde = X-X_tilde
  delta_hat = matrix(0,k,k)
  CI = array(0,c(T,k,2))
  S = matrix(0,T,k) 
  cor_F = rep(0,k)
  for(k0 in 1:F0_num){
    delta_hat[,k0] = solve(t(F_tilde)%*%F_tilde)%*%t(F_tilde)%*%F0[,k0]
    cor_F[k0] = cor(F_tilde%*%delta_hat[,k0],F0[,k0])
    for(t in 1:T){
      mid = 0
      for(i in 1:N){
        mid = mid + e_tilde[t,i]^2*as.matrix(Lambda_tilde[i,])%*%t(Lambda_tilde[i,])
      }
      mid = mid/N
      S[t,k0]  = sqrt(t(delta_hat[,k0])%*%solve(V_NT)%*%mid%*%solve(V_NT)%*%delta_hat[,k0])
    }
    CI[,k0,1] = F_tilde%*%delta_hat[,k0]-1.96*N^(-1/2)*matrix(S[,k0])
    CI[,k0,2] = F_tilde%*%delta_hat[,k0]+1.96*N^(-1/2)*matrix(S[,k0])
  }
  return(list(CI = CI, cor_F = cor_F))
}

# Calculate the 95% confidence interval of generalized dynamic factor model and the
# sample correlation coefficient between the rotated estimated factor and true factor
F_confidence_interval_GDFM <- function(X, F0, r, q){
  #---------------------------------------------------
  #
  # Input: 
  #   X: Sample matrix
  #   F0: True factors, a matrix
  #   r: Number of I(1) factors to be estimated
  #   q: Number of I(0) factors to be estimated
  # Output: 
  #   A list containing confidence interval and sample correlation coefficient
  #---------------------------------------------------
  T = dim(X)[1]
  N = dim(X)[2]
  F0_num = dim(F0)[2]
  estimationlist = estimate_F_GDFM(X, r, q)
  F_tilde = estimationlist$F_hat
  Lambda_tilde = estimationlist$Lambda_hat
  X_tilde = estimationlist$X_hat
  V_NT = estimationlist$V_NT
  e_tilde = X-X_tilde
  delta_hat = matrix(0,r+q,F0_num)
  CI = array(0,c(T,F0_num,2))
  S = matrix(0,T,r+q) 
  cor_F = rep(0,F0_num)
  for(k0 in 1:F0_num){
    delta_hat[,k0] = solve(t(F_tilde)%*%F_tilde)%*%t(F_tilde)%*%F0[,k0]
    cor_F[k0] = cor(F_tilde%*%delta_hat[,k0],F0[,k0])
    for(t in 1:T){
      mid = 0
      for(i in 1:N){
        mid = mid + e_tilde[t,i]^2*as.matrix(Lambda_tilde[i,])%*%t(Lambda_tilde[i,])
      }
      mid = mid/N
      S[t,k0]  = sqrt(t(delta_hat[,k0])%*%solve(V_NT)%*%mid%*%solve(V_NT)%*%delta_hat[,k0])
    }
    CI[,k0,1] = F_tilde%*%delta_hat[,k0]-1.96*N^(-1/2)*matrix(S[,k0])
    CI[,k0,2] = F_tilde%*%delta_hat[,k0]+1.96*N^(-1/2)*matrix(S[,k0])
  }
  return(list(CI = CI, cor_F = cor_F))
}

# Works just like Matlab's squeeze function: if anything in dim(x) equals one the corresponding dimension is removed
squeeze <- function(x) {
  d <- dim(x)
  dim(x) <- d[d>1]
  x
}
