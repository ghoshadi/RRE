#===========================================
# R codes for implementing our method
#===========================================

library("car")

# computing our rank-based estimator w/o regression adjustment:
compute.tau.hat_R_unadj <- function(v, Z){
  trtd = v[Z == 1]; ctrl = v[Z == 0]
  that = median(outer(trtd, ctrl, FUN = "-"))
  return(that) 
}

# computing regression Wilcoxon rank-sum statistic:
compute.WRS.adj <- function(Y, Z, X, tau)
  return(sum(rank(lm((Y - tau*Z) ~ X)$res)*Z)/length(Z)) #with intercept

# computing our rank-based regression adjusted estimator:
compute.tau.hat_R_adj <- function(Y, Z, X){
  m = sum(Z); N = length(Z)
  fun <- function(t) return(compute.WRS.adj(Y, Z, X, t) - (m/N)*(N+1)/2)
  return(uniroot(fun, c(0, 1), extendInt="yes")$root)
}

# auxiliary function to compute different estimates:
estimates <- function(Y, Z, X){
  m = sum(Z); n = length(Z)
  Xw = cbind(1, X)
  PX = Xw %*% solve(t(Xw) %*% Xw) %*% t(Xw)
  
  tau_unadj = compute.tau.hat_R_unadj(Y, Z)
  bh = Y - tau_unadj*Z
  bdiff = outer(bh, bh, "-")
  Ib.hat = N^(-3/2) * sum((bdiff > 0) * (bdiff <= 1/sqrt(N)))
  se_R.unadj = sqrt((12*(m/N)*(1 - m/N)*Ib.hat^2)^(-1)/N)
  
  tau_adj = compute.tau.hat_R_adj(Y, Z, X)  
  bh = (diag(N) - PX) %*% (Y - tau_adj*Z)
  bdiff = outer(bh, bh, "-")
  Jb.hat = N^(-3/2) * sum((bdiff > 0) * (bdiff <= 1/sqrt(N)))
  se_R.adj = sqrt((12*(m/N)*(1 - m/N)*Jb.hat^2)^(-1)/N)
  
  tau_dm = mean(Y[Z == 1]) - mean(Y[Z == 0])
  se_dm = sqrt(var(Y[Z == 1])/m + var(Y[Z == 0])/(N-m)) 
  
  X = scale(X)
  l_fit1 = lm(Y ~ Z + X)
  tau_ladj = l_fit1$coef[2]
  se_ladj = sqrt(diag(hccm(l_fit1, type = "hc0"))[2])
  
  l_fit = lm(Y ~ Z  + X + Z*X)
  tau_linter = l_fit$coef[2]
  se_linter = sqrt(diag(hccm(l_fit, type = "hc0"))[2])
  
  return(list(est = c(tau_unadj, tau_dm, tau_adj, tau_ladj, tau_linter), se = c(se_R.unadj, se_dm, se_R.adj, se_ladj, se_linter)))
}

# A simple example
set.seed(123)
N = 1e3; m = 0.4*N
tau0 = 2
X = runif(N, min = -4, max = 4)
a = (X + X^3)/4 + rnorm(N)
b = a - tau0

Z = as.numeric(sample(N) <= m)
Y = a*Z + b*(1-Z)

out = estimates(Y, Z, X)
Mat = cbind(out$est - tau0, out$est - 1.96*out$se, out$est + 1.96*out$se, 2*1.96*out$se)
rownames(Mat) = c("tau.hat.R", "diff-in-means", "tau.hat.R.adj", "Lin.reg.adj", "Lin.interact")
colnames(Mat) = c("bias", "ci.low", "ci.upp", "ci.length")
print(round(Mat, digits = 3))

#===========================================
# R codes for our simulations section
#===========================================

estimates.new <- function(Y, Z, X, PX, m, N, tau0 = 2){
  
  tau_unadj = compute.tau.hat_R_unadj(Y, Z)
  bh = Y - tau_unadj*Z
  bdiff = outer(bh, bh, "-")
  Ib.hat = N^(-3/2) * sum((bdiff > 0) * (bdiff <= 1/sqrt(N)))
  se_R.unadj = sqrt((12*(m/N)*(1 - m/N)*Ib.hat^2)^(-1)/N)
  
  tau_adj = compute.tau.hat_R_adj(Y, Z, X)  
  bh = (diag(N) - PX) %*% (Y - tau_adj*Z)
  bdiff = outer(bh, bh, "-")
  Jb.hat = N^(-3/2) * sum((bdiff > 0) * (bdiff <= 1/sqrt(N)))
  se_R.adj = sqrt((12*(m/N)*(1 - m/N)*Jb.hat^2)^(-1)/N)
  
  true.b.diff = outer(Y - tau0*Z, Y - tau0*Z, "-")
  I_N = N^(-3/2) * sum((true.b.diff > 0) * (true.b.diff <= 1/sqrt(N)))
  true.b.adj = (diag(N) - PX) %*% (Y - tau0*Z)
  true.b.adiff = outer(true.b.adj, true.b.adj, "-")
  J_N = N^(-3/2) * sum((true.b.adiff > 0) * (true.b.adiff <= 1/sqrt(N)))
  se_oracle.unadj = sqrt((12*(m/N)*(1 - m/N)*I_N^2)^(-1)/N)
  se_oracle.adj = sqrt((12*(m/N)*(1 - m/N)*J_N^2)^(-1)/N)
  
  tau_dm = mean(Y[Z == 1]) - mean(Y[Z == 0])
  se_dm = sqrt(var(Y[Z == 1])/m + var(Y[Z == 0])/(N-m)) 
  
  X = scale(X)
  l_fit1 = lm(Y ~ Z + X)
  tau_ladj = l_fit1$coef[2]
  se_ladj = sqrt(diag(hccm(l_fit1, type = "hc0"))[2])
  
  l_fit = lm(Y ~ Z  + X + Z*X)
  tau_linter = l_fit$coef[2]
  se_linter = sqrt(diag(hccm(l_fit, type = "hc0"))[2])
  
  return(c(tau_unadj, tau_dm, tau_adj, tau_ladj, tau_linter, se_R.unadj, se_oracle.unadj, se_dm, se_R.adj,  se_oracle.adj, se_ladj, se_linter))
}

aux <- function(a, tau0, x, Px, mlist){
  N = length(a)
  b = a - tau0
  new = c()
  for(m in mlist){ 
    Z = as.numeric(sample(N) <= m)
    Y = a*Z + b*(1-Z)
    new = c(new, estimates.new(Y, Z, x, Px, m, N))
  }
  return(new)
}

main <- function(){
  N = 1e3
  tau0 = 2
  mlist = c(0.75, 0.6, 0.5, 0.4, 0.25)*N
  
  Output = c()
  
  ## Setting 1: a = iid + eps
  
  x = runif(N, min = -4, max = 4)
  xw = cbind(1, x)
  Px = xw %*% solve(t(xw) %*% xw) %*% t(xw)
  
  a = rexp(N, rate = 1/10) + rnorm(N)
  Output = c(Output, aux(a, tau0, x, Px, mlist))
  
  a = rexp(N, rate = 1/10) + rt(N, df = 1)
  Output = c(Output, aux(a, tau0, x, Px, mlist))
  
  a = rexp(N, rate = 1/10) + rt(N, df = 3)
  Output = c(Output, aux(a, tau0, x, Px, mlist))
  
  ## Setting 2: a = 3x + eps
  
  a = 3*x + rnorm(N)
  Output = c(Output, aux(a, tau0, x, Px, mlist))
  
  a = 3*x + rt(N, df = 1)
  Output = c(Output, aux(a, tau0, x, Px, mlist))
  
  a = 3*x + rt(N, df = 3)
  Output = c(Output, aux(a, tau0, x, Px, mlist))
  
  ## Setting 3: a = (x + sqrt(x))/4 + eps
  
  x = exp(runif(N, min = -4, max = 4))
  xw = cbind(1, x)
  Px = xw %*% solve(t(xw) %*% xw) %*% t(xw)
  
  a = (x + sqrt(x))/4 + rnorm(N)
  Output = c(Output, aux(a, tau0, x, Px, mlist))
  
  a = (x + sqrt(x))/4 + rt(N, df = 1)
  Output = c(Output, aux(a, tau0, x, Px, mlist))
  
  a = (x + sqrt(x))/4 + rt(N, df = 3)
  Output = c(Output, aux(a, tau0, x, Px, mlist))
  
  return(Output) 
}

length(out <- main()) # returns a vector of length 540

# The above code was run repeatedly using parallel computing.
# The outputs were then stored in an excel file and summarized in Tables 3-5 of the main paper.
