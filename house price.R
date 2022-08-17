library("car")

compute.tau.hat_R_unadj<-function(Y, Z){
  trtd = Y[Z == 1]; ctrl = Y[Z == 0]
  that = median(outer(trtd, ctrl, FUN = "-"))
  return(that) 
}

compute.tsadj<-function(Y, Z, X, tau)
  return(sum(rank(lm((Y - tau*Z) ~ X)$res)*Z)/length(Z)) #with intercept

compute.tau.hat_R_adj<-function(Y, Z, X){
  m = sum(Z); N = length(Z)
  fun<-function(t) return(compute.tsadj(Y, Z, X, t) - (m/N)*(N+1)/2)
  return(uniroot(fun, c(0, 1), extendInt="yes")$root)
}

estimates<-function(Y, Z, X){
  X = as.matrix(X)
  m = sum(Z); N = length(Z)
  Xw = cbind(1, X)
  PX = Xw %*% solve(t(Xw) %*% Xw) %*% t(Xw)
  
  tau_unadj = compute.tau.hat_R_unadj(Y, Z)
  
  bh = Y - tau_unadj*Z
  bdiff = outer(bh, bh, "-")
  Ib.hat = N^(-3/2) * sum((bdiff >= 0) * (bdiff < 1/sqrt(N)))
  se_R.unadj = sqrt((12*(m/N)*(1 - m/N)*Ib.hat^2)^(-1)/N)
  
  tau_adj = compute.tau.hat_R_adj(Y, Z, X)  
  bh = (diag(N) - PX) %*% (Y - tau_adj*Z)
  bdiff = outer(bh, bh, "-")
  Jb.hat = N^(-3/2) * sum((bdiff >= 0) * (bdiff < 1/sqrt(N)))
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

# House price
library("rio")
Data = import("/Users/adityaghosh/Desktop/House price data/20060620_stata_datsets/analysis_sample.dta")
head(Data)
tail(Data)
names(Data)
dim(Data)

Mydat = Data[c("amt_Price", "log_price", "sale_year", "HEATED", "BEDROOMS", "AIRCOND", "BQM1", "BQM2", "AGE")]
Mydat = na.omit(Mydat)

library(moments)
kurtosis(log(Data$amt_Price))

set.seed(2027)
L = 1000; K = 1000
Dat = Mydat[sample(1:(dim(Mydat)[1]), K, replace = F), ]
system.time({
u = v = c()
for(l in 1:L){
Trt.index = sample(1:K, K/2, replace = F)
Z = rep(0, K)
Z[Trt.index] = 1
Y = log(Dat$amt_Price)

X = as.matrix(Dat[c("sale_year", "HEATED", "BEDROOMS", "AIRCOND", "BQM1", "BQM2", "AGE")])
#fit = lm(Y ~ Z +X)
#summary(fit)

q = qnorm(0.025, lower.tail = F)
out = estimates(Y, Z, X)
Mat = cbind(out$est - q*out$se, out$est + q*out$se, 2*q*out$se)
rownames(Mat) = c("R.un", "dm", "R.adj", "L.adj", "L.int")
colnames(Mat) = c("ci.low", "ci.upp", "ci.length")

u = rbind(u, apply(Mat, 1, function(x) as.numeric((x[1]<0)&&(x[2]>0))))
v = rbind(v, Mat[,3])
}
})

# coverage probs:
x = apply(u, 2, sum)/L
# avg C.I. lengths:
y = print(apply(v, 2, mean), digits = 2)

library(xtable)
print(xtable(rbind(x, y), type = "latex"), digits = 3)



set.seed(1234)
K = 1000
Dat = Mydat[sample(1:(dim(Mydat)[1]), K, replace = F), ]
Trt.index = sample(1:K, K/2, replace = F)
Z = rep(0, K)
Z[Trt.index] = 1
Y = log(Dat$amt_Price)
X = as.matrix(Dat[c("sale_year", "HEATED", "BEDROOMS", "AIRCOND", "BQM1", "BQM2", "AGE")])
fit = lm(Y ~ Z +X)
summary(fit)
q = qnorm(0.025, lower.tail = F)
out = estimates(Y, Z, X)  
M = cbind(out$est, out$se, out$est - q*out$se, out$est + q*out$se, 2*q*out$se)
colnames(M) = c("Estimate", "Std.error", "CI.lower", "CI.upper", "CI.length")
rownames(M) = c("Rosenbaum.unadj", "diff-in-means", "Rosenbaum.adj", "Lin.w/o.int", "Lin.interact")
M 

library(xtable)
print(xtable(M, type = "latex"))
