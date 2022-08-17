library("car")

compute.tau.hat_R_unadj<-function(Y, Z){
  trtd = Y[Z == 1]; ctrl = Y[Z == 0]
  that = median(outer(trtd, ctrl, FUN = "-"))
  return(that) 
}

library("rio")
Data = import("PROGRESA.csv")
Data$pri2000s = Data$pri2000s

Y = Data$pri2000s
Z = Data$treatment

R.unadj = median(outer(Y[Z == 1], Y[Z == 0], FUN = "-"))

m = sum(Z); N = length(Z)
aux<-function(tau){
  fit = lm(I(Y - tau*Z) ~ as.factor(villages) + pri1994 + pan1994 + prd1994 +
                          votos1994 + avgpoverty	+ pobtot1994, data = Data)
  return(sum(rank(fit$res)*Z)/N - (m/N)*(N+1)/2)
}

R.adj = uniroot(aux, c(0, 1), extendInt="yes")$root

bh = Y - R.unadj*Z
bdiff = outer(bh, bh, "-")
Ib.hat = N^(-3/2) * sum((bdiff >= 0) * (bdiff < 1/sqrt(N)))
se_R.unadj = sqrt((12*(m/N)*(1 - m/N)*Ib.hat^2)^(-1)/N)

fit = lm(I(Y - R.adj*Z) ~ as.factor(villages) + pri1994 + pan1994 + prd1994 +
                      votos1994 + avgpoverty	+ pobtot1994, data = Data)
bh = fit$res
bdiff = outer(bh, bh, "-")
Jb.hat = N^(-3/2) * sum((bdiff >= 0) * (bdiff < 1/sqrt(N)))
se_R.adj = sqrt((12*(m/N)*(1 - m/N)*Jb.hat^2)^(-1)/N)


dm = mean(Y[Z == 1]) - mean(Y[Z == 0])
se_dm = sqrt(var(Y[Z == 1])/m + var(Y[Z == 0])/(N-m)) 


l_fit1 = lm(pri2000s ~ treatment + as.factor(villages) + pri1994 + pan1994 + prd1994 +
             votos1994 + avgpoverty	+ pobtot1994, data=Data) 
l.adj = l_fit1$coef[2]
se_l.adj = sqrt(diag(hccm(l_fit1, type = "hc0"))[2])

vars = c('pri1994', 'pan1994', 'prd1994', 'votos1994', 'avgpoverty', 'pobtot1994')
for(v in vars){ 
  Data[,paste0(v,1)] = Data[,v] - mean(Data[,v])
}

l_fit2 = lm(pri2000s ~ treatment*(pri19941 + pan19941 + prd19941 +
                                    votos19941 + avgpoverty1 + pobtot19941)+as.factor(villages), data=Data) 
l.inter = l_fit2$coef[2]
se_l.inter = sqrt(diag(hccm(l_fit2, type = "hc0"))[2])

est = cbind(R.unadj, dm, R.adj, l.adj, l.inter)
stderr = cbind(se_R.unadj, se_dm, se_R.adj, se_l.adj, se_l.inter)
q = qnorm(1-0.025)
M.s = t(rbind(est, stderr, est - q*stderr, est + q*stderr, 2*q*stderr))
colnames(M.s) = c("estimate", "std.error", "CI.lower", "CI.upper", "CI.length")
summary(fit)$r.sq
summary(l_fit1)$r.sq
summary(l_fit2)$r.sq

M.s # outcome variable is pri2000s (support)

library(xtable)
print(xtable(M.s, type = "latex", digits = 3))
