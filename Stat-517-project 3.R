set.seed(517)
# initial values from project 2
library(MASS)
library(mixtools)
library(Rfast)
rm(list=ls())
Pi=c(0.3,1-0.3)
mu1=c(2,3)
mu0=c(0,0)
S1=rbind(c(3,0),c(0,1))
S0=rbind(c(3,0),c(0,1))
mu1.hold=c(0,0)
mu0.hold=c(0,0)
cov1.hold=rbind(c(0,0),c(0,0))
cov0.hold=rbind(c(0,0),c(0,0))
SE.sim=c()
theta.sim=c()
j=1
R=100

# low prev pi = 0.3
while (j<101){
  Di=rbinom(100,1,0.30)
  Xi_A=c(rep(0,100))
  Xi_B=c(rep(0,100))
  data = data.frame(Xi_A,Xi_B,Di)
  for (i in 1:nrow(data)){
    if (data$Di[i] == 0){
      x = mvrnorm(1,mu0,S0)
      data$Xi_A[i] = x[1]
      data$Xi_B[i] = x[2]
    }
    else{
      x = mvrnorm(1,mu1,S1)
      data$Xi_A[i] = x[1]
      data$Xi_B[i] = x[2]
    }
  }
  X_i = cbind(data$Xi_A,data$Xi_B)
  
  PI=c(0.3+0.02,0.7-0.02)
  MU = list(mu1+c(0.01,0.01),mu0+c(0.01,0.01))
  SIGMA = list(S1+rbind(c(0.01,0.01),c(0.01,0.01)),S0+rbind(c(0.01,0.01),c(0.01,0.01)))
  
  out.test1.sim = mvnormalmixEM(x = X_i, lambda = PI, mu = MU, sigma = SIGMA, k = 2, arbmean = T, arbvar = T, epsilon = 1e-08, maxit = 1000)
  
  pi_1.sim = out.test1.sim$lambda[1]
  mu_1.sim = out.test1.sim$mu[[1]]
  mu_0.sim = out.test1.sim$mu[[2]]
  cov_1.sim = out.test1.sim$sigma[[1]]
  cov_0.sim = out.test1.sim$sigma[[2]]
  
  mu1.hold = mu1.hold + mu_1.sim
  mu0.hold = mu0.hold + mu_0.sim
  cov1.hold = cov1.hold + cov_1.sim
  cov0.hold = cov0.hold + cov_0.sim
  
  theta.sim = rbind(theta.sim,cbind(mu_1.sim[1], mu_1.sim[2],mu_0.sim[1], mu_0.sim[2],cov_1.sim[1,1], cov_1.sim[2,2],cov_1.sim[1,2], cov_0.sim[1,1], cov_0.sim[2,2], cov_0.sim[1,2],pi_1.sim))
  
  # Jackknife
  Pi.sim = vector()
  Mu_1.sim = matrix(rep(0,200),nrow = 100,ncol = 2)
  Mu_0.sim = matrix(rep(0,200),nrow = 100,ncol = 2)
  var.A1.sim = vector()
  var.B1.sim= vector()
  var.AB1.sim= vector()
  var.A0.sim= vector()
  var.B0.sim= vector()
  var.AB0.sim= vector()
  
  
  for (n in 1:100){
    X_i.n.sim=X_i[-n,]
    out.test1.n.sim = mvnormalmixEM(x = X_i.n.sim, lambda = PI, mu = MU, sigma = SIGMA, k = 2, arbmean = T, arbvar = T, epsilon = 1e-04, maxit = 1000) 
    Pi.sim[n]= out.test1.n.sim$lambda[1]
    Mu_1.sim[n,] = as.vector(out.test1.n.sim$mu[[1]])
    Mu_0.sim[n,] = as.vector(out.test1.n.sim$mu[[2]])
    var.A1.sim[n] = out.test1.n.sim$sigma[[1]][1,1]
    var.B1.sim[n] = out.test1.n.sim$sigma[[1]][2,2]
    var.AB1.sim[n] = out.test1.n.sim$sigma[[1]][1,2]
    var.A0.sim[n] = out.test1.n.sim$sigma[[2]][1,1]
    var.B0.sim[n] = out.test1.n.sim$sigma[[2]][2,2]
    var.AB0.sim[n] = out.test1.n.sim$sigma[[2]][1,2]
    MLEs.sim = cbind(Mu_1.sim,Mu_0.sim,var.A1.sim, var.B1.sim, var.AB1.sim, var.A0.sim, var.B0.sim, var.AB0.sim, Pi.sim)
  }
  
  
  
  PHI.LIST.sim <- list() 
  for(i in 1:100){
    PHI.LIST.sim[[i]] = matrix(data = c(MLEs.sim[i,]- theta.sim[j,]),nrow = 11,ncol =1) %*%t(matrix(data = c(MLEs.sim[i,]- theta.sim[j,]),nrow = 11,ncol =1))
  }
  
  #PHI list will contain a list of the calculation of (PHI.jk - PHI.est) - t(PHI.jk - PHI.est)
  PHI.TOTAL.sim = matrix(c(rep(0,121)), nrow =11, ncol = 11) 
  
  for(i in 1:99){
    PHI.TOTAL.sim <- PHI.TOTAL.sim + PHI.LIST.sim[[i]]
    
  }
  
  Cov.Var.Matrix.sim  <- PHI.TOTAL.sim *99/ 100
  colnames(Cov.Var.Matrix.sim) <- c("muA1", "muB1", "muA0", "muB0", "varA1","varB1","covAB1","varA0","varB0","covAB0","Pi")
  rownames(Cov.Var.Matrix.sim) <- c("muA1", "muB1", "muA0", "muB0", "varA1","varB1","covAB1","varA0","varB0","covAB0","Pi")
  
  #Standart error for a0,a1,b0,b1,pi
  SE.sim=rbind(SE.sim,sqrt(diag(Cov.Var.Matrix.sim)))
  
  j=j+1
}

#Average of simulations
mu1.avg = mu1.hold/R
mu0.avg = mu0.hold/R
cov1.avg = cov1.hold/R 
cov0.avg = cov0.hold/R
avg=colMeans(theta.sim)
avg.SE=colMeans(SE.sim)

#Bias
param=c(2,3,0,0,3,1,0,3,1,0,0.3)
bias=avg-param

#Efficiency
V=colVars(theta.sim)
ef=avg.SE/sqrt(V)

#fisher QDA
a0=solve(S0+S1)%*%(mu1-mu0)
a0.sim = solve(cov0.avg+cov1.avg)%*%(mu1.avg-mu0.avg)
Y=t(t(a0)%*%t(X_i[1:100,1:2]))
Y.sim = t(t(a0.sim)%*%t(X_i[1:100,1:2]))
mu_y1 = t(a0)%*%mu1
mu_y1.sim = t(a0.sim)%*%mu1.avg
mu_y0 = t(a0)%*%mu0
mu_y0.sim = t(a0.sim)%*%mu0.avg
sigma_y1 = sqrt(t(a0)%*%S1%*%a0)
sigma_y1.sim = sqrt(t(a0.sim)%*%cov1.avg%*%a0.sim) 
sigma_y0 = sqrt(t(a0)%*%S0%*%a0) 
sigma_y0.sim = sqrt(t(a0.sim)%*%cov0.avg%*%a0.sim) 

#plot of cutoffs
#determin the cutpoint when we set the specificity = 0.8
t_y = 0.2 #FPR = 0.2

c= mu_y0 - qnorm(t_y)*sigma_y0
cutoff1 = t(a0)%*%as.matrix(c(308,299), nrow = 1)
cutoff2 = t(a0)%*%c(315,299)
cutoff3 = t(a0)%*%c(308,290)
cutoff4 = t(a0)%*%c(315,290)
cutoff=c(cutoff3,cutoff4,cutoff1,c,cutoff2)

c.sim= mu_y0.sim - qnorm(t_y)*sigma_y0.sim
cutoff1.sim = t(a0.sim)%*%as.matrix(c(308,299), nrow = 1)
cutoff2.sim = t(a0.sim)%*%c(315,299)
cutoff3.sim = t(a0.sim)%*%c(308,290)
cutoff4.sim = t(a0.sim)%*%c(315,290)
cutoff.sim=c(cutoff3.sim,cutoff4.sim,cutoff1.sim,c.sim,cutoff2.sim)

FPR = 1- pnorm(cutoff, mean=mu_y0,sd=sigma_y0)
TPR = 1- pnorm(cutoff, mean=mu_y1,sd=sigma_y1)

FPR.sim = 1- pnorm(cutoff.sim, mean=mu_y0.sim,sd=sigma_y0.sim)
TPR.sim = 1- pnorm(cutoff.sim, mean=mu_y1.sim,sd=sigma_y1.sim)

#construct the ROC function of t
ROC.t = function(t){
  pnorm(q = (mu_y1-mu_y0)/sigma_y1 + (sigma_y0/sigma_y1)*qnorm(t))
} 

ROC.t.sim = function(t){
  pnorm(q = (mu_y1.sim-mu_y0.sim)/sigma_y1.sim + (sigma_y0.sim/sigma_y1.sim)*qnorm(t))
} 

#calculage the AUC
AUC_y = pnorm((mu_y1-mu_y0)/sqrt(t(a0)%*%(S0+S1)%*%a0)) 
AUC_y.sim = pnorm((mu_y1.sim-mu_y0.sim)/sqrt(t(a0.sim)%*%(cov_0.sim+cov_1.sim)%*%a0.sim)) 

#ROC
plot(FPR.sim, TPR.sim, ylim = c(0,1), xlim = c(0,1),xlab="FPR",ylab="TPR", main = "ROC of simulation and true parameter", pch = 19,col ="red")
curve(ROC.t.sim(x),from = 0, to=1, add = T, col="red")

points(FPR, TPR, ylim = c(0,1), xlim = c(0,1), pch = 19)
curve(ROC.t(x),from = 0, to=1, add = T)

legend("bottomright", legend=c("True ROC", "Simulation"),col=c("black", "red"), lty=1, cex = 0.9,title="ROC Curve", text.font=4, bg='lightblue')


