sim.2d.data = function(n,beta,rho,f=0.2){
  X = rbinom(n,2,f)
  E = rbinorm(n,rho)
  Y = t(t(t(X)) %*% t(beta)) + E
  return(list(Y=Y,X=X))
}

#simulate independent 5-d data
sim.5d.data = function(n,beta,f=0.2){
  X = rbinom(n,2,f)
  Sigma = diag(1,5)
  mu = rep(0,5)
  E = rmvnorm(n,mu,Sigma)
  Y = t(t(t(X)) %*% t(beta)) + E
  return(list(Y=Y,X=X))
}

#simulate data with a single underlying "factor" lambda
#model is Y (d by n) = Lambda (d by k) F(k by n) + e (d by n)
#and rows of F are N(0,1)+beta X
#beta is a k-vector (k by 1), X is an n-vector (n by 1)
#and e_d,n is N(0,sigma^2_d)
sim.data.factor = function(n,beta,Lambda,sigma,f=0.2){
  X = rbinom(n,2,f)
  beta = cbind(beta)
  X = cbind(X)
  k = ncol(Lambda)
  F = rmvnorm(n,beta %*% t(X),diag(1,k))
  Y = Lambda %*% F + rmvnorm(n,0,diag(sigma^2))
  return(list(Y=Y,X=X,F=F))
}

#does nrep 2d simulations and applies multivariate association test to them
simset=function(n,beta,rho,f=0.2,nrep=100){
s = list()
t.04=list()
t.01=list()
m=list()
for(i in 1:nrep){
  s[[i]] = sim.2d.data(n,beta,rho,f)
  t.04[[i]] = logBF.rankone.matrix(s[[i]]$X,t(s[[i]]$Y),sigmaa=0.4)
  t.01[[i]] = logBF.rankone.matrix(s[[i]]$X,t(s[[i]]$Y),sigmaa=0.1)
  m[[i]] = summary(manova(t(s[[i]]$Y) ~ s[[i]]$X))$stats[1,6]
}
return(list(s=s,t.01=t.01,t.04=t.04,m=m))
}

#does nrep 5d simulations and applies multivariate association test to them
simset5=function(nrep,n,beta,f=0.2){
s = list()
t.04=list()
t.01=list()
pca=list()
m=list()
for(i in 1:nrep){
  #s[[i]] = sim.5d.data(n,beta,f)
  s[[i]] = sim.data.factor(n,beta,cbind(c(1,1,1,1,1),c(2,-2,1,-1,3)),c(1,1,1,1,1))
  t.04[[i]] = logBF.rankone.matrix(s[[i]]$X,t(s[[i]]$Y),sigmaa=0.4)
  pca[[i]] = prcomp.test(s[[i]]$X,t(s[[i]]$Y),0.4)
  t.01[[i]] = logBF.rankone.matrix(s[[i]]$X,t(s[[i]]$Y),sigmaa=0.1)
  #m[[i]] = summary(manova(t(s[[i]]$Y) ~ s[[i]]$X))$stats[1,6]
}
return(list(s=s,t.04=t.04,t.01=t.01,pca=pca))
}

