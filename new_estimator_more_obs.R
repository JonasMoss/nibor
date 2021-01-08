new_estimator_1dim = function(xx,yy,bound=3) {
  # A relatively efficient estimator in the 2d-case. 
  # Works by online updating of each LM.
  
  n = length(xx)
  loglik  = rep(0,n-2*bound) 
  logliks = matrix(rep(NA,(n-2*bound)*2),ncol=2)
  beta0   = matrix(rep(NA,(n-2*bound)*2),ncol=2)
  beta1   = matrix(rep(NA,(n-2*bound)*2),ncol=2)
  sigma   = matrix(rep(NA,(n-2*bound)*2),ncol=2)
  
  
  for (i in (1:(n-2*bound))) {
    mod1 = lm(yy[1:(i+bound)] ~ xx[1:(i+bound)])
    mod2 = lm(yy[(i+bound+1):n] ~ xx[(i+bound+1):n])
    loglik[i]   = as.numeric(logLik(mod1) + logLik(mod2))
    logliks[i,] = c(as.numeric(logLik(mod1)),as.numeric(logLik(mod2)))
    beta0[i,]   = c(coef(mod1)[1],coef(mod2)[1])
    beta1[i,]   = c(coef(mod1)[2],coef(mod2)[2])
    sigma[i,]   = c(summary(mod1)$sigma,summary(mod2)$sigma)
  }
  
  which_maximises = which.max(loglik)
  return_coefs = c(beta0 = beta0[which_maximises,],beta1 = beta1[which_maximises,],
                   sigma = sigma[which_maximises,], N = which_maximises + 10, 
                   logliks = logliks[which_maximises,])
  res = return_coefs
}

store = new_estimator_1dim(xx,yy,bound=10)

new_estimator = function(xx, yy, k=3, maxIter = 100) {
  # A new estimator. Works for k > 2.
  n = length(xx)
  indices = c(1,floor(N/k)*(1:(k-1)),N)
  loglik = rep(NA,k)
  runs = matrix(rep(NA,(11+k+1)*maxIter),ncol=11+k+1)
  for (i in 1:maxIter) {
    index = sample(2:k,1)
    xx_new = xx[indices[index-1]:indices[index+1]]
    yy_new = yy[indices[index-1]:indices[index+1]]
    run_1d = new_estimator_1dim(xx_new,yy_new)
    indices[index] = run_1d[7] + indices[index-1]
    loglik[index-1] = run_1d[8]
    loglik[index+1] = run_1d[9]
    runs[i,] = c(run_1d,index,sum(loglik),indices)
  }
  runs = data.frame(runs)
  names(runs) = c("beta01","beta02","beta11","beta12","sigma1","sigma2",
                  "N","loglik_low","loglik_upper","index","loglik")
  list(runs = runs,splits = indices,loglik=sum(loglik))
}


ind1 = sort(sample(1:1000,700))
ind2 = sort(sample(1:1000,700))
ind3 = sort(sample(1:1000,700))
s = new_estimator(xx,yy,k=3,maxIter=25)
s1 = new_estimator(xx[ind1],yy[ind1],k=5,maxIter=10)
s2 = new_estimator(xx[ind2],yy[ind2],k=5,maxIter=10)
s3 = new_estimator(xx[ind3],yy[ind3],k=5,maxIter=10)





new_estimator(xx,yy,k=3,maxIter=50)$loglik
new_estimator(xx,yy,k=4,maxIter=50)$loglik
new_estimator(xx,yy,k=5,maxIter=50)$loglik
new_estimator(xx,yy,k=6,maxIter=50)$loglik
new_estimator(xx,yy,k=7,maxIter=50)$loglik