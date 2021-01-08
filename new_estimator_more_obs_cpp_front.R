new_estimator = function(xx, yy, k=4, maxIter = 100, bound=10) {
  # A new estimator. Works for k > 2.
  n = length(xx)
  
  indices = c(1,floor(n/k)*(1:(k-1)),n+1)
  splits = matrix(rep(NA,(k+1)*maxIter),ncol=(k+1))
  sumlogliks = rep(NA,maxIter)
  
  index = NA
  loglik = rep(NA,k)
  run_order = 1:(k-1)
  
  for (i in 1:(maxIter)) {
    
    for (index in run_order) {
      xx_new = xx[indices[index]:(indices[index+2]-1)]
      yy_new = yy[indices[index]:(indices[index+2]-1)]
      
      run_1 = cluster_2cl(xx_new, yy_new, bound)
      run_2 = cluster_2cl(rev(xx_new),rev(yy_new),bound)
      
      best_index = which.max(run_1[,5] + rev(run_2[,5]))
      indices[index+1] = best_index - 1 + bound + indices[index]
      
      loglik[index] = run_1[best_index,5]
      loglik[index+1] = rev(run_2[,5])[best_index]
    }
    
    splits[i,] = indices
    sumlogliks[i] = sum(loglik)
    run_order = sample(run_order)
    
  }
  best_run = which.max(sumlogliks)
  return(list(splits = splits[best_run,], loglikelihood = sumlogliks[best_run]))
}

xx = rnorm(1000)
yy = 5*xx + 1 + rnorm(1000*xx^2)
new_estimator(xx,yy,k = 10, maxIter = 10, bound = 10)


yy = ndf$Spread
xx = ndf$Median


ks = 3:60
logliks = sapply(ks,function(k) new_estimator(xx,yy,k=k)$loglikelihood)
AICs = 2*logliks - 2*(ks + 3)
BICs = 2*logliks - log(length(xx))*(ks + 3)
plot(ks,logliks)
#plot(ks,AICs)
#plot(ks,BICs)