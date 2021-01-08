new_estimator = function(xx,yy) {
  n = length(xx)
  loglik = rep(0,n-20) 
  beta0 = matrix(rep(NA,(n-20)*2),ncol=2)
  beta1 = matrix(rep(NA,(n-20)*2),ncol=2)
  sigma = matrix(rep(NA,(n-20)*2),ncol=2)
  for (i in (1:(n-20))) {
    mod1 = lm(yy[1:(i+10)] ~ xx[1:(i+10)])
    mod2 = lm(yy[(i+10+1):n] ~ xx[(i+10+1):n])
    loglik[i] = as.numeric(logLik(mod1) + logLik(mod2))
    beta0[i,] = c(coef(mod1)[1],coef(mod2)[1])
    beta1[i,] = c(coef(mod1)[2],coef(mod2)[2])
    sigma[i,] = c(summary(mod1)$sigma,summary(mod2)$sigma)
  }
  
  which_maximises = which.max(loglik)
  return_coefs = c(beta0 = beta0[which_maximises,],beta1 = beta1[which_maximises,],
                   sigma = sigma[which_maximises,], N = which_maximises + 10)
  res = list(loglik = loglik, beta0 = beta0, beta1 = beta1, 
             sigma = sigma, return_coefs = return_coefs)
}

store = new_estimator(xx,yy)


ggplot(data.frame(store$loglik),aes(x = store$loglik, y = store$sigma[,2], colour = tt[10:989])) +
  geom_point(shape = 18,size=1.6) +    
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint = median(tt[10:989]), name = "Time")