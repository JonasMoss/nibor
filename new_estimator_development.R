new_estimator_2cl = function(xx, yy, bound = 10) {
  n = length(xx)

  # Housekeeping variables.
  xx1 = xx[1:(bound)]
  xx2 = xx[(bound+1):n]
  
  yy1 = yy[1:(bound)] 
  yy2 = yy[(bound+1):n]

  sumx1 = sum(xx1)
  sumx2 = sum(xx2)

  sumy1 = sum(yy1)
  sumy2 = sum(yy2)

  sumxy1 = sum(xx1*yy1)
  sumxy2 = sum(xx2*yy2)

  sumxsq1 = sum(xx1^2)
  sumxsq2 = sum(xx2^2)

  sumysq1 = sum(yy1^2)
  sumysq2 = sum(yy2^2)

  # Variables of interest outside the scope of this function.
  n1 = (bound)
  n2 = n - (bound)

  beta11 = cov(xx1,yy1)/var(xx1)
  beta12 = cov(xx2,yy2)/var(xx2)

  beta01 = mean(yy1) - beta11*mean(xx1)
  beta02 = mean(yy2) - beta12*mean(xx2) 

  sigmasq1 = (sumysq1 - 2*beta01*sumy1 - 2*beta11*sumxy1 + 
              sumx1*2*beta01*beta11 + beta11^2*sumxsq1 + n1*beta01^2)/n1
  
  sigmasq2 = (sumysq2 - 2*beta02*sumy2 - 2*beta12*sumxy2 + 
              sumx2*2*beta02*beta12 + beta12^2*sumxsq2 + n2*beta02^2)/n2

  loglik1 = -1/2*(bound)*(1+log(sigmasq1)+log(2*pi))
  loglik2 = -1/2*(n-bound)*(1+log(sigmasq2)+log(2*pi))

  score = (loglik1 + loglik2)/n

  # Definition of matrices for house keeping and what's more.

  coef_df = data.frame(matrix(rep(NA,(n-2*bound)*10),ncol=10))
  names(coef_df) = c("Beta0_lower","Beta1_lower","Sigma_lower","Loglik_lower",
                     "Beta0_upper","Beta1_upper","Sigma_upper","Loglik_upper",
                     "Split","Score")

  coef_df[1,] = c(beta01,beta11,sqrt(sigmasq1),loglik1,
                beta02,beta12,sqrt(sigmasq2),loglik2,
                bound+1,score)


  # Updating functions. 
  update_up = function(j){
    x = xx[j]
    y = yy[j]
    n1         <<- n1 + 1
    sumx1      <<- sumx1 + x
    sumy1      <<- sumy1 + y
    sumxsq1    <<- sumxsq1 + x^2 
    sumysq1    <<- sumysq1 + y^2
    sumxy1     <<- sumxy1 + x*y
    covXY      = 1/n1*(sumxy1 - 1/n1*sumx1*sumy1)
    varX       = 1/n1*(sumxsq1 - 1/n1*sumx1^2)
    meanX      = 1/n1*sumx1
    meanY      = 1/n1*sumy1
    beta1      = covXY/varX
    beta0      = meanY - beta1*meanX
    sigmasq    = 1/n1*(sumysq1 - 2*beta0*sumy1 - 
                       2*beta1*sumxy1 + sumx1*2*beta0*beta1 + 
                       beta1^2*sumxsq1 + n1*beta0^2)
    c(beta0 = beta0, beta1 = beta1, sigma = sqrt(sigmasq), 
              loglik = -1/2*n1*(1+log(sigmasq)+log(2*pi)))
  }

  update_down = function(j){
    x = xx[j]
    y = yy[j]
    n2         <<- n2 - 1
    sumx2      <<- sumx2 - x
    sumy2      <<- sumy2 - y
    sumxsq2    <<- sumxsq2 - x^2 
    sumysq2    <<- sumysq2 - y^2
    sumxy2     <<- sumxy2 - x*y
    covXY      = 1/n2*(sumxy2 - 1/n2*sumx2*sumy2)
    varX       = 1/n2*(sumxsq2 - 1/n2*sumx2^2)
    meanX      = 1/n2*sumx2
    meanY      = 1/n2*sumy2
    beta1      = covXY/varX
    beta0      = meanY - beta1*meanX
    sigmasq    = 1/n2*(sumysq2 - 2*beta0*sumy2 - 
                       2*beta1*sumxy2 + sumx2*2*beta0*beta1 + 
                       beta1^2*sumxsq2 + n2*beta0^2)
    c(beta0 = beta0, beta1 = beta1, sigmasq = sqrt(sigmasq), 
      loglik = -1/2*n2*(1+log(sigmasq)+log(2*pi)))
  }

  for (j in (1:(n-2*bound))) {
    coef_df[j+1,] = c(update_up(bound+j),update_down(bound+j),j+bound+1,NA)
  }

  coef_df$Score = 1/n*(coef_df$Loglik_lower + coef_df$Loglik_upper)

  coef_df[which.max(coef_df$Score),]
  
  return(list(coefdf = coef_df,maximiser = coef_df[which.max(coef_df$Score),]))
}

library(ggplot2)

ggplot(coef_df,aes(x = Score, y = Sigma_upper, colour = Split)) +
  geom_point(shape = 18,size=1.6) +    
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint = median(coef_df$Split), name = "Splits")

ggplot(coef_df,aes(x = Score, y = Sigma_lower, colour = Split)) +
  geom_point(shape = 18,size=1.6) +    
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint = median(coef_df$Split), name = "Splits")


ggplot(coef_df,aes(x = Split, y = Score, colour = Sigma_lower)) +
  geom_point(shape = 18,size=1.6) +    
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint = median(coef_df$Sigma_lower), name = "Sigma")