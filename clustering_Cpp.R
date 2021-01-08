library("Rcpp")
sourceCpp("cpp/clustering_2.cpp")

yy = ndf$Spread
xx = ndf$Median

n = length(xx)
bound = 10

df = cluster_2cl(xx,yy,bound)
df2 = cluster_2cl(rev(xx),rev(yy),bound)

totliks = df[,5] + rev(df2[,5])
plot(bound:(dim(df)[1]+bound-1),totliks)
abline(v = which.max(totliks) + bound)

new_estimator_2cl(xx,yy,bound=10)

totsigma= df[,3] + rev(df2[,3])
plot(bound:(dim(df)[1]+bound-1),totsigma)

colnames(df) = c("beta1","beta0","sigmasq","rsq","loglik")

head(df)
sapply(3:16,function(i) logLik(lm(yy[1:i] ~ xx[1:i])))
sapply(3:16,function(i) summary(lm(yy[1:i] ~ xx[1:i]))$r.squared)

plot(1:dim(df)[1],df[,1])
plot(1:dim(df)[1],df[,2])
plot(1:dim(df)[1],df[,3])
plot(1:dim(df)[1],df[,4])
plot(1:dim(df)[1],df[,5])

tail(df)
logLik(lm(yy[1:990]~xx[1:990]))

head(totliks)
logLik((lm(yy[1:11]~xx[1:11]))) + logLik((lm(yy[12:n]~xx[12:n])))
