library("Rcpp")
sourceCpp("cpp/clustering_all_2.cpp")

bound = 3
df_ = cluster_all(xx,yy,bound)
df = cluster_2cl(xx,yy,bound)

head(cbind(df_[1,],df[,5]))
df_[3,]
cluster_2cl(xx[3:N],yy[3:N],bound)[,5]

optim_it = function(i,j) {
  lower = c(0,df_[,j+1])
  which.max(df_[i,] + lower)
}

plot(df_[1,] + lower)

lower = df_[i,]
upper = c(df_[,j][-1],NA)
sapply(bound:bound:(50-bound),function(i) lower[i] + upper[i+1])
plot(lower + upper)
points(10:(dim(df)[1]+10-1),totliks,col="red")
logLik(lm(yy[1:50]~xx[1:50]))

df2 = rev(cluster_2cl(rev(xx),rev(yy),bound)[,5])