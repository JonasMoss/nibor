tt = seq(0.01,10,by=0.01)
N = length(tt)
set.seed(17932561417)
splits = c(1,sort(sample(2:1000,2)),1001)
split_lengths = sapply(1:3,function(i) splits[i+1] - splits[i])
split_lengths

beta0s = runif(3,-3,3)
beta0s = unlist(sapply(1:3,function(i) rep(beta0s[i],split_lengths[i])))

beta1s = runif(3,1,7)
beta1s = unlist(sapply(1:3,function(i) rep(beta1s[i],split_lengths[i])))

sigmas = runif(3,1,6)
sigmas = unlist(sapply(1:3,function(i) rep(sigmas[i],split_lengths[i])))

mus = runif(3,-10,10)
mus = unlist(sapply(1:3,function(i) rep(mus[i],split_lengths[i])))

mu_sigmas = runif(3,0,6)
mu_sigmas = unlist(sapply(1:3,function(i) rep(mu_sigmas[i],split_lengths[i])))

xx = rnorm(N,mus,mu_sigmas)
yy = beta0s + xx*beta1s + rnorm(N,0,sigmas)
plot(xx,yy)

ggplot(data.frame(xx),aes(x = xx, y = yy, colour = tt)) +
  geom_point(shape = 18,size=1.6) +    
  geom_smooth(method = "lm",col = "purple",size = 0.5) +
  scale_colour_gradient2(low="red", high="blue", mid="pink",
                         midpoint = median(tt), name = "Time")