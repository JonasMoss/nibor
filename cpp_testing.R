library("Rcpp")
library("RcppArmadillo")
settings=getPlugin("Rcpp")
settings$env$PKG_CXXFLAGS='-std=c++11'

sourceCpp("cpp/testing.cpp")

library('Rcpp')
library('inline')

rcpp_inc <- '
using namespace Rcpp;
using namespace arma;
'

m1 <- matrix(1:16, nr=4)
m2 <- matrix(17:32, nr=4)
v1 <- 1:10
v2 <- 11:20

src <- '
mat m1 = as<mat>(m1in);
mat m2 = as<mat>(m2in);
mat out = join_cols(m1, m2);
return(wrap(out));
'

fn <- cxxfunction(signature(m1in="numeric", m2in="numeric"), src, plugin='RcppArmadillo', rcpp_inc)
res <- fn(m1, m2)
test <- rbind(m1, m2)
all.equal(test, res)
