#include <cmath>
#include <vector>
#include <iostream>
#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace std;

// [[Rcpp::export]]
int cool (int j) {
    auto k = 5.3;
    return k+j;
}
