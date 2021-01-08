/* Contains functions related to the Gaussian copula KDE of Jones & Henderson,
 * specifically the randomizers. */
// [[Rcpp::depends(BH)]]

#include <cmath>
#include <vector>
#include <ctime>
#include <Rcpp.h>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/math/special_functions/erf.hpp>

double probit(double p){
    return sqrt(2)*boost::math::erf_inv(2*p-1);
}


double phi(double x){
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
 
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
 
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return 0.5*(1.0 + sign*y);
}

// [[Rcpp::export]]
Rcpp::NumericVector crnorm(long n, Rcpp::NumericVector mean, 
                           Rcpp::NumericVector sd, int seed) {
  const int maxMean = mean.size();
  const int maxSd = sd.size();
  boost::random::normal_distribution <double> p(0,1);
  boost::mt19937 randGen(seed);
  boost::variate_generator < boost:: mt19937&, boost::random::normal_distribution < double > >
         vargen (randGen, p);
         
  Rcpp::NumericVector variates(n);
  for (int i = 0;i < n;i++){
    variates[i] = vargen()*sd[i%maxSd]+mean[i%maxMean];
  }
  return(variates);
}

// [[Rcpp::export]]
Rcpp::NumericVector crgc(long n, Rcpp::NumericVector X, double rho, int seed) {
  const double mean = 0;
  const double sd   = sqrt(1-rho*rho);
  const int max     = X.size();
  
  for (int i=0; i < max;i++){
    X[i] = rho*probit(X[i]);
  }

  boost::random::normal_distribution < double > normal_variate(mean,sd);
  boost::mt19937 randGen(seed);
  boost::variate_generator < boost::mt19937&, boost::random::normal_distribution < double > >
    vargenNormal (randGen, normal_variate);
  
  Rcpp::NumericVector variates (n);
  
  for (int i = 0;i < n;i++){
    variates[i] = phi(vargenNormal()+X[i%max]);
  }
  
  return(variates);
}

