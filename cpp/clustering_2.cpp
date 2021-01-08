#include <cmath>
#include <vector>
#include <Rcpp.h>
#include <iostream>

using namespace std;

// [[Rcpp::export]]
Rcpp::NumericMatrix cluster_2cl(Rcpp::NumericVector xx, Rcpp::NumericVector yy, int bound){
    /* Calculates the ML splits, coefficients, loglikelihood, etc. when k = 2.
     * It's the foundation of the general DP algorithm. "bound" gives the
     * minimal number of elements in one cluster. 
     * 
     * n = length of xx
     * i is connected to xx[0,i-1], a proposed number of elements belonging to the first
     * 
     * The algorithm works iteratively and onliny. First it calculates the
     * coefficients for i = 10. Then it exhausts all other possibilities, up to n - i.
     * 
     * The function returns a (n-2*bounds) x 12-matrix containg the following columns:
     * Beta0_lower Beta1_lower Sigma_lower Loglik_lower Rsq_lower ...
     * Beta0_upper Beta1_upper Sigma_upper Loglik_upper Rsq_upper ...
     * Split Score
     */
    
    int n = xx.size();
    Rcpp::NumericMatrix value_matrix(n - 2*bound + 1,5);
    
    /*
    vector < vector < double > > value_matrix;
    value_matrix.resize(rows);
    for (int i = 0; i < rows; i++){
        value_matrix[i].resize(columns);
    } 
    */

    /* A numerically stable online algorithm for the computation of the
     * covariance and variance. See the webpage
     * http://prod.sandia.gov/techlib/access-control.cgi/2008/086212.pdf 
     * (Philippe Pebay, 2008).*/
     
    /* This loop initialises the values meanx, meany, Mx, My, an Cxy,
     * if there are no initial values provided to the function. */
     
    double meanx = 0, deltax_prev = 0, deltax_now = 0, Mx = 0;
    double meany = 0, deltay_prev = 0, deltay_now = 0, My = 0, Cxy = 0;
    
    for (int i = 0; i < bound; i++){
        deltax_prev = xx[i] - meanx; // xx[i] - mean(xx[0:(i-1)]
        deltay_prev = yy[i] - meany; // yy[i] - mean(yy[0:(i-1)]
        
        meanx += deltax_prev/(i+1);  // updates the means for x and y
        meany += deltay_prev/(i+1);
        
        deltax_now = xx[i] - meanx;  // xx[i] - mean(xx[0:i])
        deltay_now = yy[i] - meany;  // yy[i] - mean(yy[0:i])
        
        Mx  += deltax_prev*deltax_now;
        My  += deltay_prev*deltay_now;
        Cxy += i*deltax_prev*deltay_prev/(i+1);
    }
    
    value_matrix(0,0) = Cxy/Mx;                                                // beta1             
    value_matrix(0,1) = meany - value_matrix(0,0)*meanx;                       // beta0
    value_matrix(0,2) = My/bound - Cxy*Cxy/(bound*Mx);                         // sigmasq
    value_matrix(0,3) = Cxy*Cxy/(Mx*My);                                       // rsq
    value_matrix(0,4) = -0.5*bound*(1+log(value_matrix(0,2))+log(2*PI));       // loglik
    
    int i;
    for (int j = 1; j <= (n - 2*bound); j++){
        i = j + bound - 1;
        
        deltax_prev = xx[i] - meanx;
        deltay_prev = yy[i] - meany; 
        
        meanx += deltax_prev/(i+1); 
        meany += deltay_prev/(i+1);
        
        deltax_now = xx[i] - meanx; 
        deltay_now = yy[i] - meany; 
        
        Mx  += deltax_prev*deltax_now;
        My  += deltay_prev*deltay_now;
        Cxy += i*deltax_prev*deltay_prev/(i+1);
        
        value_matrix(j,0) = Cxy/Mx;                                                // beta1             
        value_matrix(j,1) = meany - value_matrix(j,0)*meanx;                       // beta0
        value_matrix(j,2) = My/(i+1) - Cxy*Cxy/((i+1)*Mx);                         // sigmasq
        value_matrix(j,3) = Cxy*Cxy/(Mx*My);                                       // rsq
        value_matrix(j,4) = -0.5*(i+1)*(1+log(value_matrix(j,2))+log(2*PI));       // loglik
    }
    
    
    
    return(value_matrix);

}
