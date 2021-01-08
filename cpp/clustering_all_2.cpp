#include <cmath>
#include <vector>
#include <Rcpp.h>
#include <iostream>

using namespace std;

// [[Rcpp::export]]
Rcpp::List cluster_all(Rcpp::NumericVector xx, Rcpp::NumericVector yy, int bound){
    int n = xx.size();
    
    Rcpp::NumericMatrix value_matrix(n, n, Rcpp::NumericVector::get_na());
    Rcpp::List return_list;

    for (int row_index = 0; row_index < (n-1); row_index++) {
        
        double meanx = 0, deltax_prev = 0, deltax_now = 0, Mx = 0, sigmasq;
        double meany = 0, deltay_prev = 0, deltay_now = 0, My = 0, Cxy = 0;
    
        int index;
        int length = 0;
            
        for (int j = 0; j < (bound-1); j++){
            index = j + row_index;
            length = j;
        
            deltax_prev = xx[index] - meanx;
            deltay_prev = yy[index] - meany; 
        
            meanx += deltax_prev/(length+1); 
            meany += deltay_prev/(length+1);
        
            deltax_now = xx[index] - meanx; 
            deltay_now = yy[index] - meany; 
        
            Mx  += deltax_prev*deltax_now;
            My  += deltay_prev*deltay_now;
            Cxy += length*deltax_prev*deltay_prev/(length+1);
        }
        
        for (int j = (bound-1); j < (n - row_index); j++){
            index = j + row_index;
            length = j;
        
            deltax_prev = xx[index] - meanx;
            deltay_prev = yy[index] - meany; 
        
            meanx += deltax_prev/(length+1); 
            meany += deltay_prev/(length+1);
        
            deltax_now = xx[index] - meanx; 
            deltay_now = yy[index] - meany; 
        
            Mx  += deltax_prev*deltax_now;
            My  += deltay_prev*deltay_now;
            Cxy += length*deltax_prev*deltay_prev/(length +1);
        
            sigmasq = My/(length+1) - Cxy*Cxy/((length+1)*Mx);          
            value_matrix(row_index,j+row_index) = -0.5*(length+1)*(1+log(sigmasq)+log(2*PI));
        }
        
    }
    
    return_list = Rcpp::List::create(value_matrix, 
                                     value_matrix);
    return(return_list);

}
