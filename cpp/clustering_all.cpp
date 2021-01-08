#include <cmath>
#include <vector>
#include <Rcpp.h>
#include <iostream>

using namespace std;

// [[Rcpp::export]]
Rcpp::NumericMatrix cluster_all(Rcpp::NumericVector xx, Rcpp::NumericVector yy, int bound){
    int n = xx.size();
    Rcpp::NumericMatrix value_matrix(n - 2*bound + 1, n - 2*bound + 1);

    
    for (int row_index = 0; row_index < n - 2*bound + 1; row_index++) {
        
        double meanx = 0, deltax_prev = 0, deltax_now = 0, Mx = 0, sigmasq;
        double meany = 0, deltay_prev = 0, deltay_now = 0, My = 0, Cxy = 0;
    
        int index;
        int length;
        
        for (int i = row_index; i < row_index+bound; i++){
            length = i - row_index;
            
            deltax_prev = xx[i] - meanx; // xx[i] - mean(xx[0:(i-1)]
            deltay_prev = yy[i] - meany; // yy[i] - mean(yy[0:(i-1)]
            
            meanx += deltax_prev/(length +1);  // updates the means for x and y
            meany += deltay_prev/(length +1);
            
            deltax_now = xx[i] - meanx;  // xx[i] - mean(xx[0:i])
            deltay_now = yy[i] - meany;  // yy[i] - mean(yy[0:i])
            
            Mx  += deltax_prev*deltax_now;
            My  += deltay_prev*deltay_now;
            Cxy += length*deltax_prev*deltay_prev/(length +1);
        }
    
        sigmasq  = My/bound - Cxy*Cxy/(bound*Mx);             
        value_matrix(row_index,0) = -0.5*bound*(1+log(sigmasq)+log(2*PI));  

        for (int j = 1; j <= (n - 2*bound - row_index); j++){
            index = j + bound - 1 + row_index;
            length = j + bound - 1;
        
            deltax_prev = xx[index] - meanx;
            deltay_prev = yy[index] - meany; 
        
            meanx += deltax_prev/(length+1); 
            meany += deltay_prev/(length+1);
        
            deltax_now = xx[index] - meanx; 
            deltay_now = yy[index] - meany; 
        
            Mx  += deltax_prev*deltax_now;
            My  += deltay_prev*deltay_now;
            Cxy += length*deltax_prev*deltay_prev/(length +1);
        
            sigmasq = My/(length +1) - Cxy*Cxy/((length+1)*Mx);          
            value_matrix(row_index,j) = -0.5*(length+1)*(1+log(sigmasq)+log(2*PI));
        }
        
    }
    
    return(value_matrix);

}
