#include <cmath>
#include <vector>
#include <Rcpp.h>
#include <iostream>

using namespace std;


// [[Rcpp::export]]
int scpp_exact_alg(bool real_hist,bool l2,Rcpp::NumericVector data,int len,double k,double lim){
    
    vector < vector < vector < int > > > estimates;
    estimates.resize(len+2);
    for (int i = 0; i <= len+1; i++){
        estimates[i].resize(k);
        for (int j = 0;j<(k-1);j++){
            estimates[i][j].resize(j+1);
        }
    }

    
    vector < vector < double > > objective;
    objective.resize(len+2);
    for (int i = 0; i <= len+1; i++){
        objective[i].resize(k-1);
    }
    
    vector <int> est;
    est.resize(1);

    return(0);

}
