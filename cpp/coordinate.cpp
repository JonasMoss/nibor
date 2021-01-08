#include <cmath>
#include <vector>
#include <Rcpp.h>

using namespace std;

/* Functions for combos of real and l2*. These pointers are needed in order to
 * avoid many repeated if calls and make the code more readable. */

double (*residual)(int,int,double*,int,double,double);
double (*residual_final)(int,double*,int,double,double);

/* These correspond to Kullback-Leibler and equal weights. */
double rkf(int i,int j,double* data,int len,double k,double lim){
  double a = -(j-i);
  double b = data[j]-data[i];
  if (b < lim) return -0.1/0.0;
  else return(a*log(b));
}

double rkf_final(int i,double* data,int len,double k,double lim){
  double a = -(len-i);
  double b = 1-data[i];
  if (b < lim) return -0.1/0.0;
  else return(a*log(b));
}

/* These correspond to L2 and equal weights. */    
double rlf(int i,int j,double* data,int len,double k,double lim){
  double a = 2*(j-i)/((double) len)-1/k;
  double b = data[j]-data[i];
  if (b < lim) return -0.1/0.0;
  else return(a/b);
}

double rlf_final(int i,double* data,int len,double k,double lim){
  double a = 2*(len-i)/((double) len)-1/k;
  double b = 1-data[i];
  if (b < lim) return -0.1/0.0;
  else return(a/b);
}
    
/* ... and these correspond to KL weights and splits! */
double rkt(int i,int j,double* data,int len,double k,double lim){
  double b = data[j]-data[i];
  if (b < lim) return -0.1/0.0;
  else return((j-i)*(log(j-i)-log(len)-log(b)));
}

double rkt_final(int i,double* data,int len,double k,double lim){
  double b = 1-data[i];
  if (b < lim) return -0.1/0.0;
  else return((len-i)*(log(len-i)-log(len)-log(b)));
}

/* Finally, L2 weights and splits. */
double rlt(int i,int j,double* data,int len,double k,double lim){
  double a = pow(j-i,2);
  double b = data[j]-data[i];
  if (b < lim) return -0.1/0.0;
  else return(a/b);
}

double rlt_final(int i,double* data,int len,double k,double lim){
  double a = pow(1-i,2);
  double b = 1-data[i];
  if (b < lim) return -0.1/0.0;
  else return(a/b);
}


// [[Rcpp::export]]
Rcpp::NumericVector cpp_greedy_alg(bool real_hist,bool l2,Rcpp::NumericVector data,
                                   int len,double k,int modulator,Rcpp::NumericVector init,
                                   double lim){
    
    /* residual and residual_final point to the functions needed in maximization. 
     * The definition of these functions vary as L2 and real varies, every other aspect
     * of the algorithm stays constant.
     */
    
    if (real_hist) {
        if (l2) {
            residual = &rlt;
            residual_final = &rlt_final;         
        }
        else {
            residual = &rkt;
            residual_final = &rkt_final;
        }
    }

    else {
        if (l2) {
            residual = &rlf;
            residual_final = &rlf_final;
         
        }
        else {
            residual = &rkf;
            residual_final = &rkf_final;
        }    
    }   

    /* We define the pretty "matrix" estimatse indices, with the correct dimensions.
     * This matrix contains the ML / L2 estimate indices, as j-ary vectors. The (i,j)-th 
     * element corresponds to ML-estimate with k=j and data[0:j]. */
    
    vector <int> estimates;
    vector <int> test_estimates;
    
    estimates.resize((int) k+1);
    test_estimates.resize((int) k+1);
    estimates[0] = 0;
    estimates[k] = 1;
    
    for(int i = 1;i<k;i++){
      estimates[i] = init[i-1];
    }
    
    int ended = 0;
    int over;
    int under;
    double max;
    double temp;
    for (int j = 0;j<modulator*k;j++){
      
      test_estimates = estimates;
      
      /* The loop takes care of all the values except the final. */
      
      for (int i=1;i<(k-1);i++){
        over = estimates[i+1];
        under = estimates[i-1];
        max = -0.1/0.0;
        for (int p = under; p<over;p++){
          if (data[p]-data[under]>lim && data[over]-data[p]>lim){
            temp = residual(under,p,data.begin(),len,k,lim)+residual(p,over,data.begin(),len,k,lim);
          }
          else temp = -0.1/0.0;
          if (max < temp) {
            max = temp;
            estimates[i] = p;
          }
        }
        
      }
      
      /* And now is the time for the last value. */
      int i = (k-1);
      under = estimates[i-1];
      max = -0.1/0.0;
      for (int p = under; p<(len+1);p++){
        if (data[p]-data[under]>lim && 1-data[p]>lim){
            temp = residual(under,p,data.begin(),len,k,lim)+residual_final(p,data.begin(),len,k,lim);
          }
        else temp = -0.1/0.0;
        if (max < temp) {
            max = temp;
            estimates[i] = p;
          }
      }
      
      /* We test the break condition. */
      
      if (test_estimates == estimates){
        ended = j;
        break;
      }
      
    }
    
    Rcpp::NumericVector xx((int) k);
  
    for (int i=0;i<k-1;i++){
        xx[i] = estimates[i+1];
      }
      
    xx[k-1] = ended;

    return(xx);

}

