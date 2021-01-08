// [[Rcpp::export]]
long fibc1(const int x) {
  if (x == 0) return(0); 
  if (x == 1) return(1);
  return (fibc1(x - 1)) + fibc1(x - 2);
}

// [[Rcpp::export]]
long fibc2(const int n) {
  long val1 = 0;
  long val2 = 1;
  long temp = 0;
  
  for (int i = 0;i < n;i++){
    temp = val2;
    val2 = val1 + val2;
    val1 = temp;
  }
  
  return (val1);
}