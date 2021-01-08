new_estimator = function(xx, yy, k=3, maxIter = 100) {
  # A new estimator. Works for k > 2.
  n = length(xx)
  indices = c(1,floor(n/k)*(1:(k-1)),n+1)
  runs = matrix(rep(NA,(10+k+1+1)*maxIter),ncol=10+k+1+1)
  index = NA
  for (i in 1:maxIter) {
    index = sample(setdiff(2:k,index),1)
    xx_new = xx[indices[index-1]:(indices[index+1]-1)]
    yy_new = yy[indices[index-1]:(indices[index+1]-1)]
    run_1d = new_estimator_2cl(xx_new,yy_new)
    indices[index] = run_1d$maximiser$Split + indices[index-1]
    runs[i,] = unlist(c(run_1d$maximiser,index,indices))
  }
  runs
}
