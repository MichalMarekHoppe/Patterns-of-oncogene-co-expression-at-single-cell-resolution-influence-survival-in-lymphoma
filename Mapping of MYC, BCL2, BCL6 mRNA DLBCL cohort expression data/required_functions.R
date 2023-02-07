
############################
## Kernel smoother function
## (simpler than LOESS)

gauss.kernel.smooth = function(x, y, newx=x, k.width) {
  
  ## x and y are the observed data
  ## newx is provided at grid points of x as the user wishes
  
  yest = y
  n = length(y)
  min.x = min(x)
  max.x = max(x)
  
  ### Iterate through each data point
  for(i in 1:n) {
    ### Weight coefficients for all data points with respect to x_i
    wt = (x - x[i]) / k.width  
    wt = dnorm(wt, 0, 1) 
    ### Weighted average at x_i
    yest[i] = sum(wt * y, na.rm=TRUE) / sum(wt, na.rm=TRUE)
  }
  fitTab = data.frame(x=x, y=yest)
  
  ### for prediction
  nx = length(newx)
  ypred = rep(NA, nx)
  for(i in 1:nx) {
    wt = (x - newx[i]) / k.width 
    wt = dnorm(wt, 0, 1)
    ypred[i] = sum(wt * y, na.rm=TRUE) / sum(wt, na.rm=TRUE)
  }
  ypred[newx < min.x] = 0
  ypred[newx >= max.x] = 1
  predTab = data.frame(x=newx, y=ypred)
  
  list(fit=fitTab, pred=predTab)  
  
}


############################
## Function to estimate empirical CDF
## -- not using stepfun environment as the built-in
## This function will no longer contain the original data
## but keeps the eCDF at unique, ordered coordinates of X

eCDF = function(x) {
  x = x[!is.na(x)]
  xx = x[order(x)]
  uxx = unique(xx)  ## ordered unique coordinates
  nxx = length(uxx)
  ecdf = rep(NA, nxx)
  for(i in 1:nxx) {
    ecdf[i] = mean(xx <= uxx[i])
  }
  res = data.frame(x=uxx, y=ecdf)
  res
}

############################
## Smoother for the eCDF

smooth.eCDF = function(X) {
  ### X is the eCDF object from above
  X = X[!is.na(X)]
  rg = range(X)
  wd = diff(rg) / 100
  x.cdf = eCDF(X)
  x.fit = gauss.kernel.smooth(x=x.cdf$x, y=x.cdf$y, newx=x.cdf$x, k.width=wd*2)
  x.cdf$Smooth = x.fit$pred$y
  x.cdf  ### returns augmented eCDF table
}

############################
## Wrapper to report eCDF at 
## the raw data coordinates

predict.eCDF = function(X) {
  ### X is the eCDF object from above
  X = X[!is.na(X)]
  rg = range(X)
  wd = diff(rg) / 100
  x.cdf = eCDF(X)
  x.fit = gauss.kernel.smooth(x=x.cdf$x, y=x.cdf$y, newx=X, k.width=wd*2)
  x.fit$pred$y   ### returns predicted eCDF table at original data points
}




