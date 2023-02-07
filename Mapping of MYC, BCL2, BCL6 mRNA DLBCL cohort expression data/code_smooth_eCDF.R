## from: https://github.com/MichalMarekHoppe/Mapping-of-MYC-BCL2-BCL6-mRNA-DLBCL-cohort-expression-data-into-extent
## author: Choi Hyung Won
## maintained by: Michal Marek Hoppe (mmlhoppe@gmail.com)
## source: https://www.medrxiv.org/content/10.1101/2020.10.20.20216101

# Mapping of MYC, BCL2, BCL6 mRNA DLBCL cohort expression data into % extent

#########################################
### Empirical CDF
### without relying on stepfun format
#########################################
source("required_functions.R")


### Function to map mRNA values to the correpsonding percentile points in IHC scale
mapGE2IHC = function(ihc, mrna, plot=TRUE) {
  
  ### First estimate eCDF's and smooth them
  ecdf.ihc = smooth.eCDF(ihc)
  ecdf.mrna = smooth.eCDF(mrna)

  ### Plot if you like  
  if(plot) {
    par(mfrow=c(1,2))
    plot(ecdf.ihc$x, ecdf.ihc$y, cex=.3)
    lines(ecdf.ihc$x, ecdf.ihc$Smooth, col=2)
    plot(ecdf.mrna$x, ecdf.mrna$y, cex=.3)
    lines(ecdf.mrna$x, ecdf.mrna$Smooth, col=2)
  }
  
  ### Predicted eCDF for observed mRNA values
  pred.mrna = predict.eCDF(mrna)

  ### Find corresponding IHC values by interpolation
  nm = length(mrna)
  ni = nrow(ecdf.ihc)
  predIHC = rep(NA, nm) ### place holder
  
  ### loop through every mRNA measurement
  for(i in 1:nm) {

    nevent = sum(ecdf.ihc$Smooth > pred.mrna[i])
    
    ### if the observed mRNA value is smaller than the first observation point in the eCDF function
    if(nevent == 0) {
      predIHC[i] = ecdf.ihc$x[ni]
    }
    ### if the observed mRNA value is greater than the last observation point in the eCDF function
    if(nevent == ni) {
      predIHC[i] = ecdf.ihc$x[1]
    }
    ### if the observed mRNA value is in between, then we interpolate the value 
    if(nevent > 0 & nevent < ni) {
      ## finding the two neighboring positions on the smoothed eCDF function to interpolate from
      wid = which(ecdf.ihc$Smooth > pred.mrna[i])
      wid = min(wid)
      y1 = ecdf.ihc$Smooth[wid-1]  ## nearest two values
      y2 = ecdf.ihc$Smooth[wid]
      predIHC[i] = ecdf.ihc$x[wid]  ## If CDF values are identical at the two coordinates, end it here
      
      if(y1 != y2) {
        ## if the observed mRNA does not match the exact position on eCDF coordinates
        ## then interpolate
        w1 = abs(pred.mrna[i]-y1) / abs(y2 - y1)
        w2 = abs(y2-pred.mrna[i]) / abs(y2 - y1)
        grid1 = ecdf.ihc$x[wid - 1]
        grid2 = ecdf.ihc$x[wid]
        predIHC[i] = w1 * grid1 + w2 * grid2  #interpolation
      }
    }
  }
  
  res = data.frame(mRNA=mrna, predIHC=predIHC, stringsAsFactors = FALSE)
  res
}


###### Data analysis
library(readr)
idata <- as.data.frame(read_csv("IHC.csv"))
mdata <- as.data.frame(read_csv("mRNA_example.csv"))

### MYC
ihc = idata$MYC
mrna = mdata$MYC
myc.res = mapGE2IHC(ihc, mrna, plot=TRUE)
par(mar = c(1, 1, 1, 1)) #adjust plot margins
par(mfrow=c(2,2))
hist(ihc, breaks=50, main="IHC - MYC")
hist(mrna, breaks=50, main="mRNA - MYC")
plot(myc.res$mRNA, myc.res$predIHC, cex=.2, ylim=c(0,100),
     xlab="Observed", ylab="Transformed")
hist(myc.res$predIHC, breaks=50, main="Transformed - MYC")

### BCL2
ihc = idata$BCL2
mrna = mdata$BCL2
bcl2.res = mapGE2IHC(ihc, mrna, plot=TRUE)
par(mfrow=c(2,2))
hist(ihc, breaks=50, main="IHC - BCL2")
hist(mrna, breaks=50, main="mRNA - BCL2")
plot(bcl2.res$mRNA, bcl2.res$predIHC, cex=.2, ylim=c(0,100),
     xlab="Observed", ylab="Transformed")
hist(bcl2.res$predIHC, breaks=50, main="Transformed - BCL2")

### BCL6
ihc = idata$BCL6
mrna = mdata$BCL6
bcl6.res = mapGE2IHC(ihc, mrna, plot=TRUE)
par(mfrow=c(2,2))
hist(ihc, breaks=50, main="IHC - BCL6")
hist(mrna, breaks=50, main="mRNA - BCL6")
plot(bcl6.res$mRNA, bcl6.res$predIHC, cex=.2, ylim=c(0,100),
     xlab="Observed", ylab="Transformed")
hist(bcl6.res$predIHC, breaks=50, main="Transformed - BCL6")






