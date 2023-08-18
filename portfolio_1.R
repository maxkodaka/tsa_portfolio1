library(raster)
library(caTools)
library(Kendall)

filename = 'GIMMS_Iberian'
filename_class = 'clc8km_clip_recode_Iberian_2classes'
#setwd('/home/maxim/Documents/coursework/time-series-analysis/data/')
setwd('U:\\Time Series Analysis\\s5')
cube = read.ENVI(filename,headerfile=paste(filename,'.hdr',sep='')) #from caTools
luc<-read.ENVI(filename,headerfile=paste(filename_class,'.hdr',sep='')) #

dim(cube)
dim(luc)
years = dim(cube)[3]/12
lin = 1:12
# Fourier Polynomials
x1 = sin(2*pi*1/12*(lin))
x2 = cos(2*pi*1/12*(lin))
x3 = sin(2*pi*2/12*(lin))
x4 = cos(2*pi*2/12*(lin))
fp = cbind(lin,x1,x2,x3,x4)
cor(fp)

nrows=nrow(cube)
ncols=ncol(cube)

# Initialize all arrays
dim = c(years,ncol(cube),nrow(cube))
t = array(NA,dim)                       # dimensions: 32, 152, 93 This is not used for anything??
meanNDVI = aperm(t)                     # dimensions: 93, 152, 32

magnitude = array(NA, 32)                   # dimensions: 32      Do we still need to initialize if its so small now?
peaktime = array(NA, 32)                    # dimensions: 32      Do we still need to initialize if its so small now?
integral = array(NA, 32)                    # dimensions: 32      Do we still need to initialize if its so small now?

# Maybe we should not allocate these because each pixel has a list of 5 model attributes.
# At least this initialization does not work.
#model_integral = array(NA, c(nrows, ncols)) # dimensions: 93, 152
#model_amp = array(NA, c(nrows, ncols))      # dimensions: 93, 152
#model_peak = array(NA, c(nrows, ncols))     # dimensions: 93, 152

pintegral = array(NA, 32) 
pamp = array(NA, 32)      
ppeak = array(NA, 32)   

syndrome<-matrix(data <- -999, nrow <- dim(cube)[1],ncol <- dim(cube)[2])

## Some helpful functions

## increasing
incr <- function(x, sig){
  abs(x) < sig & x > 0
}

## decreasing
decr <- function(x, sig){
  abs(x) < sig & x < 0
}

## no change
nought <- function(x, sig){
  abs(x) > sig
}

pb = txtProgressBar(title='progress bar',min=0,max =nrows) #run together with loop
for (i in 1:nrows) {
  for (j in 1:ncols) {
    setTxtProgressBar(pb, i, title = round(i / nrow(cube) * 100, 0))
    GIMMS = cube[i, j,]
    
    if (var(GIMMS) > 0) {
      time = 1:12
      for (k in 1:years) {
        # Fit Fourier Polynomials to monthly NDVI data derive the following phenological
        # properties in each year: a) NDVI peaking times, b) NDVI magnitudes, c) NDVI integrals.
        model = lm(GIMMS[time] ~ fp)
        ##################################################### The following have been changed to single vectors because:
        magnitude[k] = max(model$fit) - min(model$fit)      # We may as well use these as temporary variables?
        peaktime[k] = which(model$fit == max(model$fit))    # We may as well use these as temporary variables?
        integral[k] = sum(GIMMS[time])                      # We may as well use these as temporary variables?
        
        meanNDVI[i, j, k] = mean(GIMMS[time])
        # Fit a trend model (Seasonal Kendall test) to derive monotonic trends in the three phenological time-series.
        
        time = time + 12
      }
      
      model_integral = MannKendall(ts(integral))     # Should be per pixel? Or just temporary variables
      model_amp = MannKendall(ts(magnitude))         # Should be per pixel? Or just temporary variables
      model_peak = MannKendall(ts(peaktime))         # Should be per pixel Or just temporary variables
      
      # Need to dissect the previous results to figure out what to do here:
      # Content of the MannKendall object: tau, sl, S, D, varS
      
      pintegral = model_integral$sl * as.numeric(sign(model_integral$tau))
      pamp = model_amp$sl * as.numeric(sign(model_amp$tau))
      ppeak = model_peak$sl * as.numeric(sign(model_peak$tau))
      
      
      ## Syndromes - part B of Project
      
      # arable land    
      sig <- 0.1 # pseudo p-value
      if (incr(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1) 
      {syndrome[i,j] <- 3} 
      # +oo = increasing biomass
      #...
      # semi-natural      
      
      if (abs(pintegral) < sig & pintegral > 0 & abs(ppeak) < sig & ppeak > 0 & abs(pamp) > sig &  luc[i,j] == 2)
      {syndrome[i,j] <- 1} 
      # ++o =  decreasing biomass
      # ...        
      
    }
  }
}
close(pb)

model_integral = MannKendall(ts(magnitude[10,10,])) # inspecting
model_integral                             
                             
x=raster(per_annual)