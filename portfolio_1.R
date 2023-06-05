library(raster)
library(caTools)
library(Kendall)

filename = 'GIMMS_Iberian'
setwd('U:\\Time Series Analysis\\s4')
cube = read.ENVI(filename,headerfile=paste(filename,'.hdr',sep='')) #from caTools

dim(cube)
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


dim = c(years,ncol(cube),nrow(cube))
t = array(NA,dim)
magnitude = aperm(t)
peaktime = magnitude
meanNDVI = magnitude
integral = array(NA, c(nrows, ncols))

pb = txtProgressBar(title='progress bar',min=0,max =nrows) #run together with loop
for (i in 1:nrows) {
  for (j in 1:ncols) {
    setTxtProgressBar(pb, i, title=round(i/nrow(cube)*100,0))
    GIMMS = cube[i,j,]
    
    if (var(GIMMS) > 0) {
      time=1:12
      for (k in 1:years) {
        # Fit Fourier Polynomials to monthly NDVI data derive the following phenological
        # properties in each year: a) NDVI peaking times, b) NDVI magnitudes, c) NDVI integrals.
        model = lm(GIMMS[time] ~ fp)
        magnitude[i,j,k] = max(model$fit) - min(model$fit)
        peaktime[i,j,k] = which(model$fit == max(model$fit))
        meanNDVI[i,j,k] = mean(GIMMS[time])
        integral[k] = sum(GIMMS[time])
        # Fit a trend model (Seasonal Kendall test) to derive monotonic trends in the three phenological time-series.
        
        time = time + 12
      }
      model_integral = MannKendall(ts(integral)) 
      model_amp = MannKendall(ts(magnitude))
      model_peak = MannKendall(ts(peaktime))
      
      pintegral = model_integral$sl * as.numeric(sign(model_integral))
      pamp = model_amp$sl * as.numeric()
      ppeak = model_peak # finish later see image
      
    }
  }
  
}
close(pb)

model_integral = MannKendall(ts(magnitude[10,10,])) # inspecting
model_integral                             
                             
x=raster(per_annual)