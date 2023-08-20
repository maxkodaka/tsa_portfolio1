library(raster)
library(caTools)
library(Kendall)

filename = 'GIMMS_Iberian'
filename_class = 'clc8km_clip_recode_Iberian_2classes'
file_out = '/home/maxim/Documents/coursework/time-series-analysis/s6/syndrome_raster.tif'
setwd('/home/maxim/Documents/coursework/time-series-analysis/data/')
#setwd('U:\\Time Series Analysis\\s5')
cube = read.ENVI(filename,headerfile=paste(filename,'.hdr',sep='')) #from caTools
luc<-read.ENVI(filename_class,headerfile=paste(filename_class,'.hdr',sep='')) #

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

#syndrome<-matrix(data <- -999, nrow <- dim(cube)[1],ncol <- dim(cube)[2])
syndrome<-matrix(data <- NA, nrow <- dim(cube)[1],ncol <- dim(cube)[2])
numsigtrends005 = syndrome
numsigtrends010 = syndrome
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
      
      model_integral = MannKendall(ts(integral))     
      model_amp = MannKendall(ts(magnitude))         
      model_peak = MannKendall(ts(peaktime))         
      
      # Content of the MannKendall object: tau, sl, S, D, varS
      
      # There is some controversy about the following. It was unclear why it's multiplied by the sign.
      # I guess because when we calculate the syndromes, it matters whether it is increasing or decreasing.
      # But then why not just use the sign of tau in the boolean conditional statements?
      
      pintegral = model_integral$sl * as.numeric(sign(model_integral$tau))
      pamp = model_amp$sl * as.numeric(sign(model_amp$tau))
      ppeak = model_peak$sl * as.numeric(sign(model_peak$tau))
      
      # Part A.3: Identify those pixels with at least two significant (α=5% and 10%) trends concerning
      # NDVI, peaking time or magnitude.
      # We count how many significant trends there are for each pixel:
      # For significance of 0.05
      
      # Here I use the model_integral, model_amp and model_peak, the quantities that were not multiplied by the sign of tau.
      # For the conditionals as stated here the range needs to be between zero and one.
      s = 0
      
      if(!is.na(model_amp$sl) & model_amp$sl < 0.05){
        s = s + 1 
      }
      
      if(!is.na(model_peak$sl) & model_peak$sl < 0.05){
        s = s + 1 
      } 
      
      if(!is.na(model_integral$sl) & model_integral$sl < 0.05){
        s = s + 1 
      }
      
      numsigtrends005[i,j] = s
      
      # For significance of 0.10
      
      s = 0
      
      if(!is.na(model_amp$sl) & model_amp$sl < 0.10){
        s = s + 1 
      }
      
      if(!is.na(model_peak$sl) & model_peak$sl < 0.10){
        s = s + 1 
      } 
      
      if(!is.na(model_integral$sl) & model_integral$sl < 0.10){
        s = s + 1 
      }
      
      numsigtrends010[i,j] = s      

      ## Syndromes - part B of Project
      
      ## increasing biomass - 1
      ## decreasing biomass - 2
      ## changing crop - 3
      ## increasing productivity - 4
      ## decreasing productivity - 5
      ## increasing irrigation - 6
      ## decreasing irrigation - 7
      
      ## luc values: Is this correct??
      ## arable land: 1
      ## semi-natural: 2
      
      sig <- 0.1 # pseudo p-value
      
      
      # semi-natural    
      
      # increasing biomass = 1
      
      # +oo
      if (incr(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 2) 
      {syndrome[i,j] <- 1}
      
      # +o+
      if (incr(pintegral,sig) & nought(ppeak,sig) & incr(pamp,sig) &  luc[i,j] == 2)
      {syndrome[i,j] <- 1} 
      
      # ++o
      if (incr(pintegral,sig) & incr(ppeak,sig) & nought(pamp,sig) &  luc[i,j] == 2)
      {syndrome[i,j] <- 1} 

      # +-o
      if (incr(pintegral,sig) & decr(ppeak,sig) & nought(pamp,sig) &  luc[i,j] == 2)
      {syndrome[i,j] <- 1} 
      
      # +o-
      if (incr(pintegral,sig) & nought(ppeak,sig) & decr(pamp,sig) &  luc[i,j] == 2)
      {syndrome[i,j] <- 1} 
      
      # decreasing biomass = 2
      
      # ---
      if (decr(pintegral,sig) & decr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 2)
      {syndrome[i,j] <- 2}
      
      # -oo
      if (decr(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 2)
      {syndrome[i,j] <- 2}
      
      # -o-
      if (decr(pintegral,sig) & nought(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 2)
      {syndrome[i,j] <- 2}
      
      # --o
      if (decr(pintegral,sig) & decr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 2)
      {syndrome[i,j] <- 2}
      
      # arable land (luc = 1)
      
      # increasing productivity = 4
      
      # +oo
      if (incr(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 4}
      
      # +o+
      if (incr(pintegral,sig) & nought(ppeak,sig) & incr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 4}
      
      # decreasing productivity = 5
      
      # -oo
      if (decr(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 5}
      
      
      # -o-
      if (decr(pintegral,sig) & nought(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 5}
      
      
      # expansion of irrigated arable land = 6
      
      # +++
      if (incr(pintegral,sig) & incr(ppeak,sig) & incr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 6}
      
      # ++o
      if (incr(pintegral,sig) & incr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 6}
      
      # ++-
      if (incr(pintegral,sig) & incr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 6}
      
      
      # decline of irrigated arable land = 7
      
      # ---
      if (decr(pintegral,sig) & decr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 7}
      
      
      # --o
      if (decr(pintegral,sig) & decr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 7}
      
      # o-o
      if (nought(pintegral,sig) & decr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 7}
      
      # o--
      if (nought(pintegral,sig) & decr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 7}
      
      
      # crop change = 3
      
      # oo-
      if (nought(pintegral,sig) & nought(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 3}
      
      # oo+
      if (nought(pintegral,sig) & nought(ppeak,sig) & incr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 3}
      
      # o+-
      if (nought(pintegral,sig) & incr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 3}
      
      # o+o
      if (nought(pintegral,sig) & incr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 3}
      
      # o++
      if (nought(pintegral,sig) & incr(ppeak,sig) & incr(pamp,sig) & luc[i,j] == 1)
      {syndrome[i,j] <- 3}
      
      # no change
      
      # ooo
      if (nought(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig))
      {syndrome[i,j] <- 0}

      

      
    }
  }
}
close(pb)

cube = brick(filename)
syndrome_raster=raster(syndrome,xmn=extent(cube)[1], xmx=extent(cube)[2], ymn=extent(cube)[3], ymx=extent(cube)[4], crs=projection(cube))
writeRaster(syndrome_raster, file_out, format="GTiff", overwrite=TRUE)
KML (syndrome_raster, file ='syndrome',  col = hcl.colors(8, "YlOrRd", rev = TRUE) ,
     maxpixels = ncell ( raster (pval_raster)),overwrite="TRUE")