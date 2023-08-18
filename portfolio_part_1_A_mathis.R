library(raster)
library(Kendall) # for Mann-Kendall test
library(caTools) # to import/export ENVI files
library(bfast)   # for break detection

# Step 1: Import raster data cube
filename = "GIMMS_Iberian"
cube = read.ENVI(filename, headerfile = paste(filename, ".hdr", sep = ""))

# create fourier polynomials
years = dim(cube)[3]/12 # cube[3] --> time in months (divide by 12 for year)

x1 = sin(2*pi*1/12*(1:12))
x2 = cos(2*pi*1/12*(1:12))
x3 = sin(2*pi*2/12*(1:12))
x4 = cos(2*pi*2/12*(1:12))
lin = 1:12

fp = cbind(lin, x1, x2, x3, x4) # fourier polynomial
#cor(fp)

nrows = nrow(cube)
ncols = ncol(cube)

# define empty arrays for parameters
dim = c(years, ncols, nrows) # 3D array
t = array(NA, dim)           # empty array as template

magnitude = aperm(t) # aperm = transpose
peaktime = magnitude
meanNDVI = magnitude

dim = c(ncols, nrows) # 2D array
t = array(NA, dim)

res_mag = aperm(t)
res_peak = res_mag
res_mean = res_mag

tau_mag = aperm(t)
tau_peak = tau_mag
tau_mean = tau_mag

pval_mag = aperm(t)
pval_peak = pval_mag
pval_mean = pval_mag

sig = pval_mag

# Initialize progress bar
pb = winProgressBar(title = "Progress bar", min = 0, max = nrows, width = 300)

start.time = Sys.time()

# Loop through each pixel and calculate test statistics
for (i in 1:nrows) {
  
  setWinProgressBar(pb, i, title = paste(round(i/nrows*100, 0),"% done"))
  
  for (j in 1:ncols) {
    
    # Extract time-series for current pixel
    GIMMS = cube[i, j, ]
    #GIMMS_ts = ts(GIMMS)
    
    # Calculate Seasonal Mann-Kendall test statistics
    if (var(GIMMS) > 0) {
      
      time = 1:12
      
      for (k in 1:years){
        #modelling the time series data to calculate magnitude, peaktime and meanNdvi
        model = lm(GIMMS[time] ~ fp)
        magnitude[i,j,k] = max(model$fit) - min(model$fit)
        peaktime[i,j,k] = which(model$fit == max(model$fit))
        meanNDVI[i,j,k] = mean(GIMMS[time])
        time = time + 12
      }
      
      # Mann Kendall Test: to dêtrmine if there is any changes of these 3 values(dêtcting trends)
      res_mag = SeasonalMannKendall(ts(magnitude[i,j,]))
      res_peak = SeasonalMannKendall(ts(peaktime[i,j,]))
      res_mean = SeasonalMannKendall(ts(meanNDVI[i,j,]))
      
      # write magnitude, peaktime and mean NDVi values in arrays
      tau_mag[i,j] = res_mag$tau
      tau_peak[i,j] = res_peak$tau
      tau_mean[i,j] = res_mean$tau
      
      # do the same for p-values 
      # MK - This part doesn't make sense!! Because in the next part everything with negative trend will be counted as siginifant!
      # W
      pval_mag[i,j] = res_mag$sl * as.numeric(sign(tau_mag[i,j]))
      pval_peak[i,j] = res_peak$sl * as.numeric(sign(tau_peak[i,j]))
      pval_mean[i,j] = res_mean$sl * as.numeric(sign(tau_mean[i,j]))
      
      s = 0
      
      if(!is.na(pval_mag[i,j]) & pval_mag[i,j] < 0.05){
       s = s + 1 
      }
      
      if(!is.na(pval_peak[i,j]) & pval_peak[i,j] < 0.05){
       s = s + 1 
      } 
      
      if(!is.na(pval_mean[i,j]) & pval_mean[i,j] < 0.05){
       s = s + 1 
      }
      
      sig[i,j] = s
    }  
  }
}

close(pb)

end.time = Sys.time()
time.taken = end.time - start.time
time.taken

# Step 3: Export results as tif files
# Create output raster files

cube = brick(filename)

tau_mag_raster = raster(tau_mag, xmn=extent(cube)[1], xmx=extent(cube)[2], ymn=extent(cube)[3], ymx=extent(cube)[4], crs=projection(cube))
tau_peak_raster = raster(tau_peak, xmn=extent(cube)[1], xmx=extent(cube)[2], ymn=extent(cube)[3], ymx=extent(cube)[4], crs=projection(cube))
tau_mean_raster = raster(tau_mean, xmn=extent(cube)[1], xmx=extent(cube)[2], ymn=extent(cube)[3], ymx=extent(cube)[4], crs=projection(cube))

pval_mag_raster = raster(pval_mag, xmn=extent(cube)[1], xmx=extent(cube)[2], ymn=extent(cube)[3], ymx=extent(cube)[4], crs=projection(cube))
pval_peak_raster = raster(pval_peak, xmn=extent(cube)[1], xmx=extent(cube)[2], ymn=extent(cube)[3], ymx=extent(cube)[4], crs=projection(cube))
pval_mean_raster = raster(pval_mean, xmn=extent(cube)[1], xmx=extent(cube)[2], ymn=extent(cube)[3], ymx=extent(cube)[4], crs=projection(cube))

sig_raster = raster(sig, xmn=extent(cube)[1], xmx=extent(cube)[2], ymn=extent(cube)[3], ymx=extent(cube)[4], crs=projection(cube))


# Set names for output files
tau_mag_file = "tau_mag.tif"
tau_peak_file = "tau_peak.tif"
tau_mean_file = "tau_mean.tif"

pval_mag_file = "pval_mag.tif"
pval_peak_file = "pval_peak.tif"
pval_mean_file = "pval_mean.tif"

sig_file = "sig.tif"


# Write output files
writeRaster(tau_mag_raster, tau_mag_file, format="GTiff", overwrite=TRUE)
writeRaster(tau_peak_raster, tau_peak_file, format="GTiff", overwrite=TRUE)
writeRaster(tau_mean_raster, tau_mean_file, format="GTiff", overwrite=TRUE)

writeRaster(pval_mag_raster, pval_mag_file, format="GTiff", overwrite=TRUE)
writeRaster(pval_peak_raster, pval_peak_file, format="GTiff", overwrite=TRUE)
writeRaster(pval_mean_raster, pval_mean_file, format="GTiff", overwrite=TRUE)

writeRaster(sig_raster, sig_file, format="GTiff", overwrite=TRUE)
