##########################
######## IMPORTS ########
########################

library(raster)
library(caTools)
library(Kendall)
library(sf)
library(geojsonio)
library(flextable)

########################
######## PATHS ########
######################

setwd('/home/maxim/Documents/coursework/time-series-analysis/data/')


#### input paths ####
####################

# path to current files in LAEA projection

filename = 'GIMMS_Iberian'
filename_class = 'clc8km_clip_recode_Iberian_2classes'

# original version in WGS84 projection

#filename = '/home/maxim/Documents/coursework/time-series-analysis/data/test/GIMMS_CLC_lamb_conf_NEW/GIMMS_Iberian'
#filename_class = '/home/maxim/Documents/coursework/time-series-analysis/data/test/GIMMS_CLC_lamb_conf_NEW/GIMMS_CLC_lamb_conf_NEW/clc8km_clip_recode_Iberian_2classes'

path_portugal = 'gadm_PRT_0.json'
path_spain = 'gadm_ESP_0.json'


#### output paths ####
#####################

# original output data path
#outdir = '../output'

# cropped data path
outdir = '../output/seasonalmk/'

## syndrome maps##
#################

file_out = paste(outdir,'syndrome_raster.tif')
file_out_syndrome_natural = paste(outdir,'syndrome_natural.tif')
file_out_syndrome_arable = paste(outdir,'syndrome_arable.tif')

## significant trend maps ##
###########################

file_out_numsigtrends005 = paste(outdir,'numsigtrends_005.tif')
file_out_numsigtrends010 = paste(outdir,'numsigtrends_010.tif')

file_out_twosigtrends005 = paste(outdir,'twosigtrends_005.tif')
file_out_twosigtrends010 = paste(outdir,'twosigtrends_010.tif')

## syndrome stats ##
###################

file_out_syndrome_stats_natural = paste(outdir,'syndrome_stats_natural.png')
file_out_syndrome_stats_arable = paste(outdir,'syndrome_stats_arable.png')

###############################
#### Input and conversion ####
#############################

# Case 1: read it in just as an array
cube = read.ENVI(filename,headerfile=paste(filename,'.hdr',sep='')) #from caTools

# Case 2: read it in as a data cube
# The array format in Case 1 is not readable by the crop function. But can it be iterated in the same way?
# -Yes it can with some minor tweaks.
cube_brick = brick(filename)


# import spain and portugal national boundaries geojson data
geom_prt = st_read(path_portugal)
geom_sp = st_read(path_spain)

# combine to create iberia geojson
geom_iberia = st_union(geom_prt$geometry,geom_sp$geometry)

# Convert to a spatial polygon object
geom_iberia <- st_as_sf(geom_iberia)

geom_iberia = st_transform(geom_iberia, crs = crs(cube_brick))

#crop cube to Iberian peninsula
cube_brick = crop(cube_brick, geom_iberia)


# convert back to array as in Prof. Udelhoven's original implementation.
cube2 = as.array(cube_brick)

# binary classification map
# Case 1: Read in simply as an array
luc<-read.ENVI(filename_class,headerfile=paste(filename_class,'.hdr',sep='')) #

# Case 2: Read in as a raster so that we have the crs for cropping purposes.
luc_rast = raster(filename_class)

# crop to Iberian peninsula
luc_rast = crop(luc_rast, geom_iberia)

# Convert back to a matrix as in Prof. Udelhoven's original implementation.
luc = as.matrix(luc_rast)


##############################
#### Fourier Polynomials ####
############################

lin =1:12
x1 = sin(2*pi*1/12*(lin))
x2 = cos(2*pi*1/12*(lin))
x3 = sin(2*pi*2/12*(lin))
x4 = cos(2*pi*2/12*(lin))
fp = cbind(lin,x1,x2,x3,x4)



#########################
#### INITIALIZATION ####
#######################

years = dim(cube)[3]/12

nrows=nrow(cube)
ncols=ncol(cube)

dim = c(years,ncol(cube),nrow(cube))

meanNDVI = aperm(array(NA,dim))                     # dimensions: 93, 152, 32

magnitude = array(NA, 32)                         
peaktime = array(NA, 32)                          
integral = array(NA, 32)                         

pintegral = array(NA, 32) 
pamp = array(NA, 32)      
ppeak = array(NA, 32)   

syndrome<-matrix(data <- NA, nrow <- dim(cube)[1],ncol <- dim(cube)[2])
syndrome_natural<-matrix(data <- NA, nrow <- dim(cube)[1],ncol <- dim(cube)[2])
syndrome_arable<-matrix(data <- NA, nrow <- dim(cube)[1],ncol <- dim(cube)[2])
numsigtrends005 = syndrome
numsigtrends010 = syndrome
twosigtrends005 = syndrome
twosigtrends010 = syndrome

## significant trend stats tables for p = 0.05 and 0.10

sigtrendstats005 = data.frame(
  s = c(0,1,2,3),
  arable_cnt = array(0, 4),
  seminatural_cnt = array(0, 4)
)

sigtrendstats010 = data.frame(
  s = c(0,1,2,3),
  arable_cnt = array(0, 4),
  seminatural_cnt = array(0, 4)
)


##########################
#### TREND FUNCTIONS ####
########################


## increasing
incr <- function(x, sig){
  if(!is.na(x)){
    abs(x) < sig & x > 0
  }
  else{
    FALSE
  }
}

## decreasing
decr <- function(x, sig){
  if(!is.na(x)){
    abs(x) < sig & x < 0
  }
  else{
    FALSE
  }
}

## no change
nought <- function(x, sig){
  if(is.na(x)){
    TRUE
  }
  else{
    abs(x) > sig
  }
}

####################
#### MAIN LOOP ####
##################

pb = txtProgressBar(title='progress bar',min=0,max =nrows)
for (i in 1:nrows) {
  for (j in 1:ncols) {
    setTxtProgressBar(pb, i, title = round(i / nrow(cube) * 100, 0))
    GIMMS = cube[i, j,] 
    
    if (var(GIMMS) > 0 & !is.na(GIMMS[1])) { #added !is.na(GIMMS) because I am using a different type of datacube now
      time = 1:12
      for (k in 1:years) {
        
        ##################################
        #### FIT FOURIER POLYNOMIALS ####
        ################################
        
        # Fit Fourier Polynomials to monthly NDVI data to derive the following phenological
        # properties in each year: a) NDVI peaking times, b) NDVI magnitudes, c) NDVI integrals.
        
        model = lm(GIMMS[time] ~ fp)
        
        magnitude[k] = max(model$fit) - min(model$fit)      
        peaktime[k] = which(model$fit == max(model$fit))  
        integral[k] = sum(GIMMS[time])                      
        
        meanNDVI[i, j, k] = mean(GIMMS[time])
        
        time = time + 12
      }
      
      #######################
      #### MANN KENDALL ####
      #####################
      
      # Fit a trend model (Seasonal Kendall test) to derive monotonic trends in the three phenological time-series.
      
      
      model_integral = SeasonalMannKendall(ts(integral))     
      model_amp = SeasonalMannKendall(ts(magnitude))         
      model_peak = SeasonalMannKendall(ts(peaktime))         
      
      # Content of the MannKendall object: tau, sl, S, D, varS
      
      ##########################
      #### PSEUDO P-VALUES ####
      ########################
      
      pintegral = model_integral$sl * as.numeric(sign(model_integral$tau))
      pamp = model_amp$sl * as.numeric(sign(model_amp$tau))
      ppeak = model_peak$sl * as.numeric(sign(model_peak$tau))
      
      ##################################
      #### SIGNIFICANT TREND STATS ####
      ################################
      
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
      
      if(s>=2){
        twosigtrends005[i,j] = 1
      } else(
        twosigtrends005[i,j] = 0
      )
      
      # Counts of trends by land cover class for significance of 0.05
      if (!is.na(luc[i,j])){
        sigtrendstats005[s+1,luc[i,j]+1] = sigtrendstats005[s+1,luc[i,j]+1] + 1
      }
      
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
      
      if(s>=2){
        twosigtrends010[i,j] = 1
      } else(
        twosigtrends010[i,j] = 0
      )
      
      # Counts of trends by land cover class for significance of 0.10
      if (!is.na(luc[i,j])){
        sigtrendstats010[s+1,luc[i,j]+1] = sigtrendstats010[s+1,luc[i,j]+1] + 1
      }

      ####################
      #### SYNDROMES ####
      ##################
      
      #part B of Project
      
      ## Syndrome definitions
      
      ## increasing biomass - 1
      ## decreasing biomass - 2
      ## changing crop - 3
      ## increasing productivity - 4
      ## decreasing productivity - 5
      ## increasing irrigation - 6
      ## decreasing irrigation - 7
      
      
      ## land cover classes
      
      ## arable land: 1
      ## semi-natural: 2
      
      sig <- 0.1 # pseudo p-value
      
      
      # semi-natural    
      
      # increasing biomass = 1
      if(!is.na(luc[i,j]))
      {
        # +oo
        if (incr(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 2) 
        {syndrome_natural[i,j] <- 1}
        
        # +o+
        if (incr(pintegral,sig) & nought(ppeak,sig) & incr(pamp,sig) &  luc[i,j] == 2)
        {syndrome_natural[i,j] <- 1} 
        
        # ++o
        if (incr(pintegral,sig) & incr(ppeak,sig) & nought(pamp,sig) &  luc[i,j] == 2)
        {syndrome_natural[i,j] <- 1} 
        
        # +-o
        if (incr(pintegral,sig) & decr(ppeak,sig) & nought(pamp,sig) &  luc[i,j] == 2)
        {syndrome_natural[i,j] <- 1} 
        
        # +o-
        if (incr(pintegral,sig) & nought(ppeak,sig) & decr(pamp,sig) &  luc[i,j] == 2)
        {syndrome_natural[i,j] <- 1} 
        
        # decreasing biomass = 2
        
        # ---
        if (decr(pintegral,sig) & decr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 2)
        {syndrome_natural[i,j] <- 2}
        
        # -oo
        if (decr(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 2)
        {syndrome_natural[i,j] <- 2}
        
        # -o-
        if (decr(pintegral,sig) & nought(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 2)
        {syndrome_natural[i,j] <- 2}
        
        # --o
        if (decr(pintegral,sig) & decr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 2)
        {syndrome_natural[i,j] <- 2}
        
        # arable land (luc = 1)
        
        # increasing productivity = 4
        
        # +oo
        if (incr(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 4}
        
        # +o+
        if (incr(pintegral,sig) & nought(ppeak,sig) & incr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 4}
        
        # decreasing productivity = 5
        
        # -oo
        if (decr(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 5}
        
        
        # -o-
        if (decr(pintegral,sig) & nought(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 5}
        
        
        # expansion of irrigated arable land = 6
        
        # +++
        if (incr(pintegral,sig) & incr(ppeak,sig) & incr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 6}
        
        # ++o
        if (incr(pintegral,sig) & incr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 6}
        
        # ++-
        if (incr(pintegral,sig) & incr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 6}
        
        
        # decline of irrigated arable land = 7
        
        # ---
        if (decr(pintegral,sig) & decr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 7}
        
        
        # --o
        if (decr(pintegral,sig) & decr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 7}
        
        # o-o
        if (nought(pintegral,sig) & decr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 7}
        
        # o--
        if (nought(pintegral,sig) & decr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 7}
        
        
        # crop change = 3
        
        # oo-
        if (nought(pintegral,sig) & nought(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 3}
        
        # oo+
        if (nought(pintegral,sig) & nought(ppeak,sig) & incr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 3}
        
        # o+-
        if (nought(pintegral,sig) & incr(ppeak,sig) & decr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 3}
        
        # o+o
        if (nought(pintegral,sig) & incr(ppeak,sig) & nought(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 3}
        
        # o++
        if (nought(pintegral,sig) & incr(ppeak,sig) & incr(pamp,sig) & luc[i,j] == 1)
        {syndrome_arable[i,j] <- 3}
        
        # no change
        
        # ooo
        if (nought(pintegral,sig) & nought(ppeak,sig) & nought(pamp,sig))
        {
          syndrome_arable[i,j] <- 0
          syndrome_natural[i,j] <- 0
        
        }
        
      }

      

      
    }
  }
}
close(pb)

###################################################
#### SIGNIFICANT TREND CUMULATIVE PERCENTAGES ####
#################################################

for(i in 4:1){
  sigtrendstats005[i,'arable_cumpercent'] = sum(sigtrendstats005['arable_cnt'][sigtrendstats005['s']>=i-1])/sum(sigtrendstats005['arable_cnt'])
  sigtrendstats010[i,'arable_cumpercent'] = sum(sigtrendstats010['arable_cnt'][sigtrendstats010['s']>=i-1])/sum(sigtrendstats010['arable_cnt'])
  
  sigtrendstats005[i,'seminatural_cumpercent'] = sum(sigtrendstats005['seminatural_cnt'][sigtrendstats005['s']>=i-1])/sum(sigtrendstats005['seminatural_cnt'])
  sigtrendstats010[i,'seminatural_cumpercent'] = sum(sigtrendstats010['seminatural_cnt'][sigtrendstats010['s']>=i-1])/sum(sigtrendstats010['seminatural_cnt'])
}


##############################
#### SYNDROME STATISTICS ####
############################

syndrome_stats = data.frame(
  syndrome_index = c(1,2,3,4,5,6,7,8),
  syndrome = c(
    'increasing_biomass',
    'decreasing_biomass',
    'crop change',
    'increasing productivity',
    'decreasing productivity',
    'expansion of irrigated arable land',
    'decline of irrigated arable land',
    'no change'),
  pixel_count = array(0,8)
)

syndrome_stats_natural = data.frame(
  syndrome = c(
    'increasing_biomass',
    'decreasing_biomass',
    'no change'),
  pixel_count = array(0,3)
)

syndrome_stats_arable = data.frame(
  syndrome = c(
    'crop change',
    'increasing productivity',
    'decreasing productivity',
    'expansion of irrigated arable land',
    'decline of irrigated arable land',
    'no change'),
  pixel_count = array(0,6)
)

for(i in 1:8){
  if(i<=2){
    syndrome_stats[i,'pixel_count'] = length(which(syndrome_natural==i))
  }
  else if(i<=7){
    syndrome_stats[i,'pixel_count'] = length(which(syndrome_arable==i))
  }
  else{
    syndrome_stats[i,'pixel_count'] = dim(luc)[1]*dim(luc)[2] - sum(syndrome_stats[1:7,'pixel_count'])
  }

}

for(i in 1:3){
  if(i<=2){
    syndrome_stats_natural[i,'pixel_count'] = length(which(syndrome_natural==i))
  }
  else{
    syndrome_stats_natural[i,'pixel_count'] = length(which(luc==2)) - sum(syndrome_stats_natural[1:2,'pixel_count'])
  }
  
}
syndrome_stats_natural['percent'] = 100*syndrome_stats_natural['pixel_count']/sum(syndrome_stats_natural['pixel_count'])

for(i in 3:8){
  if(i<=7){
    syndrome_stats_arable[i-2,'pixel_count'] = length(which(syndrome_arable==i))
  }
  else{
    syndrome_stats_arable[i-2,'pixel_count'] = length(which(luc==1)) - sum(syndrome_stats_arable[1:5,'pixel_count'])
  }
  
}
syndrome_stats_arable['percent'] = 100*syndrome_stats_arable['pixel_count']/sum(syndrome_stats_arable['pixel_count'])

#########################
#### WRITING OUTPUT ####
#######################

## Raster to Geotiff ##
######################


# Syndromes #
############
syndrome_raster=raster(syndrome,xmn=extent(cube_brick)[1], xmx=extent(cube_brick)[2], ymn=extent(cube_brick)[3], ymx=extent(cube_brick)[4], crs=projection(cube_brick))
writeRaster(syndrome_raster, file_out, format="GTiff", overwrite=TRUE)

syndrome_natural_raster=raster(syndrome_natural,xmn=extent(cube_brick)[1], xmx=extent(cube_brick)[2], ymn=extent(cube_brick)[3], ymx=extent(cube_brick)[4], crs=projection(cube_brick))
writeRaster(syndrome_natural_raster, file_out_syndrome_natural, format="GTiff", overwrite=TRUE)

syndrome_arable_raster=raster(syndrome_arable,xmn=extent(cube_brick)[1], xmx=extent(cube_brick)[2], ymn=extent(cube_brick)[3], ymx=extent(cube_brick)[4], crs=projection(cube_brick))
writeRaster(syndrome_arable_raster, file_out_syndrome_arable, format="GTiff", overwrite=TRUE)

# significant trends #
#####################

numsigtrends005_raster=raster(numsigtrends005,xmn=extent(cube_brick)[1], xmx=extent(cube_brick)[2], ymn=extent(cube_brick)[3], ymx=extent(cube_brick)[4], crs=projection(cube_brick))
numsigtrends010_raster=raster(numsigtrends010,xmn=extent(cube_brick)[1], xmx=extent(cube_brick)[2], ymn=extent(cube_brick)[3], ymx=extent(cube_brick)[4], crs=projection(cube_brick))

twosigtrends005_raster=raster(twosigtrends005,xmn=extent(cube_brick)[1], xmx=extent(cube_brick)[2], ymn=extent(cube_brick)[3], ymx=extent(cube_brick)[4], crs=projection(cube_brick))
twosigtrends010_raster=raster(twosigtrends010,xmn=extent(cube_brick)[1], xmx=extent(cube_brick)[2], ymn=extent(cube_brick)[3], ymx=extent(cube_brick)[4], crs=projection(cube_brick))

writeRaster(numsigtrends005_raster, file_out_numsigtrends005, format="GTiff", overwrite=TRUE)
writeRaster(numsigtrends010_raster, file_out_numsigtrends010, format="GTiff", overwrite=TRUE)

writeRaster(twosigtrends005_raster, file_out_twosigtrends005, format="GTiff", overwrite=TRUE)
writeRaster(twosigtrends010_raster, file_out_twosigtrends010, format="GTiff", overwrite=TRUE)



## Convert dataframes to flextables and output to docx or png ##
###############################################################

set_flextable_defaults(
  font.size = 10, #theme_fun = theme_vanilla,
  padding = 6,
  )

## saving as docx (results don't look good)
#save_as_docx(set_caption(flextable(syndrome_stats_arable),path=file_out_syndrome_stats_arable),'Arable Land Change Syndromes')
#save_as_docx(set_caption(flextable(syndrome_stats_natural),path=file_out_syndrome_stats_natural),'Semi-Natural Land Change Syndromes')

ft_a = flextable(syndrome_stats_arable)
ft_n = flextable(syndrome_stats_natural)

ft_a = set_caption(ft_a,'Arable Land Change Syndromes')
ft_n = set_caption(ft_n,'Semi-Natural Land Change Syndromes')

save_as_image(ft_a,path=paste(outdir,'syndrome_stats_arable.png'))
save_as_image(ft_n,path=paste(outdir,'syndrome_stats_natural.png'))

## Iberia Geometry ##
####################
geojson_write(st_sf(geometry=geom_iberia),file=paste(outdir,'gadm_IBR_0.json'))
