#===================================================================================================================#
#=====================================    Spatial analysis & plotting functions     ================================#


# Function to process shape file and clustered dataframe to plot semi-variogram # 

process_plot_variogram_func <- function(shp, clustered_data, variable_interest) {
  
  grid <- makegrid(shp, cellsize = 0.0083) # sample  point locations within spatial object (shp) with given cellsize
  
  coordinates(grid) <- ~x1 + x2 # make spatial object with dataframe of sample points in grid 
  
  geodata <- clustered_data
  
  coordinates(geodata) <- ~x + y # make spatial object with lat and long from dataset with livestock data
  
  if(variable_interest == "pg_h"){
    variable_interest = geodata$Pg_h
  }
  
  if(variable_interest == "pg_2"){
    variable_interest = geodata$Pg_2
  }
  
  vgm = variogram(variable_interest ~ 1, geodata) # make a variogram (plotting spatial autocorrelation: https://stats.idre.ucla.edu/r/faq/how-do-i-generate-a-variogram-for-spatial-data-in-r/) of housholds with pigs
  
  fit = fit.variogram(vgm, model = vgm(0.03, "Sph", 2, 0.01)) # model fit variogram fit (?: https://r-spatial.org/r/2016/02/14/gstat-variogram-fitting.html)
  
  vmg <- plot(vgm, fit) # plot a semi-variogram (spatial autocorr) data vs. model fit (line)
  
  return(list(geodata, grid, fit, vmg))
  
}


# Function: kriging number of variable of interest 
# e.g. (kriging = estimating avg. variable of interest
# at locations where dont have measurements : https://rpubs.com/nabilabd/118172)

kriging_func <- function(spatial_object, admin, variable_interest) {
  geodata = spatial_object[[1]] 
  grid = spatial_object[[2]] 
  fit = spatial_object[[3]]
  
  if(variable_interest == "pg_h"){
    variable_interest = geodata$Pg_h
  }
  
  if(variable_interest == "pg_2"){
    variable_interest = geodata$Pg_2
  }
  
  kriged = krige(variable_interest ~ 1, geodata, grid, model = fit) # kriging = estimating avg. variable of interest 
  # at locations where dont have measurements : https://rpubs.com/nabilabd/118172)
  # n.b. this takes some time to run
  
  # proportion of households with pigs (predicted with kriging)
  kriged_plot <-
    ggplot(data = kriged@data, aes(x = kriged$x1, y = kriged$x2)) +
    geom_tile(aes(fill = kriged$var1.pred)) +
    coord_equal() +
    scale_fill_gradient(low = "yellow", high = "red") +
    geom_polygon(data = admin, aes(x = long, y = lat, group = group), color = "black", alpha = 0) +
    theme_bw()
  
  return(list(kriged, kriged_plot))
  
}

# specify spatial data and create raster as being gridded function  
spatially_process_variable_function <- function(kriged){
  
  gridded(kriged) <- T 
  
  variable_spatially_processed <- raster(kriged , "var1.pred")
  
  return(variable_spatially_processed)
}


# This function assesses flock sizes from households who have reported owning pigs : using IDW for spatial interpolation

pig_flock_size_HH_interpolated_func <- function(admin_processed, clustered_data, spatial_object, variable_interest) {
  
  geodata = spatial_object[[1]] 
  grid = spatial_object[[2]] 
  
  if(variable_interest == "pg_hs"){
    variable_interest = geodata$Pg_hs
  }
  
  pig_flock_plot1 <- 
    ggplot() + 
    geom_polygon(data = admin_processed, aes(x = long, y = lat, group = group), fill = "grey80", colour = "grey90", alpha = 1) + 
    geom_point(data = clustered_data, aes(x = x, y = y,color = Pg_hs) ) + 
    scale_colour_gradientn(colours = heat.colors(12))+    
    coord_equal(ratio = 1)
  
  idw = idw(Pg_hs ~ 1, geodata[!is.na(geodata$Pg_hs),], grid, maxdist = Inf, idp = 2) # interpolation of flock/herd sizes using [inverse distance weighted interpolation]
  
  pig_flock_plot2 <- 
    ggplot(data = idw@data, aes(x = idw$x1, y = idw$x2)) + 
    geom_tile(aes(fill = idw$var1.pred)) + 
    coord_equal() +
    scale_fill_gradient(low = "yellow", high = "red") +
    geom_polygon(data = admin_processed, aes(x = long, y = lat, group = group), color = "black", alpha = 0) + 
    theme_bw()
  
  return(list(idw, pig_flock_plot1, pig_flock_plot2))
  
}


# Function to compute pig density - based on the following calculation : 
# 2011: human population (from world pop 2015 unadjusted) x percent of household keeping pigs x flock size
# 2006: human population (from world pop 2006 unadjusted) x percent of household keeping pigs x flock size
# word pop unadjusted link: https://www.worldpop.org/geodata/listing?id=76 

compute_pig_denisty_func <- function(pop, kriged, shp, idw) {
  
  gridded(kriged) <- T 
  
  Pg_h_r <- raster(kriged , "var1.pred") # create raster layer object with HH with pigs data
  
  Pg_h_r <- mask(Pg_h_r , shp) # mask NA values in Raster object
  
  gridded(idw) <- T
  
  Pg_hs_r <- raster(idw , "var1.pred") # create raster layer object with pig flock size data
  
  Pg_hs_r <- mask(Pg_hs_r , shp) # mask NA values in Raster object 
  
  pop <- resample(pop , Pg_h_r) # resample raster object: resample UGA pop data spatial data, by re-sampling to HH with pigs data
  
  Pg_den <- pop * 10 * Pg_h_r * Pg_hs_r # Equation as above for pig density: adjust pop * % HH keeping pigs * flock size (note: 10 because of the resolution change)
  
  #Pg_den <- pop * Pg_h_r * Pg_hs_r # Equation as above for pig density: adjust pop * % HH keeping pigs * flock size (note: 10 because of the resolution change)
  
  Pg_den_co <- Pg_den
  
  Pg_den_co[Pg_den_co > 100] <- NA # % over 100 is NA?
  
  return(list(Pg_den, Pg_den_co))

}

