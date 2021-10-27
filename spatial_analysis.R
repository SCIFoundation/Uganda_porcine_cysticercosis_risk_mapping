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


# This function assesses flock sizes from households who have reported owning pigs 

pig_flock_size_HH_kriged_func <- function(admin_processed, clustered_data, spatial_object, variable_interest) {
  
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
  
  idw = idw(Pg_hs ~ 1, geodata[!is.na(geodata$Pg_hs),], grid, maxdist = Inf, idp = 2) # krigin flock/herd sizes using [inverse distance weighted interpolation]
  
  pig_flock_plot2 <- 
    ggplot(data = idw@data, aes(x = idw$x1, y = idw$x2)) + 
    geom_tile(aes(fill = idw$var1.pred)) + 
    coord_equal() +
    scale_fill_gradient(low = "yellow", high = "red") +
    geom_polygon(data = admin_processed, aes(x = long, y = lat, group = group), color = "black", alpha = 0) + 
    theme_bw()
  
  return(list(idw, pig_flock_plot1, pig_flock_plot2))
  
}