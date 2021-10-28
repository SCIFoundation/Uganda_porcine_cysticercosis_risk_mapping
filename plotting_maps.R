#===================================================================================================================#
#============================    Plotting maps of variables of interest (cluster level)       ======================#

# plot cluster-level datapoints for each variable #
plot_basic_distribution_variables_func <- function(admin_processed, clustered_data, variable_interest){
  
  if(variable_interest == "pg_h"){
    variable_interest = clustered_data$Pg_h
    legend_nm = c("proportion of HH with pigs by cluster")
  }
  
  if(variable_interest == "pg_2"){
    variable_interest = clustered_data$Pg_h
    legend_nm = c("proportion of free-roaming systems")
  }
  
  
  visualise_cluster_data <- 
    ggplot() + 
    geom_polygon(data = admin_processed, aes(x = long, y = lat, group = group), fill = "grey80", colour = "grey90", alpha = 1)+ 
    geom_point(data = clustered_data, aes(x = x, y = y,color = variable_interest) )+ 
    scale_colour_gradientn(name = legend_nm, colours = heat.colors(12))+    
    coord_equal(ratio = 1)
  
  return(visualise_cluster_data)
  
}


# this function plots the Robinson et al. pig pop density map & stores plot object

plot_pig_density_map_func <- function(ppop, admin_processed) {
  
  ppop <- raster('ppop.tif') # need to make sure the ppop.tif file is in the working directory (pig pop density image): ok to keep this as 2006 data?
  
  breakpoints <- c(0, 1, 2, 5, 10, 25, 75, 250,1000, 3000)
  
  colors <- c("gray90", rainbow_hcl(7))
  
  plot(ppop, breaks = breakpoints, col = colors)
  
  p = plot(admin_processed, add = T)
  
  return(p)

}

