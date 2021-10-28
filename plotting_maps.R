#===================================================================================================================#
#============================    Plotting maps of variables of interest (cluster level)       ======================#

# plot cluster-level datapoints for each variable #
plot_basic_distribution_variables_func <- function(admin_processed, clustered_data, variable_interest){
  
  variable_interest_tocheck <- variable_interest
  
  if(variable_interest_tocheck == "pg_h"){
    variable_interest = clustered_data$Pg_h
    legend_nm = c("proportion of HH with pigs by cluster")
  }
  
  if(variable_interest_tocheck == "pg_2"){
    variable_interest = clustered_data$Pg_2
    legend_nm = c("proportion of free-roaming systems")
  }
  
  if(variable_interest_tocheck == "wc2"){
    variable_interest = clustered_data$wc2
    legend_nm = c("proportion of HH w/ low sanitation")
  }
  
  if(variable_interest_tocheck == "wc1"){
    variable_interest = clustered_data$wc1
    legend_nm = c("proportion of HH w/ low sanitation 
(inclduing uncovered)")
  }
  
  if(variable_interest_tocheck == "w1"){
    variable_interest = clustered_data$w1
    legend_nm = c("proportion of HH in lowest 
socio-economic quintile")
  }
  
  if(variable_interest_tocheck == "w1&w2"){
    variable_interest = c(clustered_data$w1 + clustered_data$w2)
    legend_nm = c("proportion of HH in lowest two
socio-economic quintile") # percentage of poor in the 40% percent poorest category (poorest and poor) i.e. w1 and w2
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

