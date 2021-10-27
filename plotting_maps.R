#===================================================================================================================#
#============================    Plotting maps of variables of interest (cluster level)       ======================#

# plot cluster-level datapoints for each variable #
plot_basic_distribution_variables_func <- function(admin_processed, clustered_data, variable_interest){
  
  if(variable_interest == "pg_h"){
    variable_interest = clustered_data$Pg_h
    legend_nm = c("proportion of HH with pigs by cluster")
  }
  
  
  visualise_cluster_data <- 
    ggplot() + 
    geom_polygon(data = admin_processed, aes(x = long, y = lat, group = group), fill = "grey80", colour = "grey90", alpha = 1)+ 
    geom_point(data = clustered_data, aes(x = x, y = y,color = variable_interest) )+ 
    scale_colour_gradientn(name = legend_nm, colours = heat.colors(12))+    
    coord_equal(ratio = 1)
  
  return(visualise_cluster_data)
  
}