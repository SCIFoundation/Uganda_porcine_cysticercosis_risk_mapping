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


# function for plotting raster layer in ggplot2
# # see: https://stackoverflow.com/questions/47116217/overlay-raster-layer-on-map-in-ggplot2-in-r 

gplot_data <- function(x, maxpixels = 50000)  {
  
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  
  ## Extract values
  
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  
  names(dat) <- c('value', 'variable')
  
  dat <- tibble::as_tibble(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}

# function for plotting overlay of risk factors (combined risk factor score) #

plotting_overlays_func <- function(risk_factor1, risk_factor2, risk_factor3, admin_processed, admin, year){

  # set-up plot breakpoints in base R
  breakpoints <- c(0, 0.99, 1.0099, 1.099, 2.0099, 2.099, 2.1099, 3.1099, 3.12)
  
  breakpoints2 <- c('all low', 'A' , 'B', 'C', 'AB', 'AC', 'BC', 'ABC')
  
  colors <- c( "gray90", "gold", "steelblue1", "red", "springgreen", 'orange', 'plum1', "brown")
  
  # combined metric with 3 risk factors (e.g. rf1 + rf2 + rf3)
  
  # * A proportion of households bad sanitation high =1 is low = 0  
  # * B proportion of is high =1.1 is low =0 
  # * C proportion of poor 40% is high = 1.01 is low zero 

  over <- risk_factor1 + risk_factor2 * 1.1 + risk_factor3 * 1.01
  over <- mask(over , admin)
  
  # plot in base r 
  plot(over , breaks = breakpoints, col = colors, legend = FALSE)
  legend("bottomleft", inset = .02, title = "classes", legend = breakpoints2
       , fill = colors, horiz = F, cex = 0.8)
  plot(admin, add = T)
  
  # plot risk map in ggplot2 #
  riskfact_df <- gplot_data(over)
  
  riskfact_df <- 
    riskfact_df %>%
    mutate(risk_fact_bins = cut(value,
                              breaks = breakpoints, right = FALSE,
                              labels = c('all low','A' , 'B', 'C', 'AB', 'AC', 'BC', 'ABC')))
  
  # continuous risk factor plot #
  a <- ggplot() +
    
    geom_tile(data = riskfact_df, 
            aes(x = x, y = y, fill = value)) +
    geom_polygon(data=admin_processed, aes(x=long, y=lat, group=group), color="black", alpha=0) +
    scale_fill_gradient("classes",
                      low = 'yellow', high = 'blue',
                      na.value = NA) +
    coord_equal() +
    theme_bw()
  
  # FINAL plot: discrete risk factor plot # 
  # discrete color gradient from continous color gradient in ggplot: https://stackoverflow.com/questions/64762377/create-r-ggplot2-discrete-colour-palette-for-gradient-map-with-continuous-values
  # convert cont to discrete : https://www.r-bloggers.com/2020/09/how-to-convert-continuous-variables-into-categorical-by-creating-bins/
  
  value_rf <- c("grey90", "gold", "royalblue1", "red1", "springgreen", "darkorange1", "pink1", "tan4")
  
  b <- ggplot() +
    geom_tile(data = riskfact_df, 
            aes(x = x, y = y, fill = risk_fact_bins)) +
    geom_polygon(data=admin_processed, aes(x=long, y=lat, group=group), color="black", alpha=0) +
    scale_fill_manual(name = "Class",
                    values = value_rf, 
                    na.value = NA, na.translate = FALSE)+
    coord_equal() +
    theme_bw() +
    theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(),
        panel.background = element_blank()) +
    panel_border(remove = TRUE) +
    annotate("text", label = year, x = 29.75, y = 3.8, size = 7, colour = "black")
  
  return(list(a, b))


}