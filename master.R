#====================================================================================================================#
#===================================              MASTER SCRIPT            ==========================================#
#====================================================================================================================#

rm(list = ls())

#=============================================#
#       Load data files                       #
#=============================================#

#=======================#                                                                                             
# Initiatie sub-scripts #                                                                                             
source('libraries.R')
source('data_processing_extraction.R')
source('Spatial_analysis.R')
source('plotting_maps.R')

#===========================#
#===========================#
# Load data & mapping files #

# 2011-2012 risk map underlying data sources #

DHS_2011_data <- read.spss('~/Uganda_porcine_cysticercosis_risk_mapping/data/UGA_2012_60/HR.SAV', to.data.frame = T,use.value.labels = FALSE)

geo_2011 <- readShapePoints('~/Uganda_porcine_cysticercosis_risk_mapping/data/UGA_2012_60/GE.shp') # shape file of datapoints from DHS

shp_2011 <- readShapePoly('~/Uganda_porcine_cysticercosis_risk_mapping/data/UGA_2012_60/sdr_subnational_boundaries.shp')

admin_2011 <- readShapePoly('~/Uganda_porcine_cysticercosis_risk_mapping/data/admin-gdam/UGA_adm1.shp') # source: https://gadm.org/download_country_v3.html

human_pop <- raster('UGA_ppp_v2b_2015.tif') # human pop data - from world pop 2015 unadjusted; loads and create raster layer object

pig_pop <- raster('ppop.tif') # pig pop data - robinson 2006 (pig densty map); loads and create raster layer object

#======================================#
#==== Data cleaning & processing ======#

admin_2011_processed <- processing_admin_data_func(admin = admin_2011) # process admin (district) spatial file
shp_2011_processed <- processing_subnational_shp_func(shp = shp_2011) # process subnational boundaries shp file

#We are extracting the DHS standard survey of 2011. We have extracted
# 1. toilet that a household use
# 2. number of livestock owned by the households
# 3. a couple of household charateristics
My_cluter_data_2011 <- DHS2011_extraction_tocluster_func(data = DHS_2011_data, geo = geo_2011)


#==========================================================================#
#========  Spatial (spatial autocorrelation & kriging analysis ============#

#============================#
# 1) households with pigs    #

# proportion of HH with pigs by cluster (0-1)
hh_pigs_cluster_plot <- plot_basic_distribution_variables_func(admin_processed = admin_2011_processed,
                                                               clustered_data = My_cluter_data_2011, 
                                                               variable_interest = "pg_h")
hh_pigs_cluster_plot

# plot semi-variogram (spatial autocorrelation) of variable of interest (model vs data)
spatial_pigHH <- process_plot_variogram_func(shp = shp_2011, 
                                             clustered_data = My_cluter_data_2011, 
                                             variable_interest = "pg_h")

plot(spatial_pigHH[[4]]) # plot model fit to variogram data

# spatial interopolation using ordinary kriging (with the underlying variogram) - this may take a few minutes
kriged_pigHH <- kriging_func(spatial_object = spatial_pigHH, 
                             admin = admin_2011_processed, 
                             variable_interest = "pg_h")

kriged_pigHH[[2]] # plot predicted variable across map

# specific for pig population : look at flock sizes (variable = pg_hs)
# this uses Inverse Distance Weighing (deterministic, not a geospatial approach) for spatial interopolation
# as we do not compute the varigogram

pig_flocksize_interpolated <- pig_flock_size_HH_interpolated_func(admin_processed = admin_2011_processed, 
                                                      clustered_data = My_cluter_data_2011, 
                                                      spatial_object = spatial_pigHH,
                                                      variable_interest = "pg_hs")

pig_flocksize_interpolated[[2]] # plot avg. pig herd size spatial datapoint
pig_flocksize_interpolated[[3]] # plot predicted pig flock size 

# Compute and plot pig density (approach #1)
# based on human population (from world pop, closest year) x percent of household keeping pigs x flock size
pig_pop_density <- compute_pig_denisty_func(pop = human_pop, kriged = kriged_pigHH[[1]], 
                                            shp = shp_2011, idw = pig_flocksize_interpolated[[1]])

plot(pig_pop_density[[1]]) # pig density (shows minimal variation)
plot(pig_pop_density[[2]]) # pig density (> 50 or 100 values excluded re high density urban areas to show more variation)

# Pig density maps do not show much variation because there is an overestimation of pig density in urban areas. 
# So we re-run it showing taking the extreme value out, i.e values above 50. Then a pattern appears.
# Next, we compare it with the FAO Robinson pig distribution map - approach #2

pig_pop_density_FAOmap <- plot_pig_density_map_func(ppop = pig_pop, admin_processed = admin_2011_processed)
pig_pop_density_FAOmap

#==========================================================================#
# 2) presence of extensive (pig roaming = high potential for transmission) #

# This these are places where people let there pigs roam and therefore there is high potential for diseases 
# pg_2 variable = household with 1 or 2 pigs and RURAL
# same steps as above - A) plot cluster-level datapoints, B) process variogram (spatial autocorrelation), 
# C) spatial interpolation (ordinary kriging - using underlying variogram) & plot D) process - gridded + rasterize

# A)
free_range_pigs_cluster_plot <- plot_basic_distribution_variables_func(admin_processed = admin_2011_processed,
                                                               clustered_data = My_cluter_data_2011, 
                                                               variable_interest = "pg_2")
free_range_pigs_cluster_plot

# B)
spatial_pigroam <- process_plot_variogram_func(shp = shp_2011, 
                                             clustered_data = My_cluter_data_2011, 
                                             variable_interest = "pg_2")

plot(spatial_pigroam[[4]]) # plot model fit to variogram data


# C)
kriged_pigroam <- kriging_func(spatial_object = spatial_pigroam, 
                             admin = admin_2011_processed, 
                             variable_interest = "pg_2")

kriged_pigroam[[2]] # plot predicted variable across map

# D)
pig_roam_processed <- spatially_process_variable_function(kriged = kriged_pigroam[[1]]) # final object: process - gridded object + rasterize


