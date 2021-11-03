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

#======================================================================#
#  Sanitation - percentage of household with low sanitation facilities #
#  The first map includes household who have no facilities at all      #

# A)
sanitation_cluster_plot <- plot_basic_distribution_variables_func(admin_processed = admin_2011_processed,
                                                                       clustered_data = My_cluter_data_2011, 
                                                                       variable_interest = "wc2")
sanitation_cluster_plot
# B)
spatial_sanitation <- process_plot_variogram_func(shp = shp_2011, 
                                               clustered_data = My_cluter_data_2011, 
                                               variable_interest = "wc2")
plot(spatial_sanitation[[4]]) # plot model fit to variogram data
# C)
kriged_sanitation <- kriging_func(spatial_object = spatial_sanitation, 
                               admin = admin_2011_processed, 
                               variable_interest = "wc2")
kriged_sanitation[[2]] # plot predicted variable across map
# D)
sanitation_processed <- spatially_process_variable_function(kriged = kriged_sanitation[[1]]) # final object: process - gridded object + rasterize

#======================================================================#
#  Sanitation - percentage of household with low sanitation facilities #
#  includes no sanitation or uncovered facilities                      #

# A)
sanitation2_cluster_plot <- plot_basic_distribution_variables_func(admin_processed = admin_2011_processed,
                                                                  clustered_data = My_cluter_data_2011, 
                                                                  variable_interest = "wc1")
sanitation2_cluster_plot
# B)
spatial_sanitation2 <- process_plot_variogram_func(shp = shp_2011, 
                                                  clustered_data = My_cluter_data_2011, 
                                                  variable_interest = "wc1")
plot(spatial_sanitation2[[4]]) # plot model fit to variogram data
# C)
kriged_sanitation2 <- kriging_func(spatial_object = spatial_sanitation2, 
                                  admin = admin_2011_processed, 
                                  variable_interest = "wc1")
kriged_sanitation2[[2]] # plot predicted variable across map
# D)
sanitation2_processed <- spatially_process_variable_function(kriged = kriged_sanitation2[[1]]) # final object: process - gridded object + rasterize

#=================================================================================##
##         poverty level                                                          ##
##  based on percentage of household in the poorest 20% (lowest socio-econ quintile)

# A)
poverty_cluster_plot <- plot_basic_distribution_variables_func(admin_processed = admin_2011_processed,
                                                                   clustered_data = My_cluter_data_2011, 
                                                                   variable_interest = "w1")
poverty_cluster_plot
# B)
spatial_poverty <- process_plot_variogram_func(shp = shp_2011, 
                                                   clustered_data = My_cluter_data_2011, 
                                                   variable_interest = "w1")
plot(spatial_poverty[[4]]) # plot model fit to variogram data
# C)
kriged_poverty <- kriging_func(spatial_object = spatial_poverty, 
                                   admin = admin_2011_processed, 
                                   variable_interest = "w1")
kriged_poverty[[2]] # plot predicted variable across map
# D)
poverty_processed <- spatially_process_variable_function(kriged = kriged_poverty[[1]]) # final object: process - gridded object + rasterize

#=================================================================================##
##         poverty level                                                          ##
##  based on percentage of household in bottom two socio-economic quintiles       ##

# A)
poverty2_cluster_plot <- plot_basic_distribution_variables_func(admin_processed = admin_2011_processed,
                                                               clustered_data = My_cluter_data_2011, 
                                                               variable_interest = "w1&w2")
poverty2_cluster_plot
# B)
spatial_poverty2 <- process_plot_variogram_func(shp = shp_2011, 
                                               clustered_data = My_cluter_data_2011, 
                                               variable_interest = "w1&w2")
plot(spatial_poverty2[[4]]) # plot model fit to variogram data
# C)
kriged_poverty2 <- kriging_func(spatial_object = spatial_poverty2, 
                               admin = admin_2011_processed, 
                               variable_interest = "w1&w2")
kriged_poverty2[[2]] # plot predicted variable across map
# D)
poverty2_processed <- spatially_process_variable_function(kriged = kriged_poverty2[[1]]) # final object: process - gridded object + rasterize

#===================================================================================================#
# Investigating the distribution of the key variables to define meaningful breaks into high and low #

# 1) livestock ownership distribution
# for the final map we are using the 3rd quantile, 
# i.e. the map shows were the 25% clusters with most households with pigs are. 
livestock_ownership <- distribution_variables_definecutoff_func(spatial_variable = pig_pop_density[[3]],
                                                                cutoff_value = 0.23)
plot(livestock_ownership[[1]]) # histogram
livestock_ownership[[2]] # distribution properties
plot(livestock_ownership[[3]]) # plot high vs low
plot(admin_2011_processed,add = T)

# 2) pig population 
# based on DHS computation
pig_populationDHS <- distribution_variables_definecutoff_func(spatial_variable = pig_pop_density[[1]],
                                                                cutoff_value = 1)
plot(pig_populationDHS[[1]])
pig_populationDHS[[2]]
plot(pig_populationDHS[[3]])
plot(admin_2011_processed,add = T)

# 3) based on FAO Robinson layer (cannot view historgram etc so use diff function)
pig_populationFAO <- distribution_variables_definecutoff_func2(spatial_variable = pig_pop,
                                                              cutoff_value = 1)
plot(pig_populationFAO[[1]]) # cannot review histogram & summarise distribution
plot(admin_2011_processed,add = T)

# 4) pig extensive rural system 
pig_roam <- distribution_variables_definecutoff_func(spatial_variable = pig_roam_processed,
                                                              cutoff_value = 0.153) # cut-off indicated as 0.153 value (not 3rd quartile)
plot(pig_roam[[1]])
pig_roam[[2]]
plot(pig_roam[[3]])
plot(admin_2011_processed,add = T)

# 5) poverty the 40% poorest
poverty_lowest40 <- distribution_variables_definecutoff_func(spatial_variable = poverty2_processed,
                                                     cutoff_value = 0.666333) # 3rd quartile cut-off value (0.666 in original)
plot(poverty_lowest40[[1]])
poverty_lowest40[[2]]
plot(poverty_lowest40[[3]])
plot(admin_2011_processed,add = T)

# 6) poverty the 20% poorest
poverty_lowest20 <- distribution_variables_definecutoff_func(spatial_variable = poverty_processed,
                                                             cutoff_value = 0.43828) # 3rd quartile cut-off value (0.438 used in original)
plot(poverty_lowest20[[1]])
poverty_lowest20[[2]]
plot(poverty_lowest20[[3]])
plot(admin_2011_processed,add = T)

# 7) uncovered sanitation
sanitation <- distribution_variables_definecutoff_func(spatial_variable = sanitation2_processed,
                                                             cutoff_value = 0.40869) # 3rd quartile cut-off value (0.408 used in original)
plot(sanitation[[1]])
sanitation[[2]]
plot(sanitation[[3]])
plot(admin_2011_processed,add = T)

#===================================================================================================#
#                            Overlays of risk factors                                               #

# Combined risk factor scoring system (risk factor 1 + (risk factor 2 * 1.01)+ (risk factor 1 * 1.1)) :
# 1. 0 = ABC low
# 2. 1 = A high (yellow)
# 3. 1.01 = B high (blue)
# 4. 1.1 = C high (red)
# 5. 2.01 = AB high (green)
# 6. 2.1 = AC high (orange)
# 7. 2.11 = BC high (purple)
# 8. 3.11 = ABC high (brown)

# 1) proportion of households with bad sanitation, proportion of extensive households and poverty 
overlay1 <- plotting_overlays_func(risk_factor1 = sanitation[[3]], risk_factor3 = poverty_lowest40[[3]], risk_factor2 = pig_roam[[3]],
                       admin_processed = admin_2011_processed, admin = admin_2011, year = "2011") 
overlay1[[2]] # ggplot map with discrete risk factors

# 2) proportion of household with bad sanitaiton (1), pig distribution (Robinson; 1.01) & proportion of extensive households (1.1)  
overlay2 <- plotting_overlays_func(risk_factor1 = sanitation[[3]], risk_factor2 = pig_populationFAO[[1]], risk_factor3 = pig_roam[[3]],
                                   admin_processed = admin_2011_processed, admin = admin_2011, year = "2011") 
overlay2[[2]]

# 3) proportion of household with bad sanitaiton (1),poorest 20% (1.01) and proportion of extensive households (1.1) 
overlay3 <- plotting_overlays_func(risk_factor1 = sanitation[[3]], risk_factor2 = poverty_lowest20[[3]], risk_factor3 = pig_roam[[3]],
                                   admin_processed = admin_2011_processed, admin = admin_2011, year = "2011") 
overlay3[[2]]

# 4) proportion of HH with bad sanitation (1), pig distribution (DHS data; 1.01) and proportion of extensive households (1.1) 
overlay4 <- plotting_overlays_func(risk_factor1 = sanitation[[3]], risk_factor2 = pig_populationDHS[[3]], risk_factor3 = pig_roam[[3]],
                                   admin_processed = admin_2011_processed, admin = admin_2011, year = "2011") 
overlay4[[2]]

# 5) FINAL RISK MAP; proportion of HH with bad sanitation (1), pig FAO (Robinson) distribution (1.01) and poorest 40% (1.1)  
overlay5 <- plotting_overlays_func(risk_factor1 = sanitation[[3]], risk_factor2 = pig_populationFAO[[1]], risk_factor3 = poverty_lowest40[[3]],
                                   admin_processed = admin_2011_processed, admin = admin_2011, year = "2011") 
overlay5[[2]]

