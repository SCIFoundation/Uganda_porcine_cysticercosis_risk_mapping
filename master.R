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
My_cluster_data_2011 <- DHS_extraction_tocluster_func(data = DHS_2011_data, geo = geo_2011, year = 2011)

#==========================================================================#
#========  Spatial (spatial autocorrelation & kriging analysis ============#

#============================#
# 1) households with pigs    #

# proportion of HH with pigs by cluster (0-1)
hh_pigs_cluster_plot <- plot_basic_distribution_variables_func(admin_processed = admin_2011_processed,
                                                               clustered_data = My_cluster_data_2011, 
                                                               variable_interest = "pg_h")
hh_pigs_cluster_plot

# plot semi-variogram (spatial autocorrelation) of variable of interest (model vs data)
spatial_pigHH <- process_plot_variogram_func(shp = shp_2011, 
                                             clustered_data = My_cluster_data_2011, 
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
                                                      clustered_data = My_cluster_data_2011, 
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
                                                               clustered_data = My_cluster_data_2011, 
                                                               variable_interest = "pg_2")
free_range_pigs_cluster_plot
# B)
spatial_pigroam <- process_plot_variogram_func(shp = shp_2011, 
                                             clustered_data = My_cluster_data_2011, 
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
                                                                       clustered_data = My_cluster_data_2011, 
                                                                       variable_interest = "wc2")
sanitation_cluster_plot
# B)
spatial_sanitation <- process_plot_variogram_func(shp = shp_2011, 
                                               clustered_data = My_cluster_data_2011, 
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
                                                                  clustered_data = My_cluster_data_2011, 
                                                                  variable_interest = "wc1")
sanitation2_cluster_plot
# B)
spatial_sanitation2 <- process_plot_variogram_func(shp = shp_2011, 
                                                  clustered_data = My_cluster_data_2011, 
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
                                                                   clustered_data = My_cluster_data_2011, 
                                                                   variable_interest = "w1")
poverty_cluster_plot
# B)
spatial_poverty <- process_plot_variogram_func(shp = shp_2011, 
                                                   clustered_data = My_cluster_data_2011, 
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
                                                               clustered_data = My_cluster_data_2011, 
                                                               variable_interest = "w1&w2")
poverty2_cluster_plot
# B)
spatial_poverty2 <- process_plot_variogram_func(shp = shp_2011, 
                                               clustered_data = My_cluster_data_2011, 
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
overlay5_2011 <- plotting_overlays_func(risk_factor1 = sanitation[[3]], risk_factor2 = pig_populationFAO[[1]], risk_factor3 = poverty_lowest40[[3]],
                                   admin_processed = admin_2011_processed, admin = admin_2011, year = "2011") 
overlay5_2011[[2]]


#====================================================================================================================#
#                                       2001 PCC risk mapping                                                        #
#====================================================================================================================#

#===========================#
# Load data & mapping files #
DHS_2001_data <- read.spss('~/Uganda_porcine_cysticercosis_risk_mapping/data/DHS 2001/UGHR41FL.SAV', to.data.frame = T,use.value.labels = FALSE)
geo_2001 <- readShapePoints("~/Uganda_porcine_cysticercosis_risk_mapping/data/DHS 2001/UGGE43FL.shp") 
shp_2001 <- readShapePoly('~/Uganda_porcine_cysticercosis_risk_mapping/data/UGA_2012_60/sdr_subnational_boundaries.shp') # use 2011 boundries ??
admin_2001 <- readShapePoly('~/Uganda_porcine_cysticercosis_risk_mapping/data/admin-gdam/UGA_adm1.shp') # source: https://gadm.org/download_country_v3.html
# pig_pop <- raster('ppop.tif') # pig pop data - robinson 2006 (pig densty map); loads and create raster layer object

#======================================#
#==== Data cleaning & processing ======#

admin_2001_processed <- processing_admin_data_func(admin = admin_2001) # process admin (district) spatial file
shp_2001_processed <- processing_subnational_shp_func(shp = shp_2001) # process subnational boundaries shp file

# data extraction & cleaning for 2001 DHS data variables of interest #
# NOTE: no pig distribution data or no latrine only (wc1) data variables in 2001 DHS
My_cluster_data_2001 <- DHS_extraction_tocluster_func(data = DHS_2001_data, geo = geo_2001, year = 2001) # 

#======================================================================#
#  Sanitation - percentage of household with low sanitation facilities #
#  The first map includes household who have no facilities at all      #

# A)
sanitation_cluster_plot_2001 <- plot_basic_distribution_variables_func(admin_processed = admin_2001_processed,
                                                                  clustered_data = My_cluster_data_2001, 
                                                                  variable_interest = "wc2")
sanitation_cluster_plot_2001
# B)
spatial_sanitation_2001 <- process_plot_variogram_func(shp = shp_2001, 
                                                  clustered_data = My_cluster_data_2001, 
                                                  variable_interest = "wc2")
plot(spatial_sanitation_2001[[4]]) # plot model fit to variogram data
# C)
kriged_sanitation_2001 <- kriging_func(spatial_object = spatial_sanitation_2001, 
                                  admin = admin_2001_processed, 
                                  variable_interest = "wc2")
kriged_sanitation_2001[[2]] # plot predicted variable across map
# D)
sanitation_processed_2001 <- spatially_process_variable_function(kriged = kriged_sanitation_2001[[1]]) # final object: process - gridded object + rasterize

#=================================================================================##
##         poverty level                                                          ##
##  based on percentage of household in the poorest 20% (lowest socio-econ quintile)

# A)
poverty_cluster_plot_2001 <- plot_basic_distribution_variables_func(admin_processed = admin_2001_processed,
                                                               clustered_data = My_cluster_data_2001, 
                                                               variable_interest = "w1")
poverty_cluster_plot_2001
# B)
spatial_poverty_2001 <- process_plot_variogram_func(shp = shp_2001, 
                                               clustered_data = My_cluster_data_2001, 
                                               variable_interest = "w1")
plot(spatial_poverty_2001[[4]]) # plot model fit to variogram data
# C)
kriged_poverty_2001 <- kriging_func(spatial_object = spatial_poverty_2001, 
                               admin = admin_2001_processed, 
                               variable_interest = "w1")
kriged_poverty_2001[[2]] # plot predicted variable across map
# D)
poverty_processed_2001 <- spatially_process_variable_function(kriged = kriged_poverty_2001[[1]]) # final object: process - gridded object + rasterize

#=================================================================================##
##         poverty level                                                          ##
##  based on percentage of household in bottom two socio-economic quintiles       ##

# A)
poverty2_cluster_plot_2001 <- plot_basic_distribution_variables_func(admin_processed = admin_2001_processed,
                                                                clustered_data = My_cluster_data_2001, 
                                                                variable_interest = "w1&w2")
poverty2_cluster_plot_2001
# B)
spatial_poverty2_2001 <- process_plot_variogram_func(shp = shp_2001, 
                                                clustered_data = My_cluster_data_2001, 
                                                variable_interest = "w1&w2")
plot(spatial_poverty2_2001[[4]]) # plot model fit to variogram data
# C)
kriged_poverty2_2001 <- kriging_func(spatial_object = spatial_poverty2_2001, 
                                admin = admin_2001_processed, 
                                variable_interest = "w1&w2")
kriged_poverty2_2001[[2]] # plot predicted variable across map
# D)
poverty2_processed_2001 <- spatially_process_variable_function(kriged = kriged_poverty2_2001[[1]]) # final object: process - gridded object + rasterize

# we re-use the FAO (Robinson et al.) pig density map 
pig_pop_density_FAOmap <- plot_pig_density_map_func(ppop = pig_pop, admin_processed = admin_2001_processed)
pig_pop_density_FAOmap

#===================================================================================================#
# Investigating the distribution of the key variables to define meaningful breaks into high and low #

# 3) based on FAO Robinson layer (cannot view historgram etc so use diff function)
pig_populationFAO <- distribution_variables_definecutoff_func2(spatial_variable = pig_pop,
                                                               cutoff_value = 1)
plot(pig_populationFAO[[1]]) # cannot review histogram & summarise distribution
plot(admin_2001_processed,add = T)

# 5) poverty the 40% poorest
poverty_lowest40_2001 <- distribution_variables_definecutoff_func(spatial_variable = poverty2_processed_2001,
                                                             cutoff_value = 0.3383) # 3rd quartile cut-off value 
plot(poverty_lowest40_2001[[1]])
poverty_lowest40_2001[[2]]
plot(poverty_lowest40_2001[[3]])
plot(admin_2001_processed,add = T)

# 6) uncovered sanitation
sanitation_2001 <- distribution_variables_definecutoff_func(spatial_variable = sanitation_processed_2001,
                                                       cutoff_value = 0.384284) # 3rd quartile cut-off value 
plot(sanitation_2001[[1]])
sanitation_2001[[2]]
plot(sanitation_2001[[3]])
plot(admin_2001_processed,add = T)

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

# 5) FINAL RISK MAP; proportion of HH with bad sanitation (1), pig FAO (Robinson) distribution (1.01) and poorest 40% (1.1)  
overlay5_2001 <- plotting_overlays_func(risk_factor1 = sanitation_2001[[3]], risk_factor2 = pig_populationFAO[[1]], risk_factor3 = poverty_lowest40_2001[[3]],
                                   admin_processed = admin_2001_processed, admin = admin_2001, year = "2001") 
overlay5_2001[[2]]


#====================================================================================================================#
#                                       2006 PCC risk mapping                                                        #
#====================================================================================================================#
#===========================#
# Load data & mapping files #
DHS_2006_data <- read.spss('~/Uganda_porcine_cysticercosis_risk_mapping/data/DHS 2006/UGHR52FL.SAV', to.data.frame = T,use.value.labels = FALSE)
geo_2006 <- readShapePoints("~/Uganda_porcine_cysticercosis_risk_mapping/data/DHS 2006/UGGE53FL.shp") # use 2006 boundaries for 2001 (unable to locate 2001 boundaries)
shp_2006 <- readShapePoly('~/Uganda_porcine_cysticercosis_risk_mapping/data/UGA_2012_60/sdr_subnational_boundaries.shp') # use 2011 boundries ??
admin_2006 <- readShapePoly('~/Uganda_porcine_cysticercosis_risk_mapping/data/admin-gdam/UGA_adm1.shp') # source: https://gadm.org/download_country_v3.html
districts_2006 <- readShapePoly("~/Uganda_porcine_cysticercosis_risk_mapping/data/District boundaries 2006/Stanford -2010/vg894mz3698.shp") # need this for 2006
human_pop_2006 <- raster('~/Uganda_porcine_cysticercosis_risk_mapping/data/DHS 2006/uga_pd_2006_1km.tif') # human pop data - from world pop 2015 unadjusted; loads and create raster layer object

# pig_pop <- raster('ppop.tif') # pig pop data - robinson 2006 (pig densty map); loads and create raster layer object

#======================================#
#==== Data cleaning & processing ======#

admin_2006_processed <- processing_admin_data_func(admin = admin_2006) # process admin (district) spatial file
shp_2006_processed <- processing_subnational_shp_func(shp = shp_2006) # process subnational boundaries shp file

# data extraction & cleaning for 2001 DHS data variables of interest #
My_cluster_data_2006 <- DHS_extraction_tocluster_func(data = DHS_2006_data, geo = geo_2006, year = 2006) # 


#==========================================================================#
#========  Spatial (spatial autocorrelation & kriging analysis ============#

#============================#
# 1) households with pigs    #

# proportion of HH with pigs by cluster (0-1)
hh_pigs_cluster_plot_2006 <- plot_basic_distribution_variables_func(admin_processed = admin_2006_processed,
                                                               clustered_data = My_cluster_data_2006, 
                                                               variable_interest = "pg_h")
hh_pigs_cluster_plot_2006

# plot semi-variogram (spatial autocorrelation) of variable of interest (model vs data)
spatial_pigHH_2006 <- process_plot_variogram_func(shp = shp_2006, 
                                             clustered_data = My_cluster_data_2006, 
                                             variable_interest = "pg_h")

plot(spatial_pigHH_2006[[4]]) # plot model fit to variogram data

# spatial interopolation using ordinary kriging (with the underlying variogram) - this may take a few minutes
kriged_pigHH_2006 <- kriging_func(spatial_object = spatial_pigHH_2006, 
                             admin = admin_2006_processed, 
                             variable_interest = "pg_h")

kriged_pigHH_2006[[2]] # plot predicted variable across map

# specific for pig population : look at flock sizes (variable = pg_hs)
# this uses Inverse Distance Weighing (deterministic, not a geospatial approach) for spatial interopolation
pig_flocksize_interpolated_2006 <- pig_flock_size_HH_interpolated_func(admin_processed = admin_2006_processed, 
                                                                  clustered_data = My_cluster_data_2006, 
                                                                  spatial_object = spatial_pigHH_2006,
                                                                  variable_interest = "pg_hs")

pig_flocksize_interpolated_2006[[2]] # plot avg. pig herd size spatial datapoint
pig_flocksize_interpolated_2006[[3]] # plot predicted pig flock size 

# Compute and plot pig density (approach #1) : based on human population (from world pop, closest year) x percent of household keeping pigs x flock size
pig_pop_density_2006 <- compute_pig_denisty_func(pop = human_pop_2006, kriged = kriged_pigHH_2006[[1]], 
                                            shp = shp_2006, idw = pig_flocksize_interpolated_2006[[1]])

plot(pig_pop_density_2006[[1]]) # pig density (shows minimal variation)
plot(pig_pop_density_2006[[2]]) # pig density (> 50 or 100 values excluded re high density urban areas to show more variation)

# Next, we compare it with the FAO Robinson pig distribution map - approach #2
pig_pop_density_FAOmap_2006 <- plot_pig_density_map_func(ppop = pig_pop, admin_processed = admin_2006_processed)
pig_pop_density_FAOmap_2006

#==========================================================================#
# 2) presence of extensive (pig roaming = high potential for transmission) #

# A)
free_range_pigs_cluster_plot_2006 <- plot_basic_distribution_variables_func(admin_processed = admin_2006_processed,
                                                                       clustered_data = My_cluster_data_2006, 
                                                                       variable_interest = "pg_2")
free_range_pigs_cluster_plot_2006
# B)
spatial_pigroam_2006 <- process_plot_variogram_func(shp = shp_2006, 
                                               clustered_data = My_cluster_data_2006, 
                                               variable_interest = "pg_2")
plot(spatial_pigroam_2006[[4]]) # plot model fit to variogram data
# C)
kriged_pigroam_2006 <- kriging_func(spatial_object = spatial_pigroam_2006, 
                               admin = admin_2006_processed, 
                               variable_interest = "pg_2")
kriged_pigroam_2006[[2]] # plot predicted variable across map
# D)
pig_roam_processed_2006 <- spatially_process_variable_function(kriged = kriged_pigroam_2006[[1]]) # final object: process - gridded object + rasterize

#======================================================================#
#  Sanitation - percentage of household with low sanitation facilities #
#  The first map includes household who have no facilities at all      #

# A)
sanitation_cluster_plot_2006 <- plot_basic_distribution_variables_func(admin_processed = admin_2006_processed,
                                                                  clustered_data = My_cluster_data_2006, 
                                                                  variable_interest = "wc2")
sanitation_cluster_plot_2006
# B)
spatial_sanitation_2006 <- process_plot_variogram_func(shp = shp_2006, 
                                                  clustered_data = My_cluster_data_2006, 
                                                  variable_interest = "wc2")
plot(spatial_sanitation_2006[[4]]) # plot model fit to variogram data
# C)
kriged_sanitation_2006 <- kriging_func(spatial_object = spatial_sanitation_2006, 
                                  admin = admin_2006_processed, 
                                  variable_interest = "wc2")
kriged_sanitation_2006[[2]] # plot predicted variable across map
# D)
sanitation_processed_2006 <- spatially_process_variable_function(kriged = kriged_sanitation_2006[[1]]) # final object: process - gridded object + rasterize

#======================================================================#
#  Sanitation - percentage of household with low sanitation facilities #
#  includes no sanitation or uncovered facilities                      #

# A)
sanitation2_cluster_plot_2006 <- plot_basic_distribution_variables_func(admin_processed = admin_2006_processed,
                                                                   clustered_data = My_cluster_data_2006, 
                                                                   variable_interest = "wc1")
sanitation2_cluster_plot_2006
# B)
spatial_sanitation2_2006 <- process_plot_variogram_func(shp = shp_2006, 
                                                   clustered_data = My_cluster_data_2006, 
                                                   variable_interest = "wc1")
plot(spatial_sanitation2_2006[[4]]) # plot model fit to variogram data
# C)
kriged_sanitation2_2006 <- kriging_func(spatial_object = spatial_sanitation2_2006, 
                                   admin = admin_2006_processed, 
                                   variable_interest = "wc1")
kriged_sanitation2_2006[[2]] # plot predicted variable across map
# D)
sanitation2_processed_2006 <- spatially_process_variable_function(kriged = kriged_sanitation2_2006[[1]]) # final object: process - gridded object + rasterize

#=================================================================================##
##         poverty level                                                          ##
##  based on percentage of household in the poorest 20% (lowest socio-econ quintile)

# A)
poverty_cluster_plot_2006 <- plot_basic_distribution_variables_func(admin_processed = admin_2006_processed,
                                                               clustered_data = My_cluster_data_2006, 
                                                               variable_interest = "w1")
poverty_cluster_plot_2006
# B)
spatial_poverty_2006 <- process_plot_variogram_func(shp = shp_2006, 
                                               clustered_data = My_cluster_data_2006, 
                                               variable_interest = "w1")
plot(spatial_poverty_2006[[4]]) # plot model fit to variogram data
# C)
kriged_poverty_2006 <- kriging_func(spatial_object = spatial_poverty_2006, 
                               admin = admin_2006_processed, 
                               variable_interest = "w1")
kriged_poverty_2006[[2]] # plot predicted variable across map
# D)
poverty_processed_2006 <- spatially_process_variable_function(kriged = kriged_poverty_2006[[1]]) # final object: process - gridded object + rasterize

#=================================================================================##
##         poverty level  # 2                                                     ##
##  based on percentage of household in bottom two socio-economic quintiles       ##

# A)
poverty2_cluster_plot_2006 <- plot_basic_distribution_variables_func(admin_processed = admin_2006_processed,
                                                                clustered_data = My_cluster_data_2006, 
                                                                variable_interest = "w1&w2")
poverty2_cluster_plot_2006
# B)
spatial_poverty2_2006 <- process_plot_variogram_func(shp = shp_2006, 
                                                clustered_data = My_cluster_data_2006, 
                                                variable_interest = "w1&w2")
plot(spatial_poverty2_2006[[4]]) # plot model fit to variogram data
# C)
kriged_poverty2_2006 <- kriging_func(spatial_object = spatial_poverty2_2006, 
                                admin = admin_2006_processed, 
                                variable_interest = "w1&w2")
kriged_poverty2_2006[[2]] # plot predicted variable across map
# D)
poverty2_processed_2006 <- spatially_process_variable_function(kriged = kriged_poverty2_2006[[1]]) # final object: process - gridded object + rasterize

#===================================================================================================#
# Investigating the distribution of the key variables to define meaningful breaks into high and low #

# 1) livestock ownership distribution
# for the final map we are using the 3rd quantile, i.e. the map shows were the 25% clusters with most households with pigs are. 
livestock_ownership_2006 <- distribution_variables_definecutoff_func(spatial_variable = pig_pop_density_2006[[3]],
                                                                cutoff_value = 0.19)
plot(livestock_ownership_2006[[1]]) # histogram
livestock_ownership_2006[[2]] # distribution properties
plot(livestock_ownership_2006[[3]]) # plot high vs low
plot(admin_2006_processed,add = T)

# 2) pig population based on DHS computation
pig_populationDHS_2006 <- distribution_variables_definecutoff_func(spatial_variable = pig_pop_density_2006[[1]],
                                                              cutoff_value = 129.07)
plot(pig_populationDHS_2006[[1]])
pig_populationDHS_2006[[2]]
plot(pig_populationDHS_2006[[3]])
plot(admin_2006_processed,add = T)

# 3) based on FAO Robinson layer (cannot view historgram etc so use diff function)
pig_populationFAO_2006 <- distribution_variables_definecutoff_func2(spatial_variable = pig_pop,
                                                               cutoff_value = 1)
plot(pig_populationFAO[[1]]) # cannot review histogram & summarise distribution
plot(admin_2006_processed,add = T)

# 4) pig extensive rural system 
pig_roam_2006 <- distribution_variables_definecutoff_func(spatial_variable = pig_roam_processed_2006,
                                                     cutoff_value = 0.07703) # cut-off indicated as 0.153 value (not 3rd quartile)
plot(pig_roam_2006[[1]])
pig_roam_2006[[2]]
plot(pig_roam_2006[[3]])
plot(admin_2006_processed,add = T)

# 5) poverty the 40% poorest
poverty_lowest40_2006 <- distribution_variables_definecutoff_func(spatial_variable = poverty2_processed_2006,
                                                             cutoff_value = 0.728294) # 3rd quartile cut-off value (0.666 in original)
plot(poverty_lowest40_2006[[1]])
poverty_lowest40_2006[[2]]
plot(poverty_lowest40_2006[[3]])
plot(admin_2006_processed,add = T)

# 6) poverty the 20% poorest
poverty_lowest20_2006 <- distribution_variables_definecutoff_func(spatial_variable = poverty_processed_2006,
                                                             cutoff_value = 0.48963) # 3rd quartile cut-off value (0.438 used in original)
plot(poverty_lowest20_2006[[1]])
poverty_lowest20_2006[[2]]
plot(poverty_lowest20_2006[[3]])
plot(admin_2006_processed,add = T)

# 7) uncovered sanitation
sanitation_2006 <- distribution_variables_definecutoff_func(spatial_variable = sanitation2_processed_2006,
                                                       cutoff_value = 0.44090) # 3rd quartile cut-off value (0.408 used in original)
plot(sanitation_2006[[1]])
sanitation_2006[[2]]
plot(sanitation_2006[[3]])
plot(admin_2006_processed,add = T)

#===================================================================================================#
#                            Overlays of risk factors (3 risk factors)                             #

# NOTE: we only plot the final overlay risk map (#5) here - please use/adapt the code from 2011 risk maps if other
#       overlays are required

# FINAL RISK MAP; proportion of HH with bad sanitation (1), pig FAO (Robinson) distribution (1.01) and poorest 40% (1.1)  
overlay5_2006 <- plotting_overlays_func(risk_factor1 = sanitation_2006[[3]], risk_factor2 = pig_populationFAO_2006[[1]], 
                                        risk_factor3 = poverty_lowest40_2006[[3]],
                                   admin_processed = admin_2006_processed, admin = admin_2006, year = "2006") 
overlay5_2006[[2]]

#====================================================================================================================#
#                                       2016 PCC risk mapping                                                        #
#====================================================================================================================#
#===========================#
# Load data & mapping files #
DHS_2016_data <- read.spss('~/Uganda_porcine_cysticercosis_risk_mapping/data/DHS 2016/UGHR7BFL.SAV', to.data.frame = T,use.value.labels = FALSE)
geo_2016 <- readShapePoints("~/Uganda_porcine_cysticercosis_risk_mapping/data/DHS 2016/UGGE7AFL.shp") # use 2006 boundaries for 2001 (unable to locate 2001 boundaries)
shp_2016 <- readShapePoly('~/Uganda_porcine_cysticercosis_risk_mapping/data/UGA_2012_60/sdr_subnational_boundaries.shp') # use 2011 boundries ??
admin_2016 <- readShapePoly('~/Uganda_porcine_cysticercosis_risk_mapping/data/admin-gdam/UGA_adm1.shp') # source: https://gadm.org/download_country_v3.html
#districts_2016 <- readShapePoly("~/Uganda_porcine_cysticercosis_risk_mapping/data/District boundaries 2006/Stanford -2010/vg894mz3698.shp") # need this for 2006
human_pop_2016 <- raster('~/Uganda_porcine_cysticercosis_risk_mapping/data/DHS 2016/uga_pd_2016_1km.tif') # human pop data - from world pop 2015 unadjusted; loads and create raster layer object

# pig_pop <- raster('ppop.tif') # pig pop data - robinson 2006 (pig densty map); loads and create raster layer object

#======================================#
#==== Data cleaning & processing ======#

admin_2016_processed <- processing_admin_data_func(admin = admin_2016) # process admin (district) spatial file
shp_2016_processed <- processing_subnational_shp_func(shp = shp_2016) # process subnational boundaries shp file

# data extraction & cleaning for 2001 DHS data variables of interest #
My_cluster_data_2016 <- DHS_extraction_tocluster_func(data = DHS_2016_data, geo = geo_2016, year = 2016) # 

#==========================================================================#
#========  Spatial (spatial autocorrelation & kriging analysis ============#

#============================#
# 1) households with pigs    #

# proportion of HH with pigs by cluster (0-1)
hh_pigs_cluster_plot_2016 <- plot_basic_distribution_variables_func(admin_processed = admin_2016_processed,
                                                                    clustered_data = My_cluster_data_2016, 
                                                                    variable_interest = "pg_h")
hh_pigs_cluster_plot_2016

# plot semi-variogram (spatial autocorrelation) of variable of interest (model vs data)
spatial_pigHH_2016 <- process_plot_variogram_func(shp = shp_2016, 
                                                  clustered_data = My_cluster_data_2016, 
                                                  variable_interest = "pg_h")

plot(spatial_pigHH_2016[[4]]) # plot model fit to variogram data

# spatial interopolation using ordinary kriging (with the underlying variogram) - this may take a few minutes
kriged_pigHH_2016 <- kriging_func(spatial_object = spatial_pigHH_2016, 
                                  admin = admin_2016_processed, 
                                  variable_interest = "pg_h")

kriged_pigHH_2016[[2]] # plot predicted variable across map

# specific for pig population : look at flock sizes (variable = pg_hs)
# this uses Inverse Distance Weighing (deterministic, not a geospatial approach) for spatial interopolation
pig_flocksize_interpolated_2016 <- pig_flock_size_HH_interpolated_func(admin_processed = admin_2016_processed, 
                                                                       clustered_data = My_cluster_data_2016, 
                                                                       spatial_object = spatial_pigHH_2016,
                                                                       variable_interest = "pg_hs")

pig_flocksize_interpolated_2016[[2]] # plot avg. pig herd size spatial datapoint
pig_flocksize_interpolated_2016[[3]] # plot predicted pig flock size 

# Compute and plot pig density (approach #1) : based on human population (from world pop, closest year) x percent of household keeping pigs x flock size
pig_pop_density_2016 <- compute_pig_denisty_func(pop = human_pop_2016, kriged = kriged_pigHH_2016[[1]], 
                                                 shp = shp_2016, idw = pig_flocksize_interpolated_2016[[1]])

plot(pig_pop_density_2016[[1]]) # pig density (shows minimal variation)
plot(pig_pop_density_2016[[2]]) # pig density (> 50 or 100 values excluded re high density urban areas to show more variation)

# Next, we compare it with the FAO Robinson pig distribution map - approach #2
pig_pop_density_FAOmap_2016 <- plot_pig_density_map_func(ppop = pig_pop, admin_processed = admin_2016_processed)
pig_pop_density_FAOmap_2016

#==========================================================================#
# 2) presence of extensive (pig roaming = high potential for transmission) #

# A)
free_range_pigs_cluster_plot_2016 <- plot_basic_distribution_variables_func(admin_processed = admin_2016_processed,
                                                                            clustered_data = My_cluster_data_2016, 
                                                                            variable_interest = "pg_2")
free_range_pigs_cluster_plot_2016
# B)
spatial_pigroam_2016 <- process_plot_variogram_func(shp = shp_2016, 
                                                    clustered_data = My_cluster_data_2016, 
                                                    variable_interest = "pg_2")
plot(spatial_pigroam_2016[[4]]) # plot model fit to variogram data
# C)
kriged_pigroam_2016 <- kriging_func(spatial_object = spatial_pigroam_2016, 
                                    admin = admin_2016_processed, 
                                    variable_interest = "pg_2")
kriged_pigroam_2016[[2]] # plot predicted variable across map
# D)
pig_roam_processed_2016 <- spatially_process_variable_function(kriged = kriged_pigroam_2016[[1]]) # final object: process - gridded object + rasterize

#======================================================================#
#  Sanitation - percentage of household with low sanitation facilities #
#  The first map includes household who have no facilities at all      #

# A)
sanitation_cluster_plot_2016 <- plot_basic_distribution_variables_func(admin_processed = admin_2016_processed,
                                                                       clustered_data = My_cluster_data_2016, 
                                                                       variable_interest = "wc2")
sanitation_cluster_plot_2016
# B)
spatial_sanitation_2016 <- process_plot_variogram_func(shp = shp_2016, 
                                                       clustered_data = My_cluster_data_2016, 
                                                       variable_interest = "wc2")
plot(spatial_sanitation_2016[[4]]) # plot model fit to variogram data
# C)
kriged_sanitation_2016 <- kriging_func(spatial_object = spatial_sanitation_2016, 
                                       admin = admin_2016_processed, 
                                       variable_interest = "wc2")
kriged_sanitation_2016[[2]] # plot predicted variable across map
# D)
sanitation_processed_2016 <- spatially_process_variable_function(kriged = kriged_sanitation_2016[[1]]) # final object: process - gridded object + rasterize

#======================================================================#
#  Sanitation - percentage of household with low sanitation facilities #
#  includes no sanitation or uncovered facilities                      #

# A)
sanitation2_cluster_plot_2016 <- plot_basic_distribution_variables_func(admin_processed = admin_2016_processed,
                                                                        clustered_data = My_cluster_data_2016, 
                                                                        variable_interest = "wc1")
sanitation2_cluster_plot_2016
# B)
spatial_sanitation2_2016 <- process_plot_variogram_func(shp = shp_2016, 
                                                        clustered_data = My_cluster_data_2016, 
                                                        variable_interest = "wc1")
plot(spatial_sanitation2_2016[[4]]) # plot model fit to variogram data
# C)
kriged_sanitation2_2016 <- kriging_func(spatial_object = spatial_sanitation2_2016, 
                                        admin = admin_2016_processed, 
                                        variable_interest = "wc1")
kriged_sanitation2_2016[[2]] # plot predicted variable across map
# D)
sanitation2_processed_2016 <- spatially_process_variable_function(kriged = kriged_sanitation2_2016[[1]]) # final object: process - gridded object + rasterize

#=================================================================================##
##         poverty level                                                          ##
##  based on percentage of household in the poorest 20% (lowest socio-econ quintile)

# A)
poverty_cluster_plot_2016 <- plot_basic_distribution_variables_func(admin_processed = admin_2016_processed,
                                                                    clustered_data = My_cluster_data_2016, 
                                                                    variable_interest = "w1")
poverty_cluster_plot_2016
# B)
spatial_poverty_2016 <- process_plot_variogram_func(shp = shp_2016, 
                                                    clustered_data = My_cluster_data_2016, 
                                                    variable_interest = "w1")
plot(spatial_poverty_2016[[4]]) # plot model fit to variogram data
# C)
kriged_poverty_2016 <- kriging_func(spatial_object = spatial_poverty_2016, 
                                    admin = admin_2016_processed, 
                                    variable_interest = "w1")
kriged_poverty_2016[[2]] # plot predicted variable across map
# D)
poverty_processed_2016 <- spatially_process_variable_function(kriged = kriged_poverty_2016[[1]]) # final object: process - gridded object + rasterize

#=================================================================================##
##         poverty level  # 2                                                     ##
##  based on percentage of household in bottom two socio-economic quintiles       ##

# A)
poverty2_cluster_plot_2016 <- plot_basic_distribution_variables_func(admin_processed = admin_2016_processed,
                                                                     clustered_data = My_cluster_data_2016, 
                                                                     variable_interest = "w1&w2")
poverty2_cluster_plot_2016
# B)
spatial_poverty2_2016 <- process_plot_variogram_func(shp = shp_2016, 
                                                     clustered_data = My_cluster_data_2016, 
                                                     variable_interest = "w1&w2")
plot(spatial_poverty2_2016[[4]]) # plot model fit to variogram data
# C)
kriged_poverty2_2016 <- kriging_func(spatial_object = spatial_poverty2_2016, 
                                     admin = admin_2016_processed, 
                                     variable_interest = "w1&w2")
kriged_poverty2_2016[[2]] # plot predicted variable across map
# D)
poverty2_processed_2016 <- spatially_process_variable_function(kriged = kriged_poverty2_2016[[1]]) # final object: process - gridded object + rasterize

#===================================================================================================#
# Investigating the distribution of the key variables to define meaningful breaks into high and low #

# 1) livestock ownership distribution
# for the final map we are using the 3rd quantile, i.e. the map shows were the 25% clusters with most households with pigs are. 
livestock_ownership_2016 <- distribution_variables_definecutoff_func(spatial_variable = pig_pop_density_2016[[3]],
                                                                     cutoff_value = 0.34)
plot(livestock_ownership_2016[[1]]) # histogram
livestock_ownership_2016[[2]] # distribution properties
plot(livestock_ownership_2016[[3]]) # plot high vs low
plot(admin_2016_processed,add = T)

# 2) pig population based on DHS computation
pig_populationDHS_2016 <- distribution_variables_definecutoff_func(spatial_variable = pig_pop_density_2016[[1]],
                                                                   cutoff_value = 1671.2)
plot(pig_populationDHS_2016[[1]])
pig_populationDHS_2016[[2]]
plot(pig_populationDHS_2016[[3]])
plot(admin_2016_processed,add = T)

# 3) based on FAO Robinson layer (cannot view historgram etc so use diff function)
pig_populationFAO_2016 <- distribution_variables_definecutoff_func2(spatial_variable = pig_pop,
                                                                    cutoff_value = 1)
plot(pig_populationFAO[[1]]) # cannot review histogram & summarise distribution
plot(admin_2016_processed,add = T)

# 4) pig extensive rural system 
pig_roam_2016 <- distribution_variables_definecutoff_func(spatial_variable = pig_roam_processed_2016,
                                                          cutoff_value = 0.06343) # cut-off indicated as 0.153 value (not 3rd quartile)
plot(pig_roam_2016[[1]])
pig_roam_2016[[2]]
plot(pig_roam_2016[[3]])
plot(admin_2016_processed,add = T)

# 5) poverty the 40% poorest
poverty_lowest40_2016 <- distribution_variables_definecutoff_func(spatial_variable = poverty2_processed_2016,
                                                                  cutoff_value = 0.6747814) # 3rd quartile cut-off value (0.666 in original)
plot(poverty_lowest40_2016[[1]])
poverty_lowest40_2016[[2]]
plot(poverty_lowest40_2016[[3]])
plot(admin_2016_processed,add = T)

# 6) poverty the 20% poorest
poverty_lowest20_2016 <- distribution_variables_definecutoff_func(spatial_variable = poverty_processed_2016,
                                                                  cutoff_value = 0.50827) # 3rd quartile cut-off value (0.438 used in original)
plot(poverty_lowest20_2016[[1]])
poverty_lowest20_2016[[2]]
plot(poverty_lowest20_2016[[3]])
plot(admin_2016_processed,add = T)

# 7) uncovered sanitation
sanitation_2016 <- distribution_variables_definecutoff_func(spatial_variable = sanitation2_processed_2016,
                                                            cutoff_value = 0.8847) # 3rd quartile cut-off value (0.408 used in original)
plot(sanitation_2016[[1]])
sanitation_2016[[2]]
plot(sanitation_2016[[3]])
plot(admin_2016_processed,add = T)

#===================================================================================================#
#                            Overlays of risk factors (3 risk factors)                             #

# NOTE: we only plot the final overlay risk map (#5) here - please use/adapt the code from 2011 risk maps if other
#       overlays are required

# FINAL RISK MAP; proportion of HH with bad sanitation (1), pig FAO (Robinson) distribution (1.01) and poorest 40% (1.1)  
overlay5_2016 <- plotting_overlays_func(risk_factor1 = sanitation_2016[[3]], risk_factor2 = pig_populationFAO_2016[[1]], 
                                        risk_factor3 = poverty_lowest40_2016[[3]],
                                        admin_processed = admin_2016_processed, admin = admin_2016, year = "2016") 
overlay5_2016[[2]]


#================================================================================================#
#             Plotting risk maps together & animation through time                               #
#================================================================================================#

# 1) plot all risk factor maps (2001 - 2016) #
Riskmap_2001to2016 <- plot_all_risk_maps_func(riskmapdata1 = overlay5_2001[[3]], 
                        riskmapdata2 = overlay5_2006[[3]], 
                        riskmapdata3 = overlay5_2011[[3]], 
                        riskmapdata4 = overlay5_2016[[3]],
                        admin_processed = admin_2011_processed)
Riskmap_2001to2016[[2]]
#ggsave("risk_factormaps.tiff")
savePlot_riskmaps(Plot = Riskmap_2001to2016[[2]]) # save a pdf in working dir. 

# 2) Animate - takes a minute #
Riskmap_animated_2001to2016 <- plot_riskmaps_animation_func(riskmapdata1 = overlay5_2001[[3]], 
                                                   riskmapdata2 = overlay5_2006[[3]], 
                                                   riskmapdata3 = overlay5_2011[[3]], 
                                                   riskmapdata4 = overlay5_2016[[3]],
                                                   admin_processed = admin_2011_processed)
Riskmap_animated_2001to2016 # will save a animated gif in working dir. 

