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

# 1) households with pigs

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

# kriging 
kriged_pigHH <- kriging_func(spatial_object = spatial_pigHH, 
                             admin = admin_2011_processed, 
                             variable_interest = "pg_h")

kriged_pigHH[[2]] # plot predicted variable across map

# specific for pig population : look at flock sizes (variable = pg_hs)
pig_flocksize_kriged <- pig_flock_size_HH_kriged_func(admin_processed = admin_2011_processed, 
                                                      clustered_data = My_cluter_data_2011, 
                                                      spatial_object = spatial_pigHH,
                                                      variable_interest = "pg_hs")

pig_flocksize_kriged[[2]] # plot avg. pig herd size spatial datapoint
pig_flocksize_kriged[[3]] # plot predicted pig flock size 
