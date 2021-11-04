#===================================================================================================================#
#=====================================    Data manipulation/ processing  extraction ================================#


# ===================================== #
# process admin (district) spatial file #
processing_admin_data_func <- function(admin) {
  
  admin$NAME_1 <- as.character(admin$NAME_1) # make character variable for district names
  
  admin$label <- NA 
  
  admin$label[(admin$NAME_1 == 'Lira'|admin$NAME_1 == 'Masaka'|admin$NAME_1 == 'Mukono'|admin$NAME_1 == 'Kamuli'|
                 admin$NAME_1 =='Hoima'|admin$NAME_1 =='Moyo'|admin$NAME_1 == 'Kumi'|admin$NAME_1 == 'Apac'|
                 admin$NAME_1 =='Kaberamaido'|admin$NAME_1 == 'Soroti'|admin$NAME_1 == 'Kayunga')] <- 
    admin$NAME_1[(admin$NAME_1 == 'Lira'|admin$NAME_1 == 'Masaka'|admin$NAME_1 == 'Mukono'|admin$NAME_1 == 'Kamuli'|
                    admin$NAME_1 == 'Hoima'|admin$NAME_1 == 'Moyo'|admin$NAME_1 == 'Kumi'|admin$NAME_1 == 'Apac'|
                    admin$NAME_1 =='Kaberamaido'|admin$NAME_1 == 'Soroti'|admin$NAME_1 == 'Kayunga')]
  
  admin2 <- fortify(admin)
  
  return(admin)
}

# ======================================= #
# process subnational boundaries shp file #
processing_subnational_shp_func <- function(shp) {
  
  shp2 <- fortify(shp)
  
  return(shp2)
}

#========================================================#
# Function to extract (below) from DHS standard surveys  #

DHS_extraction_tocluster_func <- function(data, geo, year) {
  
  # Variables of interest:
  # 1. toilet that a household use
  # 2. number of livestock owned by the households
  # 3. a couple of household charateristics
  
  # recode variables of interest from DHS #
  
  # 2001 DHS extraction #
  if(year == 2001) {
    data$wc2 <- ifelse(data$HV205 == '31', 1, 0) # if no toilet facility/bush/field = 1 
    data$w1 <- ifelse(data$SH052 == '1', 1, 0) # (HV270: wealth quintile); if = 1 (poorest or lowest 20%) = 1
    data$w2 <- ifelse(data$SH052 == '2', 1, 0) # (HV270: wealth quintile); if = 2 (poorer or next to lowest 20%) = 1
    data$w3 <- ifelse(data$SH052 == '3', 1, 0) # (HV270: wealth quintile); if = 3 (middle 20%) = 1
    data$w4 <- ifelse(data$SH052 == '4', 1, 0) # (HV270: wealth quintile); if = 4 (richer or next to highest 20%) = 1
    data$w5 <- ifelse(data$SH052 == '5', 1, 0) # (HV270: wealth quintile); if = 5 (richest or highest 20%) = 1
  } 
  
  # 2006 DHS extraction #
  if(year == 2006) {
  data$Pg_num <- ifelse(data$HV246G > 95.5, NA, data$HV246G) # HG246G = Cs owns pigs: replace entries with >95 pigs with NA
  data$Pg_h <- ifelse (data$Pg_num > 0, 1, 0) # if pig number > 0, specify 1 for datapoint (household with pigs)
  data$Pg_2 <- ifelse (data$Pg_num == 1 | data$Pg_num == 2 & data$HV025 == 'Rural', 1, 0) # (HV025 = rural v urban) if pig number is 1 or 2 & classified as rural = 1
  data$wc1 <- ifelse(data$HV205 == '31' | data$HV205 == '25' | data$HV205 == '24', 1, 0) # HV205 type of toilet facility; if either of (31 = no facility/bush/field; 25 = uncovered pit latrine with slab; 24 = uncovered pit latrine no slab) = 1
  data$wc2 <- ifelse(data$HV205 == '31', 1, 0) # if no toilet facility/bush/field = 1 
  data$w1 <- ifelse(data$HV270 == '1', 1, 0) # (HV270: wealth quintile); if = 1 (poorest or lowest 20%) = 1
  data$w2 <- ifelse(data$HV270 == '2', 1, 0) # (HV270: wealth quintile); if = 2 (poorer or next to lowest 20%) = 1
  data$w3 <- ifelse(data$HV270 == '3', 1, 0) # (HV270: wealth quintile); if = 3 (middle 20%) = 1
  data$w4 <- ifelse(data$HV270 == '4', 1, 0) # (HV270: wealth quintile); if = 4 (richer or next to highest 20%) = 1
  data$w5 <- ifelse(data$HV270 == '5', 1, 0) # (HV270: wealth quintile); if = 5 (richest or highest 20%) = 1
  }
  
  
  # 2011 DHS extraction #
  if(year == 2011) {
    data$Pg_num <- ifelse(data$HV246G > 95.5, NA, data$HV246G) # HG246G = Cs owns pigs: replace entries with >95 pigs with NA
    data$Pg_h <- ifelse (data$Pg_num > 0, 1, 0) # if pig number > 0, specify 1 for datapoint (household with pigs)
    data$Pg_2 <- ifelse (data$Pg_num == 1 | data$Pg_num == 2 & data$HV025 == 'Rural', 1, 0) # (HV025 = rural v urban) if pig number is 1 or 2 & classified as rural = 1
    data$wc1 <- ifelse(data$HV205 == '31' | data$HV205 == '25' | data$HV205 == '24', 1, 0) # HV205 type of toilet facility; if either of (31 = no facility/bush/field; 25 = uncovered pit latrine with slab; 24 = uncovered pit latrine no slab) = 1
    data$wc2 <- ifelse(data$HV205 == '31', 1, 0) # if no toilet facility/bush/field = 1 
    data$w1 <- ifelse(data$HV270 == '1', 1, 0) # (HV270: wealth quintile); if = 1 (poorest) = 1
    data$w2 <- ifelse(data$HV270 == '2', 1, 0) # (HV270: wealth quintile); if = 2 (poorer) = 1
    data$w3 <- ifelse(data$HV270 == '3', 1, 0) # (HV270: wealth quintile); if = 3 (middle) = 1
    data$w4 <- ifelse(data$HV270 == '4', 1, 0) # (HV270: wealth quintile); if = 4 (richer) = 1
    data$w5 <- ifelse(data$HV270 == '5', 1, 0) # (HV270: wealth quintile); if = 5 (richest) = 1
  }
  
  # extract DHS geo data for each variable of interest #
  myClu = geo@data 
  myClu <- myClu[which(!myClu$LONGNUM == 0),] # subset (exclude Long & Lat = 0)
  
  myClu = myClu[,c("DHSCC","DHSCLUST","DHSREGCO","DHSREGNA","URBAN_RURA", "ADM1NAME", "LONGNUM","LATNUM")] # subset on these col's
  names(myClu) = c("Cc","HV001","REGCODE", "HV024","IsRural", "AdName","x","y") # rename & replace col names
  
  
  # Estimate aggregate/average number fo (variable of interest) within each DHS cluster #
  
  # pig variables #
  if(year == 2011 || year == 2006) {
    temp <- aggregate(x = data$Pg_num, by = list(data$HV001), 'mean', na.rm = T) # aggregate/ average (mean) of pig number by cluster number (HV001)
    myPosVec = match(myClu$HV001, temp$Group.1 ) # match those clusters in temp (with an average pig number) to cluster in total list
    myClu$Pg_num = temp$x[myPosVec] # add average pig number for each cluster to full table
  
    temp <- aggregate(x=ifelse(data$Pg_num == 0,NA,data$Pg_num), by = list(data$HV001), 'mean', na.rm = T) # aggregate further (??) & create new vector replacing 0 with NA for pig number (to replace 0s in average pig number by cluster)
    myPosVec = match(myClu$HV001,temp$Group.1 )
    myClu$Pg_hs = temp$x[myPosVec] # average number of pigs per HH per cluster (?????)
  
    temp <- aggregate(x = data$Pg_h, by = list(data$HV001), 'mean', na.rm = T) # proportion of HH with pigs (by HH in cluster) - aggregate/average number of HH with pigs by cluster (calc mean)
    myPosVec = match(myClu$HV001,temp$Group.1 )
    myClu$Pg_h = temp$x[myPosVec] # average number of HH with pigs per cluster
    
    temp<- aggregate(x = data$Pg_2, by = list(data$HV001), 'mean', na.rm = T) # aggregate/average number of HH (with 1 or 2 pigs) classified as rural by cluster (calc mean)
    myPosVec = match(myClu$HV001,temp$Group.1 )
    myClu$Pg_2 = temp$x[myPosVec] # average number of HH with 1 or 2 pigs in rural area (per cluster)
  }
  

  # sanitation variables #
  if(year == 2011 || year == 2006) {
    temp<- aggregate(x = data$wc1, by = list(data$HV001), 'mean', na.rm = T) # aggregate/avg number reporting no facility/bush field or uncovered latrine w/ or /wo slab per cluster (calc mean)
    myPosVec = match(myClu$HV001,temp$Group.1 )
    myClu$wc1 = temp$x[myPosVec] # avg number reporting no facility/bush field or uncovered latrine w/ or /wo slab (per cluster)
  }
  
  temp <- aggregate(x = data$wc2, by = list(data$HV001), 'mean', na.rm = T) # aggregate/avg number reporting no facility/bush field per cluster (calc mean)
  myPosVec = match(myClu$HV001,temp$Group.1 )
  myClu$wc2 = temp$x[myPosVec] # avg number reporting no facility/bush field (per cluster)
  
  
 # poverty variables #
  temp <- aggregate(x = data$w1, by = list(data$HV001), 'mean', na.rm = T) # aggregate/avg number classified as poorest (1) quintile per cluster (calc mean)
  myPosVec = match(myClu$HV001,temp$Group.1 )
  myClu$w1 = temp$x[myPosVec] # avg number classified as poorest (per cluster)
  
  temp <- aggregate(x = data$w2, by = list(data$HV001), 'mean', na.rm = T) # aggregate/avg number classified as poorer (2) quintile per cluster (calc mean)
  myPosVec = match(myClu$HV001,temp$Group.1 )
  myClu$w2 = temp$x[myPosVec]  # avg number classified as poorer (per cluster) 
  
  temp <- aggregate(x = data$w3, by = list(data$HV001), 'mean', na.rm = T) # aggregate/avg number classified as middle (3) quintile per cluster (calc mean)
  myPosVec = match(myClu$HV001,temp$Group.1 )
  myClu$w3 = temp$x[myPosVec] # avg number classified as middle (per cluster) 
  
  temp <- aggregate(x = data$w4, by = list(data$HV001), 'mean', na.rm = T) # aggregate/avg number classified as richer (4) quintile per cluster (calc mean)
  myPosVec = match(myClu$HV001,temp$Group.1 )
  myClu$w4 = temp$x[myPosVec] # avg number classified as richer (per cluster) 
  
  temp <- aggregate(x = data$w5, by = list(data$HV001), 'mean', na.rm = T) # aggregate/avg number classified as richest (5) quintile per cluster (calc mean)
  myPosVec = match(myClu$HV001,temp$Group.1 )
  myClu$w5 = temp$x[myPosVec] # avg number classified as richest (per cluster)   
  
  return(myClu)
  
}


# ============================================================================= #
# function to explore distirbution of key variables to define low/high cut-offs #

distribution_variables_definecutoff_func <- function(spatial_variable, cutoff_value){
  
  histogram <- hist(spatial_variable@data@values) # histogram of HH with pigs
  
  dist_summary <- summary(spatial_variable@data@values)
  
  spatial_variable_cutoff <- spatial_variable
  
  spatial_variable_cutoff[spatial_variable > cutoff_value] <- 1 #values over 23% = 1
  
  spatial_variable_cutoff[spatial_variable_cutoff < cutoff_value] <- 0 # values less than 23% = 0
  
  #dist_plot <- plot(spatial_variable_cutoff)
  
  #dist_plot_admin <- dist_plot + plot(admin,add=T)

return(list(histogram, dist_summary, spatial_variable_cutoff))

}

# =================================================================================== #
# same as function above but not assessing histrogram/distribution for FAO Robinson pig pop layer (as this does not work)
distribution_variables_definecutoff_func2 <- function(spatial_variable, cutoff_value){
  
  spatial_variable_cutoff <- spatial_variable
  
  spatial_variable_cutoff[spatial_variable > cutoff_value] <- 1 #values over 23% = 1
  
  spatial_variable_cutoff[spatial_variable_cutoff < cutoff_value] <- 0 # values less than 23% = 0

  return(list(spatial_variable_cutoff))
  
}