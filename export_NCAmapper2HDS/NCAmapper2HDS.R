# export NCAmapper output to HDS input
# used to parameterize HDS.
# load needed libraries
library(sf)
library(terra)
##############
if (rstudioapi::isAvailable()) {
  if (require('rstudioapi') != TRUE) {
    install.packages('rstudioapi')
    
  } else{
    library(rstudioapi) # load it
  }
  wdir <- dirname(getActiveDocumentContext()$path)
} else{
  wdir <- getwd()
}
setwd(wdir)
###########################
# read shapefile of the subbasin to calculate regional stats
subbasins <- st_read('/Users/mohamed/scratch/spatial_data/modified_shapefiles/Modified_SMMcat.shp')
subbasins$area <- st_area(subbasins)
# read depression depth raster
depr_depth_ras <- rast('/Users/mohamed/scratch/Milk_river_basin/NCAmapper_temp/depression_depth_non_zeros.tif')
#read depression watershed raster
depr_watershed_ras <- rast('/Users/mohamed/scratch/Milk_river_basin/NCAmapper_temp/depression_watershed.tif')
# column name with subbasin ID
subID_colname <- 'hru_nhm'
# remove depressions from contributing area (check)
depr_depth_ras[is.na(depr_watershed_ras)] <- NA

# extract information from raster for each subbasin
# sub_basin	depression_area_m2	depression_volume_m3	total_catchment_m2
# 1	84035800	18594986.11	385052700
HDS_parameters <- data.frame(matrix(NA, nrow = nrow(subbasins), ncol = 8))
colnames(HDS_parameters) <- c('subid', 'depressionArea_m2', 'depressionVol_m3', 'depressionCatchment_m2', 'landArea_m2','depressionDepth_m',
                              'depressionAreaFrac',	'deprCatchAreaFrac')


my_summary <- function(x, na.rm) c(mean = mean(x, na.rm=na.rm), 
                                   min = min(x, na.rm=na.rm), 
                                   max = max(x, na.rm=na.rm), 
                                   sum = sum(x, na.rm=na.rm),
                                   count = length(na.omit(x)))


# extract information for all subbasins
depth_summary <- as.data.frame(terra::extract(depr_depth_ras, subbasins, fun=my_summary, na.rm=T))
colnames(depth_summary) <- c('ID', 'mean', 'min', 'max', 'sum', 'count')

#depression watershed area
watershed_summary <- terra::extract(depr_watershed_ras, subbasins, fun=\(x) length(na.omit(x)))
# get grid area
grid_area <- res(depr_depth_ras)[1]*res(depr_depth_ras)[2]


HDS_parameters$subid <- subbasins[[subID_colname]]
HDS_parameters$depressionArea_m2 <- depth_summary$count*grid_area
HDS_parameters$depressionVol_m3 <- depth_summary$sum*grid_area
HDS_parameters$depressionCatchment_m2 <- watershed_summary$depression_watershed*grid_area
HDS_parameters$landArea_m2 <- as.numeric(subbasins$area) - HDS_parameters$depressionArea_m2
HDS_parameters$depressionDepth_m <- HDS_parameters$depressionVol_m3/HDS_parameters$depressionArea_m2
HDS_parameters$depressionAreaFrac <- as.numeric(HDS_parameters$depressionArea_m2/subbasins$area)
HDS_parameters$deprCatchAreaFrac <- as.numeric((HDS_parameters$depressionCatchment_m2 - HDS_parameters$depressionArea_m2)
                                               /HDS_parameters$landArea_m2)#subbasins$area)

# some checks
HDS_parameters$deprCatchAreaFrac[HDS_parameters$deprCatchAreaFrac>1] = 1
HDS_parameters[is.na(HDS_parameters)]=0


# write output to file
write.table(HDS_parameters, 'Milk_HDS_parameters.csv', sep = ',', row.names = F)
