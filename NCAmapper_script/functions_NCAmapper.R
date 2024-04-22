# Functions used by the NCAmapper model (main_NCAmapper)
############################
###########################
init_NCAmapper_config <- function(inp_dir){
  
  # Print the model credits to the screen
  cat('Non-Contributing Area mapping tool \n')
  cat('NCAmapper (v 1.10) \n')
  cat('© Mohamed Ismaiel Ahmed 2024 @ UCalgary \n')
  cat('Initializing NCAmapper \n')
  # set seed number
  set.seed(123456788)
  # read config file and save to list
  
  # Read lines from the file
  lines <- readLines(paste0(inp_dir,'/NCAmapper_config.ini'))
  
  # Remove comments at the end of each line
  lines <- gsub("#.*$", "", lines)
  
  # Filter out lines starting with "#"
  filtered_lines <- lines[!grepl("^#", lines)]
  
  # Split each line into elements
  config_data <- lapply(filtered_lines, function(line) unlist(strsplit(line, "\t")))
  
  # Convert the list to a data frame
  config_data <- as.data.frame(do.call(rbind, config_data), stringsAsFactors = FALSE)
  # Remove extra whitespace from values
  config_data[] <- lapply(config_data, trimws)
  # Create a named list for easy lookup with appropriate types
  key_value_map <- setNames(
    lapply(config_data$V2, function(value) type.convert(value, as.is = TRUE)),
    config_data$V1
  )
  
  # add the inp_dir to the list
  key_value_map$inp_dir <- inp_dir
  
  
  # install missing packages
  # list of dependecies (packages) for the NCAmapper modelexcept whitebox as it needs special treatment
  packages_dependencies <- c('data.table', 'raster', 'sf', 'terra', 'viridis')
  # missing packages (needs to be installed)
  new.packages <- packages_dependencies[!(packages_dependencies %in% installed.packages()[,"Package"])]
  # install the missing packages
  if(length(new.packages)) install.packages(new.packages)
  
  
  #load the required libraries
  if (!require(whitebox)){
    install.packages("whitebox")
    # install binaries for whitebox
    if(Sys.info()['sysname'] == "Linux"){
      install_whitebox(platform="linux_musl") #this is for unix
    } else{
      install_whitebox() #for windows/mac
    }
  } 
  # load needed libraries and set the number of cores for processing
  require('data.table') #faster processing of dataframes
  require('raster') #working with raster objects
  require('sf') #working with shapefiles
  require('terra') #faster cropping and raster operations
  # libraries required for plotting
  # require("ggnewscale")# use multiple scales
  # require('ggplot2')
  require('viridis')
  # require("stars")
  # require(profvis) #benchmarking and measuring code time
  # require(tictoc)
  # Set ncores used by different packages
  wbt_max_procs(max_procs = key_value_map$ncores)
  setDTthreads(threads = key_value_map$ncores)
  
  # create output directories
  dir.create(paste0(key_value_map$inp_dir, '/', key_value_map$temp_dir))
  dir.create(paste0(key_value_map$inp_dir, '/', key_value_map$out_dir))
  
  return(key_value_map)
  
}
##########################
preprocessing_wbt <- function(config_file){
  # free up unused memory
  gc()
  #preprossing using whitebox
  # define function arguments
  
  DEM_name <- paste0(config_file$inp_dir, '/', config_file$DEM_name)
  rvr_data <-  paste0(config_file$inp_dir, '/', config_file$rvr_data)
  temp_dir <- paste0(config_file$inp_dir, '/', config_file$temp_dir)
  dep_rvr_areafrac <- config_file$dep_rvr_areafrac
  dep_dpth_threshold <- config_file$dep_dpth_threshold
  buffer_dist <- config_file$buffer_dist
  
  cat('---------------------------- \n')
  cat('Preprossing the DEM','\n')
  cat('---------------------------- \n')
  start_time = Sys.time()
  
  #check if the provided DEM & shp are in projected system
  if(isLonLat(raster(DEM_name)) | st_is_longlat(st_read(rvr_data, quiet = T))){
    stop("STOP! Make sure that your raster and vector files are in projected crs")
    }
  
  
  
  #Note to run everything correctly you need to convert your data to "*.tif" & uncompressed file format in metric units
  ##############
  #run whitebox tools
  #identify depressions using sink function
  wbt_sink(input = DEM_name, wd = temp_dir, output = paste0(temp_dir,'/depression_no.tif'), verbose_mode = F)
  
  #get depressions depth
  wbt_depth_in_sink(dem = DEM_name, wd = temp_dir, output = paste0(temp_dir,'/depression_depth.tif'), verbose_mode = F)
  
  #--------#--------
  #--------#--------
  # remove cells on the main rivers identified from the NHN datasets
  dep_no_r <- raster(paste0(temp_dir,'/depression_no.tif'))
  # ##########
  # #read the river lines
  main_river_poly <- st_read(rvr_data)
  #apply buffer
  if (buffer_dist>0) {
    #apply fixed buffer
    main_river_poly <- st_buffer(x = main_river_poly, dist = buffer_dist)
  }else{ 
    # apply variable buffer
    for (i in 1:nrow(main_river_poly)) {
      main_river_poly[i,] <- st_buffer(x = main_river_poly[i,], dist = (main_river_poly$width_m[i]/2))
    }
  }
  
  # re-project the shp file to the same crs as the raster
  main_river_poly <- st_transform(main_river_poly, crs = st_crs(dep_no_r))
  # slow
  # write the re-projected shapefile
  # if(file.exists(paste0(temp_dir,'/main_river_poly.shp'))){unlink(paste0(temp_dir,'/main_river_poly.*'))}
  # shapefile(as(main_river_poly, "Spatial"), filename = paste0(temp_dir,'/main_river_poly.shp')) # doesn't support overwrite
  # #clip the raster to the polygon using WBT (faster than the sf package)
  # wbt_clip_raster_to_polygon(input = paste0(temp_dir,'/depression_no.tif'), wd = temp_dir, polygons = paste0(temp_dir,'/main_river_poly.shp'),
  #                            output = paste0(temp_dir,'/depression_no_rvr.tif'), verbose_mode = F)
  ## read the cropped raster that contains the depressions within the river
  # dep_no_rvr_r <- raster(paste0(temp_dir,'/depression_no_rvr.tif'))
  # mask raster directly in R using terra package (much faster)
  dep_no_rvr_r <- terra::mask(rast(paste0(temp_dir,'/depression_no.tif')), main_river_poly) #mask(dep_no_r, main_river_poly)
  #write the raster to file
  writeRaster(dep_no_rvr_r, paste0(temp_dir,'/depression_no_rvr.tif'), overwrite=TRUE)
  #convert it back to rasterlayer
  dep_no_rvr_r <- raster(dep_no_rvr_r)
  #get dep_no and number of cells for each dep_no
  dep_rvr <- data.table(dep_no_rvr =getValues(dep_no_rvr_r))
  dep_rvr <- dep_rvr[, .(count = .N), by = dep_no_rvr]
  dep_rvr <- na.omit(dep_rvr)
  
  # get the total number of cells for all dep_no
  dep_no <- data.table(dep_no=getValues(dep_no_r))
  dep_no <- dep_no[, .(tot_count = .N), by = dep_no]
  dep_no <- na.omit(dep_no)
  #left joint
  dep_rvr <- merge(dep_rvr, dep_no, by.x = "dep_no_rvr", by.y = "dep_no", all.x = TRUE)
  dep_rvr$frac_in_rvr <- dep_rvr$count/dep_rvr$tot_count
  # any depression with 25% of it in the river is to be removed
  dep_no_in_rvr_toremove <- dep_rvr$dep_no_rvr[dep_rvr$frac_in_rvr> dep_rvr_areafrac]
  #--------#--------
  #--------#--------
  
  
  ############
  # Remove depressions with depth <= threshold to reduce the runtime
  dep_no_r <- raster(paste0(temp_dir,'/depression_no.tif'))
  dep_dpth_r <- raster(paste0(temp_dir,'/depression_depth.tif'))
  #data.table for faster indexing
  dep_info_dt <- data.table(dep_no=getValues(dep_no_r), dep_dpth=getValues(dep_dpth_r))
  tot_ndep <- length(na.omit(unique(dep_info_dt$dep_no)))
  
  # Find unique depression IDs
  unique_dep_ids <- unique(na.omit(dep_info_dt$dep_no))
  
  # Calculate maximum depths per depression
  max_depths <- dep_info_dt[, .(max_value = max(dep_dpth)), by = dep_no]
  
  # Find depression IDs with max depth <= dep_dpth_threshold
  zero_depth_ids <- na.omit(max_depths$dep_no[max_depths$max_value <= dep_dpth_threshold])
  # Identify the cells associated with these depression IDs
  cells_to_remove <- dep_info_dt$dep_no %in% zero_depth_ids
  
  # Set the corresponding cells to NA
  dep_info_dt$dep_no[cells_to_remove] <- NA
  dep_info_dt$dep_dpth[cells_to_remove] <- NA
  
  # remove depressions within the river
  dep_in_rvr_id <- which(dep_info_dt$dep_no %in% dep_no_in_rvr_toremove)
  # Set the corresponding cells to NA
  dep_info_dt$dep_no[dep_in_rvr_id] <- NA
  dep_info_dt$dep_dpth[dep_in_rvr_id] <- NA
  
  
  #re assign depression no to go from 1:n after removing some depressions
  unique_dep_ids <- sort(na.omit(unique(dep_info_dt$dep_no)))

  # Create a lookup data.table for dep_no mapping
  dep_no_lookup <- data.table(orig_dep_no = unique_dep_ids, new_dep_no = seq_along(unique_dep_ids))
  # Merge the original data.table with the lookup table to update dep_no
  # Map the new dep_no to old ones
  dep_info_dt[dep_no_lookup, dep_no := i.new_dep_no, on = .(dep_no = orig_dep_no)]
  
  values(dep_no_r) <- dep_info_dt$dep_no
  values(dep_dpth_r) <- dep_info_dt$dep_dpth
  
  #write the depression no and depth new rasters
  writeRaster(dep_no_r,paste0(temp_dir,'/depression_no_non_zeros.tif'),overwrite=T,
              options=c("COMPRESS=NONE", paste0("NUM_THREADS=",config_file$ncores)))
  writeRaster(dep_dpth_r,paste0(temp_dir,'/depression_depth_non_zeros.tif'),overwrite=T,
              options=c("COMPRESS=NONE", paste0("NUM_THREADS=",config_file$ncores)))
  ############

  #get flow direction for filled dem (will be used to get contributing area)
  
  #filldepressions
  wbt_fill_depressions(dem = DEM_name, wd = temp_dir, output = paste0(temp_dir,'/filled_dem.tif'),verbose_mode = F)
  # #burn streams into DEM not working 
  # wbt_fill_burn(dem = DEM_name, streams = rvr_data, output = paste0(temp_dir,'/filled_dem_stream_burn.tif'), wd = temp_dir, verbose_mode = F)
  #generate d8 pointer based on the filled DEM
  wbt_d8_pointer(wd = temp_dir, dem = paste0(temp_dir,'/filled_dem.tif'), output = paste0(temp_dir,'/D8_flow_pointer_filled.tif'), verbose_mode = F)
  #get location of max flow accumulation (possible outlet)
  wbt_d8_flow_accumulation(input = paste0(temp_dir,'/filled_dem.tif'), wd = temp_dir, output = paste0(temp_dir,'/flow_acc_fill.tif'),verbose_mode = F)
  #get the watersheds of each depression using 'Watershed' function
  wbt_watershed(wd = temp_dir, d8_pntr = paste0(temp_dir,'/D8_flow_pointer_filled.tif'),pour_pts = paste0(temp_dir,'/depression_no_non_zeros.tif'),output = paste0(temp_dir,'/depression_watershed.tif'),verbose_mode = F)
  #find main stem (main channel). This might be used to remove depressions from the main channel
  # wbt_extract_streams(flow_accum = 'flow_acc_fill.tif', output = 'streams.tif', threshold = 1e4, wd = temp_dir, verbose_mode = F)
  # wbt_find_main_stem(d8_pntr = 'D8_flow_pointer_filled.tif', streams = 'streams.tif', output = 'main_stem.tif', wd = temp_dir)
  end_time = Sys.time()
  cat('Finished pre-processing the inputs in', as.numeric(end_time-start_time,units="secs"), 'sec \n')
}
##########################

get_depressions_summary <- function(config_file){
  
  # free up unused memory
  gc()
  
  cat('------------------------------------------ \n')
  cat('Generating depressions information table','\n')
  cat('------------------------------------------ \n')
  
  # function arguments
  DEM_name <- paste0(config_file$inp_dir, '/', config_file$DEM_name)
  temp_dir <- paste0(config_file$inp_dir, '/', config_file$temp_dir)
  out_dir <- paste0(config_file$inp_dir, '/', config_file$out_dir)
  
  start_time = Sys.time()
  # dir.create(out_dir)
  # read raster and get values
  facc_fill_r <- raster(paste0(temp_dir, '/flow_acc_fill.tif'))
  facc_fill_val <- getValues(facc_fill_r)

  dep_no_r <- raster(paste0(temp_dir,'/depression_no_non_zeros.tif'))
  # dep_dpth_r <- raster('./output/depression_depth_non_zeros.tif')
  filldem_r <- raster(paste0(temp_dir,'/filled_dem.tif'))
  #dep_basin_r <- raster('./output/depression_watershed.tif')
  dep_basin_val <- getValues(raster(paste0(temp_dir,'/depression_watershed.tif')))
  d8_pntr_val <- getValues(raster(paste0(temp_dir,'/D8_flow_pointer_filled.tif')))
  
  dep_no_val <- getValues(dep_no_r)
  dep_dpth_val <- getValues(raster(paste0(temp_dir,'/depression_depth_non_zeros.tif')))
  DEM_val <- getValues(raster(DEM_name))
  full_wl_val <- getValues(filldem_r)
  #get number of active_cells
  active_cells <- length(which(DEM_val>=0))
  #get the number of depressions 
  tot_ndep <- length(na.omit(unique(dep_no_val)))
  cell_size <- res(dep_no_r)
  #create depression summary datatable
  dep_no_dpth_dt <- data.table(dep_no=dep_no_val, dep_dpth=dep_dpth_val, max_wl=full_wl_val,
                               dep_basin=dep_basin_val, facc_fill = facc_fill_val)
  dep_no_val_dt <- data.table(dep_no_val)
  dep_basin_val_dt <- data.table(dep_basin_val)
  
  
  # Calculate max, mean, count, sum, and index of max_value by dep_no
  result <- dep_no_dpth_dt[, .(max_depth_m = max(dep_dpth),
                               max_wl_m =   NA,
                               mean_depth_m = mean(dep_dpth),
                               dep_area_m2 = .N*cell_size[1]*cell_size[2],
                               dep_vol_m3 = sum(dep_dpth)*cell_size[1]*cell_size[2],
                               max_d_idx = .I[which.max(dep_dpth)]
                               ),
                          by = dep_no]
  result$max_wl_m <- dep_no_dpth_dt$max_wl[result$max_d_idx]
  
  #get depression basin area
  #count number of cells in the depression_watershed
  counts <- dep_no_dpth_dt[, .(basin_area_m2 = .N*cell_size[1]*cell_size[2],
                               spill_loc_idx = .I[which.max(facc_fill)]), by = dep_basin]

  #left joint the data.tables
  result <- merge(result, counts, by.x = "dep_no", by.y = "dep_basin", all.x = TRUE)
  result <- na.omit(result)
  #####get the next downstream depression
  ###get neighbouring cells id will be added inside the loop
  spill_loc_d8flow_code <- data.table(dep_no=result$dep_no, spill_loc_idx = result$spill_loc_idx, flow_code=d8_pntr_val[result$spill_loc_idx])
  
  # get ids of neighboring cells (using raster package)
  id_neighbour_dt <- adjacent(facc_fill_r, result$spill_loc_idx, directions = 8)
  # #associate depression no with the data.table
  id_neighbour_dt <- data.table(id_neighbour_dt,dep_no=rep(result$dep_no,8))
  #sort by dep no to put all neighbors together
  id_neighbour_dt <- setorder(id_neighbour_dt, dep_no)
  #get D8 flow direction code for central cell
  id_neighbour_dt$DS_direction <- d8_pntr_val[id_neighbour_dt$from]
  #generate D8 flow direction code for neighbor cell (D8 rule: 64, 32, 16, 1, 2, 4, 128, 8; 'NW', 'W', 'SW', 'NE', 'E', 'SE', 'N', 'S')
  id_neighbour_dt$ngbr_direction <- rep(c(64, 32, 16, 1, 2, 4, 128, 8), length(result$spill_loc_idx))
  # keep rows where DS_direction = ngbr_direction (i.e., centeral cell direction to ds cell)
  id_neighbour_dt <- id_neighbour_dt[DS_direction == ngbr_direction]
  # get the next DS depression based on depression watershed and to value from id_neighbour
  id_neighbour_dt$DS_dep_no <- dep_no_dpth_dt$dep_basin[id_neighbour_dt$to]
  # add this to the result dt
  result$next_downstream_dep <- id_neighbour_dt$DS_dep_no
  # # ##**##**##**##**##**
  
  # ###### if decided to switch to terra (check v05)
  # ###get neighbouring cells id will be added inside the loop
  # spill_loc_d8flow_code <- data.table(dep_no=result$dep_no, spill_loc_idx = result$spill_loc_idx, flow_code=d8_pntr_val[result$spill_loc_idx])
  # 
  # # get ids of neighboring cells
  # id_neighbour_dt <- data.table(adjacent(facc_fill_r, result$spill_loc_idx, directions = 8, pairs=T))
  # # #associate depression no with the data.table
  # id_neighbour_dt[, dep_no := result[.SD, on = .(spill_loc_idx = from), x.dep_no]]
  # # #sort by dep no to put all neighbors together
  # # id_neighbour_dt <- setorder(id_neighbour_dt, dep_no)
  # #get D8 flow direction code for central cell
  # id_neighbour_dt$DS_direction <- d8_pntr_val[id_neighbour_dt$from]
  # #generate D8 flow direction code for neighbor cell (D8 rule: 64, 128, 1, 32, 2, 16, 8, 4; 'NW', 'N', 'NE', 'W', 'E', 'SW', 'S', 'SE')
  # id_neighbour_dt$ngbr_direction <- rep(c(64, 128, 1, 32, 2, 16, 8, 4), length(result$spill_loc_idx))
  # # keep rows where DS_direction = ngbr_direction (i.e., centeral cell direction to ds cell)
  # id_neighbour_dt <- id_neighbour_dt[DS_direction == ngbr_direction]
  # # get the next DS depression based on depression watershed and to value from id_neighbour
  # id_neighbour_dt$DS_dep_no <- dep_no_dpth_dt$dep_basin[id_neighbour_dt$to]
  # # add this to the result dt
  # result$next_downstream_dep <- id_neighbour_dt$DS_dep_no
  # 
  # # # ##**##**##**##**##**
  # this block does the same as the above but it's problematic as it
  # doesn't follow the D8 flow directions that are used to delineate the NCA
  # # old based on max_facc of neighboring cells (problematic)
  # id_neighbour_dt <- adjacent(facc_fill_r, result$spill_loc_idx, directions = 8)
  # # #associate depression no with the data.table
  # id_neighbour_dt <- data.table(id_neighbour_dt,dep_no=rep(result$dep_no,8))
  # #get facc based on all `to` cells and get max by dep_no
  # id_neighbour_dt$facc_fill <- dep_no_dpth_dt$facc_fill[id_neighbour_dt$to]
  # next_ds_dep <-  id_neighbour_dt[, .(max_next_facc_idx = id_neighbour_dt$to[.I[which.max(facc_fill)]]), by=dep_no]
  # next_ds_dep$ds_dep <- dep_no_dpth_dt$dep_basin[next_ds_dep$max_next_facc_idx]
  # # update next downstream depression id
  # result$next_downstream_dep <- next_ds_dep$ds_dep
  # ##**##**##**##**##**
  result$solve_order <- dep_no_dpth_dt$facc_fill[result$spill_loc_idx]
  #order data.table by solve order to make it from upstream to downstream
  result <- result[order(solve_order)]
  # re-number solve order from 1:tot_ndep
  result$solve_order <- 1:tot_ndep
  
  # Remove the "spill_loc" column from the data.table
  result[, spill_loc_idx := NULL]
  
  # get index (row location) of the nextds depression in the data.table
  result <- result[, next_ds_dep_index := match(next_downstream_dep, dep_no)]
  
  # check for self draining depressions or depressions that drain to an upstream depression
  result$dep_cascade_err <- 0 # 0 means no problem, 1 means either self draining or drains to an upstream depression
  result$row_num <- 1:nrow(result)
  result$dep_cascade_err[result$next_ds_dep_index<=result$row_num] <- 1
  
  # get number of depressions with problems
  ndep_cascade_error <- sum(result$dep_cascade_err)
  # fix this by setting the next ds depression to NA (i.e., drain to the river/outlet)
  result$next_downstream_dep[result$dep_cascade_err == 1] <- NA
  result$next_ds_dep_index[result$dep_cascade_err == 1] <- NA
  
  if (ndep_cascade_error >= 1) {
    cat('********************************* \n')
    cat('⚠------------WARNING-----------⚠ \n')
    cat('********************************* \n')
    cat('Found self draining or upstream draining depressions. \n') 
    cat(ndep_cascade_error, 'out of', nrow(result), '->', sprintf("%.3f",(ndep_cascade_error/nrow(result))*100), '% \n')
    cat('This is fixed by setting the next downstream depression of these depressions to NA. \n')
    cat('This means they drain directly to the river. \n')
    cat('Please check the DEM and outputs of the pre-processing if you want to change this or to identify the actual downstream depression. \n')
  }
  
  #drop extra columns
  result <- result[ ,`:=`(row_num = NULL)]
  
  
  
  # #write summary to csv file
  write.csv(result,file = paste0(out_dir,'/depressions_summary.csv'),row.names = F)
  
  end_time = Sys.time()
  cat('Finished generating the depressions summary table in', as.numeric(end_time-start_time,units="secs"), 'sec \n')

  return(result)
  
}

#########################################
plot_NCA_depr_stor_map <- function(basin_data, NCA_ras, depr_stor_ras, rain, return_period, out_dir, itime){
  
  basin_bound <- st_read(basin_data, quiet = TRUE)
  
  #base R plot (fast)
  
  # Open PNG device
  jpeg(paste0(out_dir,'/map_',itime,'.jpeg'), width = 9, height = 5, units = "in", res = 300)
  
  
  plot(NCA_ras, col = adjustcolor("red", alpha = 0.5), 
       legend=FALSE, xaxt='n', yaxt='n',
       main = paste0('Return Period = ',return_period,' years', " & Rainfall depth = ", rain, ' mm'))
  # plot(dep_ras, add=TRUE, col = terrain.colors(50), legend = FALSE)
  plot(depr_stor_ras, add=TRUE, legend=TRUE, col = viridis(10))
  mtext("Water depth (m)", side = 4, cex = 1, line = 1.05)
  plot(basin_bound$geometry, add=TRUE, border = "black", col=NA)
  
  # Add legends
  # Create a custom legend combining all elements
  legend("topright", 
         legend = c("NCA", "Boundary"), 
         col = c(adjustcolor("red", alpha = 0.5), NA),
         fill = c(adjustcolor("red", alpha = 0.5), NA),
         border = c(NA, "black"),
         title = "Legend"
  )
  
  
  dev.off()
  
  ## ggplot approach (slow)
  # # library("tidyverse")
  # # library("ggnewscale")# use multiple scales
  # # library('viridis')
  # # library("stars")
  # 
  # 
  # 
  # # depth_stack_df$type <- "depth"
  # min_depth <- 0.001
  # max_depth <- 5 #min(max(depth_stack[[nrow(RP_depth)]])@data@max, 10)
  # basin_bound <- st_read(basin_data, quiet = TRUE)
  # # get the basin boundary
  # # basin_bound <- DEM_ras
  # # basin_bound[!is.na(basin_bound)]=1
  # # basin_bound <- sf::st_as_sf(stars::st_as_stars(basin_bound), 
  # #                             as_points = FALSE, merge = TRUE) # requires the sf, sp, raster and stars packages
  # # plots <- list()
  # 
  # 
  # # convert the stack to a data.frame
  # NCA_df <- as.data.frame(rasterToPoints(NCA_ras))
  # colnames(NCA_df) <- c('x','y','val')
  # # remove NAs
  # NCA_df <- na.omit(NCA_df)
  # 
  # # convert the stack to a data.frame
  # depr_stor_df <- as.data.frame(rasterToPoints(depr_stor_ras))
  # colnames(depr_stor_df) <- c('x','y','val')
  # #remove NA
  # depr_stor_df <- na.omit(depr_stor_df)
  # 
  # 
  # map <- ggplot() +
  #   theme_void() +
  #   geom_sf(data=basin_bound, fill=NA, aes(color='Watershed boundary'))+
  #   scale_color_manual(name='', values = c('gray60'))+
  #   geom_tile(data = NCA_df,
  #             aes(x = x, y = y, fill = 'Non-Contributing Area'), alpha=0.3) +
  #   
  #   # scale_fill_viridis_c(name= 'NCA', option = "D")+
  #   scale_fill_manual(name = '', values = c('red'))+
  #   # geoms below will use another fill scale
  #   new_scale_fill() +
  #   geom_tile(data = depr_stor_df,
  #             aes(x = x, y = y, fill = val)) +
  #   scale_fill_viridis_c(name= 'Pond depth (m)', option = "D", limits = c(min_depth, max_depth))+
  #   labs(title = paste0('Return Period = ',return_period,' years', " & Rainfall depth = ", rain, ' mm'))+
  #   guides(fill = guide_colorbar(order=1),
  #          color=guide_legend(order = 4))
  # 
  # ggsave(filename = paste0(out_dir,'/map_',itime,'.jpeg'),plot = map, units = 'in', width = 9, height = 5, dpi = 300)
  # 
  # 
  # ## save the plots as gif using imagemagick from system
  # # system("convert -delay 100 -dispose Background -loop 0 *.jpeg anim_NCA.gif")
  
}


#########################################


simulate_depStorage_NCA <- function(config_file, dep_smr){
  # free up unused memory
  gc()
  
  #function arguments
  DEM_name <- paste0(config_file$inp_dir, '/', config_file$DEM_name)
  temp_dir <- paste0(config_file$inp_dir, '/', config_file$temp_dir)
  out_dir <- paste0(config_file$inp_dir, '/', config_file$out_dir)
  sim_type <- config_file$sim_type
  RP_rain <- read.csv(paste0(config_file$inp_dir, '/',config_file$RP_rain_file))
  # flag to activate the tracking of mass balance for each individual depression
  debug_flag <- FALSE
  
  # profvis({
  inc_depth <- RP_rain$rain_depth_mm/1000 # mm -> m
  return_period <- RP_rain$return_period
  # total number of depressions
  tot_ndep <- nrow(dep_smr)
  # outlet_location as cell with max flow accumulation
  out_loc <- which.max(getValues(raster(paste0(temp_dir, '/flow_acc_fill.tif'))))
  # read depression no raster to use as a template for newer rasters
  dep_no_r <- raster(paste0(temp_dir,'/depression_no.tif'))
  #depressions filling state
  dep_state <- data.frame(dep_smr$dep_no,dep_smr$dep_vol_m3,matrix(NA,nrow = tot_ndep,ncol = length(inc_depth)))
  colnames(dep_state) <- c('dep_no','dep_max_vol_m3', paste0('cur_dep_stor_m3_ts',1:length(inc_depth)))
  #define initial conditions of the depressions as empty
  ini_cond <- rep(0,tot_ndep)
  # initialize depressional storage vector
  cur_dep_stor <- rep(0,tot_ndep)
  # initialize average water level per depression
  new_elev <- rep(0,tot_ndep)
  # initialize fraction depth (pondstorage/depressionstorage) for each depression
  frac_depth <- rep(0,tot_ndep)
  # initialize outflow volume
  outflow_volume <- rep(0,length(inc_depth))
  #faster indexing of depression no
  dep_no_val <- getValues(raster(paste0(temp_dir,'/depression_no_non_zeros.tif')))
  # dep_no_val_dt <- data.table(dep_no_val)
  # 
  # #get indexes of cells within each dep no. outside loop to save time.
  # ## This is slow
  # # ind_list <- lapply(dep_smr$dep_no, function(idep_no_value) {
  # #   dep_no_val_dt[dep_no_val == idep_no_value, which = TRUE]
  # # })
  # # ########
  # # faster indexing compared to the above
  # # Create an index column to keep track of the original row number
  # dep_no_val_dt[, idx := .I]
  # # Find all matches for each dep_no using a join operation
  # result <- dep_no_val_dt[.(dep_no_val = dep_smr$dep_no), on = "dep_no_val", .(dep_no = i.dep_no_val, idx)]
  # # # Create a column to store the order of appearance of each dep_no
  # # result[, dep_no_order := seq_len(.N), by = dep_no]
  # # Convert the result to a list and maintain the same dep_no order as in dep_smr
  # ind_list <- split(result$idx, factor(result$dep_no, levels = dep_smr$dep_no))
  # rm(result)
  # #drop extra columns
  # dep_no_val_dt <- dep_no_val_dt[ ,`:=`(idx = NULL)]
  # ########
  DEM_val <- getValues(raster(DEM_name))
  cat('---------------------------- \n')
  cat('Starting the simulation','\n')
  cat('---------------------------- \n')
  
  for(itime in 1:length(inc_depth)){
    cat('Simulating event number ', itime,' with a return period of ', return_period[itime],' years and depth of ',inc_depth[itime]*1000, ' mm \n')
    start_time <- Sys.time()
    # #print(di)
    # #volume of filling for each depression
    inc_vol <- inc_depth[itime]*dep_smr$basin_area_m2
    # Mass Balance for the depressions (used for debugging)
    if(debug_flag){
      #MB_dep <- data.frame(dep_smr$dep_no,dep_smr$dep_vol_m3, matrix(NA,nrow = tot_ndep,ncol = 6))
      #colnames(MB_dep) <- c('dep_no','dep_max_vol_m3', 'ini_stor', 'final_stor', 'inflow', 'dep_outflow', 'to_river', 'MB')
      idep_ini_stor <- idep_final_stor <- idep_inflow <- idep_outflow <- idep_to_river <- idep_MB <- rep(0,tot_ndep)
      }
    #Initialize variables
    wl_new <- DEM_val #empty conditions
    non_fill_dep <- dep_no_val #this to mark depressions that are filled with NA and leave those unfilled

    
    # start of the loop
    for (idep in 1:tot_ndep){ #loop through depression from most upstream to downstream by solving order
      
      excess_vol <- 0
      ####filling depressions
      # get current depression information
      inc_voli <- inc_vol[idep]
      idep_no <- dep_smr$dep_no[idep]
      ## this is the code bottleneck (look for ways to improve this)
      ## ind <- dep_no_val_dt[dep_no_val == idep_no, which = TRUE]
      # ind <- ind_list[[idep]] #much faster
      max_vol <- dep_smr$dep_vol_m3[idep]
      inext_dep <- dep_smr$next_downstream_dep[idep]
      
      inext_dep_id <- dep_smr$next_ds_dep_index[idep] #which(dep_smr$dep_no==inext_dep)
      
      max_dep_d <- dep_smr$max_depth_m[idep]
      max_d_idx=dep_smr$max_d_id[idep] #cell id with max depth in each depression
      # calculate total volume within the depression
      tot_vol <- ini_cond[idep]+inc_voli
      
      if(tot_vol>max_vol){ # depression is full
        cur_dep_stor[idep] <- max_vol
        frac_depth[idep] <- 1 # this means that the dep is 100% full
        excess_vol <- tot_vol-max_vol #excess volume (extra volume after filling the depression)
        #move excess volume to the next downstream dep or move it as outflow if there is not downstream depression
        if(is.na(inext_dep)){ #if no downstream move excess volume to outflow
          outflow_volume[itime] <- outflow_volume[itime]+excess_vol
          if(debug_flag){idep_to_river[idep] <- excess_vol}
        }else{ #move excess volume to next downstream depression
          inc_vol[inext_dep_id] <- inc_vol[inext_dep_id]+excess_vol
          if(debug_flag){idep_to_river[idep] <- 0}
        }
        ################################
      }else{ # depression is not full
        cur_dep_stor[idep] <- tot_vol
        frac_depth[idep] <- cur_dep_stor[idep]/max_vol
      }
      # mass balance checks (by depression)
      if(debug_flag){
        idep_ini_stor[idep] <- ini_cond[idep]
        idep_final_stor[idep] <- cur_dep_stor[idep]
        idep_inflow[idep] <- inc_voli
        idep_outflow[idep] <- excess_vol
        idep_MB[idep] <- idep_inflow[idep] - idep_outflow[idep] - (idep_final_stor[idep] - idep_ini_stor[idep])
      }
      
      #new elevation of cell with max depth in each depression
      new_elev[idep]=DEM_val[max_d_idx]+(frac_depth[idep]*max_dep_d)
      # # This section makes the loop slow (moved to be processed by data.table to be faster)
      ## currently implemented outside the loop
      # if(new_elev[idep]>=dep_smr$max_wl_m[idep] || frac_depth[idep]==1){ #check to limit the filling to the max possible depth
      #   #this means that the depression is full
      #   new_elev[idep]=dep_smr$max_wl_m[idep]
      #   wl_new[ind]=new_elev[idep]
      #   #no need to loop on cells in the depression as the depression is full and wl is the max
      #   non_fill_dep[ind] <- NA #this means that this depression is filled
      # } else{
      #   #this check is done to avoid reducing dem at potholes edges (current elevation > new elevation)
      # 
      #   update_condition <- wl_new[ind] <= new_elev[idep]
      #   # Update wl_new for cells where the condition is true
      #   wl_new[ind[update_condition]] <- new_elev[idep]
      # }
    }
    
    ##################
    # perform depressional storage flood (depth) mapping
    # create a data.table with new_elev and frac_depth
    calc_wl_frac <- data.table(dep_no=dep_smr$dep_no, new_elev, frac_depth, max_wl_m=dep_smr$max_wl_m)
    rasters_data <- data.table(dep_no_val, DEM_val, wl_new=DEM_val, non_fill_dep=dep_no_val)
    # Create a temporary column containing row indices to preserve cells order
    rasters_data[, index := .I]
    
    # Calculate new_elev and condition in calc_wl_frac
    calc_wl_frac[, new_elev := ifelse(new_elev >= max_wl_m | frac_depth == 1, max_wl_m, new_elev)]
    calc_wl_frac[, condition := new_elev >= max_wl_m | frac_depth == 1]
    # Left join raster_data with calc_wl_frac on dep_no_val and dep_no
    rasters_data <- merge(rasters_data, calc_wl_frac[, .(dep_no, new_elev, condition)], by.x = "dep_no_val", by.y = "dep_no", all.x = TRUE)
    # Reorder raster_data based on the temporary index column
    setkey(rasters_data, index)
    
    # Update wl_new and non_fill_dep in raster_data based on condition
    rasters_data[condition == TRUE, `:=` (wl_new = new_elev, non_fill_dep = NA)]
    rasters_data[condition == FALSE, `:=` (wl_new = new_elev)]
    # rasters_data[, wl_new := new_elev]
    # rasters_data[condition == TRUE, `:=` (non_fill_dep = NA)]
    
    # this check is done to avoid reducing dem at potholes edges (current elevation > new elevation)
    # Update wl_new if wl_new < DEM_val
    rasters_data[wl_new < DEM_val, wl_new := DEM_val]
    # Replace NA values in non_fill_dep with zero
    rasters_data[is.na(non_fill_dep), non_fill_dep := 0]
    ###add output cell to the non_fill_dep
    ###assume that max flow accumulation is the outlet
    rasters_data$non_fill_dep[out_loc] <- tot_ndep+1e6
    ##################
    #bookkeeping
    dep_state[,2+itime] <- cur_dep_stor
    # # replace NA with zero so that WBT can run it
    # non_fill_dep[is.na(non_fill_dep)] <- 0
    # ###add output cell to the non_fill_dep
    # ###assume that max flow accumulation is the outlet
    # non_fill_dep[out_loc] <- tot_ndep+1e6
    
    
    #-----------------------
    ######################
    #note that the generated water raster file is approximation (not actual storage)
    #it was generated using interpolation (frac_depth)
    #####################
    #mass balance check (second check done at the watershed scale)
    #stored water in depressions
    cur_basn_storage <- sum(cur_dep_stor)
    chng_in_storage <- cur_basn_storage-sum(ini_cond)
    # applied volume
    app_vol <- sum(inc_depth[itime]*dep_smr$basin_area_m2)
    # Mass balance error in meters
    MB <- (app_vol-(chng_in_storage+outflow_volume[itime]))/sum(dep_smr$basin_area_m2)
    #cat('Finished event number',di,'\n')
    cat('Mass Balance(In-ds-out)=',sprintf("%.3f",MB), 'm\n')
    if(debug_flag){cat('Mass Balance(In-ds-out) per depression=',sprintf("%.3f",sum(idep_MB)), 'm^3\n')}
    # cat('---------------------------- \n')
    #----------------------
    
    # bookeeping of the initial conditions
    if(sim_type == "timeseries"){
      # assign the initial condition as the final state of current time step
      ini_cond <- cur_dep_stor #dep_state[,2+di]
      # otherwise leave the ini_cond unchanged (()zero)
    }
    
    
    #create new raster for the wl_new
    #dreate dummy rasters
    wl_new_ras <- dep_no_r
    dep_new_ras <- dep_no_r
    nonfilled_dep_ras <- dep_no_r
    
    values(wl_new_ras) <- rasters_data$wl_new #wl_new
    values(dep_new_ras) <-  rasters_data$wl_new - rasters_data$DEM_val #wl_new-DEM_val
    values(nonfilled_dep_ras) <-rasters_data$non_fill_dep  #non_fill_dep
    dep_new_ras <- reclassify(dep_new_ras, cbind(-Inf, 0.001, NA)) #dep_new_ras[dep_new_ras<0.001] <- NA #remove small depths to show only depressions
    
    #write new filled depression wl file
    writeRaster(wl_new_ras,
                paste0(temp_dir, '/wl_',itime,'_',return_period[itime],'yr_',sprintf("%.3f",inc_depth[itime]*1000),'_mm','.tif'),
                overwrite=T,options=c("COMPRESS=NONE", paste0("NUM_THREADS=",config_file$ncores)))
    writeRaster(dep_new_ras,
                paste0(out_dir,'/depr_stor_dep_',itime,'_',return_period[itime],'yr_',sprintf("%.3f",inc_depth[itime]*1000),'_mm','.tif'),
                overwrite=T,options=c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9', paste0("NUM_THREADS=",config_file$ncores)))
    writeRaster(nonfilled_dep_ras,
                paste0(temp_dir,'/3-nonfill_depressions_',itime,'_',return_period[itime],'yr_',sprintf("%.3f",inc_depth[itime]*1000),'_mm','.tif'),
                overwrite=T,options=c("COMPRESS=NONE", paste0("NUM_THREADS=",config_file$ncores)))
    
    #delineate watesheds based on the new depressions (calc NCA)
    wbt_watershed(wd = temp_dir, d8_pntr = paste0(temp_dir,'/D8_flow_pointer_filled.tif'),
                  pour_pts = paste0(temp_dir,'/3-nonfill_depressions_',itime,'_',return_period[itime],'yr_',sprintf("%.3f",inc_depth[itime]*1000),'_mm','.tif'),
                  output = paste0(temp_dir,'/watershed_',itime,'_',return_period[itime],'yr_',sprintf("%.3f",inc_depth[itime]*1000),'_mm','.tif'),verbose_mode = F)
    #quantify the contributing area (cell have 'tot_ndep+1e6' value in the watershed file)
    NCA_ras <- raster(paste0(temp_dir,'/watershed_',itime,'_',return_period[itime],'yr_',sprintf("%.3f",inc_depth[itime]*1000),'_mm','.tif'))
    CA_ras <- NCA_ras # copy the same raster as contributing to preserve extents
    # slow
    # NCA_ras[NCA_ras != tot_ndep+1e6 & !is.na(NCA_ras) ] <- 0 #non-contributing
    # NCA_ras[NCA_ras == tot_ndep+1e6] <- NA #contributing
    NCA_ras <- reclassify(NCA_ras, cbind(tot_ndep+1e6, NA)) # Replace values of tot_ndep+1e6 with NA
    NCA_ras <- reclassify(NCA_ras, cbind(-Inf, tot_ndep+1e6, 0)) # Replace all other values with 0
    
    # Create a CA raster
    # slow
    # CA_ras[CA_ras != tot_ndep+1e6 & !is.na(CA_ras) ] <- NA #non-contributing
    # CA_ras[CA_ras == tot_ndep+1e6] <- 1 #contributing
    CA_ras <- reclassify(CA_ras, cbind(-Inf, tot_ndep+1e6-1, NA)) # Replace all values other than tot_ndep+1e6 with NA
    CA_ras <- reclassify(CA_ras, cbind(tot_ndep+1e6, 1)) # Replace values of tot_ndep+1e6 with 1
    #
    #write the NCA_ras & CA_ras
    writeRaster(NCA_ras,
                paste0(out_dir,'/NCA_',itime,'_',return_period[itime],'yr_',sprintf("%.3f",inc_depth[itime]*1000),'_mm','.tif'),
                overwrite=T,options=c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9', paste0("NUM_THREADS=",config_file$ncores)))
    writeRaster(CA_ras,
                paste0(out_dir,'/CA_',itime,'_',return_period[itime],'yr_',sprintf("%.3f",inc_depth[itime]*1000),'_mm','.tif'),
                overwrite=T,options=c('COMPRESS=DEFLATE', 'PREDICTOR=2', 'ZLEVEL=9', paste0("NUM_THREADS=",config_file$ncores)))
    
    # plot the NCA and depressional storage
    if(config_file$map_plot){
      plot_NCA_depr_stor_map(basin_data=paste0(config_file$inp_dir, '/', config_file$basin_bound), 
                             NCA_ras = NCA_ras, depr_stor_ras = dep_new_ras, 
                             rain = RP_rain$rain_depth_mm[itime], 
                             return_period = RP_rain$return_period[itime], out_dir = out_dir, itime = itime)
      # free up unused memory
      gc()
    }
    end_time = Sys.time()
    cat('Finished simulating event number ', itime,'in', as.numeric(end_time-start_time,units="secs"), 'sec \n')
    cat('---------------------------- \n')
    }
  write.csv(dep_state, paste0(out_dir,'/depression_state.csv'), row.names = F)
  
  return(dep_state)
  # })
}
#####################################################
