###########################################################################################
#
#        this script is for creating canopy height model from sfm or lidar point clouds
#
#    --- Last updated:  2023.02.15 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("lidR", 'ggplot2', 'future', 'rgdal', 'terra')  
# check for dependencies and install if needed
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, dependencies=c("Depends", "Imports",
                                                                       "LinkingTo"))
version_requirements <- c("3.3.2")
if (!packageVersion("ggplot2") >= version_requirements[1]) {
  remotes::install_version(package="ggplot2", version=paste0(">= ", version_requirements), 
                           dependencies=c("Depends", "Imports", "LinkingTo"), upgrade="ask",
                           quiet=TRUE)
}
# load libraries
invisible(lapply(list.of.packages, library, character.only = TRUE))
#*****************************************************************************************#

#************************************ user parameters ************************************#
out.dir <- "/Users/darylyang/Desktop/gdrive/github/lidar_tools/output"
### Create output folders
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
### create a temporary folder
temp.dir <- paste0(out.dir, '/', 'temp')
if (! file.exists(temp.dir)) dir.create(temp.dir,recursive=TRUE)

### define point cloud density for generating canopy height and other stats. this is 
### necessary when super dense cloud is not necessary, but not affect the result
reducer = 'YES'
if (reducer == 'YES')
{
  dens = 30 # point/m2
}

### define if you need to classify the point clouds. YES for classification, NO for don't
### classify (which means the point cloud is already classified)
classify = 'YES'

### define output resolution of terrain and canopy height raster
reso = 0.5
#*****************************************************************************************#

#*************************************** load data ***************************************#
### load all laz files as catalog
data.dir <- '/Users/darylyang/Desktop/gdrive/github/lidar_tools/data'
ctg <- readLAScatalog(data.dir)

### set up parallel processing, please change these as necessary
future::plan(multisession) 
set_lidr_threads(5)
opt_chunk_size(ctg) <- 20
opt_chunk_buffer(ctg) <- 0

### split laz files into new tiles. this is not necessary if the original laz files are
### already tiles
opt_output_files(ctg) <- paste0(temp.dir, "/retile_{XLEFT}_{YBOTTOM}")
newctg = catalog_retile(ctg)

### reload the new tiles and check
ctg <- readLAScatalog(temp.dir)
las_check(ctg)
#*****************************************************************************************#

#***************************************** get chm ***************************************#
### reduce point cloud density if needed
if (reducer == 'YES')
{
  opt_output_files(ctg) <- paste0(temp.dir, "/{*}_thinned")
  thinned_ctg <- decimate_points(ctg, homogenize(dens, 1))
} else 
{
  thinned_ctg <- ctg
}

### classify point cloud if needed
if (classify == 'YES')
{
  ### user can tweak the parameters as needed
  ws <- seq(3, 30, 6)
  th <- seq(0.1, 0.3, length.out = length(ws))
  opt_output_files(thinned_ctg) <- paste0(temp.dir, "/{*}_classified")
  classified_ctg <- classify_ground(thinned_ctg, algorithm = pmf(ws = ws, th = th))
} else
{
  # reassign classified ctg to ctg
  classified_ctg <- thinned_ctg
}

### create digital terrain model raster ile
opt_output_files(classified_ctg) <-  paste0(temp.dir, "/{*}_dtm")
opt_stop_early(classified_ctg) <- FALSE
dtm <- rasterize_terrain(classified_ctg, reso, tin())

### create canopy height laz file
opt_output_files(classified_ctg) <-  paste0(temp.dir, "/{*}_norm")
opt_stop_early(classified_ctg) <- FALSE
ctg_norm <- normalize_height(classified_ctg, dtm)

### create canopy height raster file
opt_output_files(ctg_norm) <- paste0(temp.dir, "/chm_{*}")
opt_stop_early(ctg_norm) <- FALSE
opt_merge(ctg_norm) <- TRUE
chm <- rasterize_canopy(ctg_norm, reso, p2r(0.15), overwrite=TRUE)

dtm.filename <- paste0(out.dir, '/', 'digital_terrain_model.tif')
writeRaster(dtm, dtm.filename, overwrite = TRUE)

chm.filename <- paste0(out.dir, '/', 'canopy_height_model.tif')
writeRaster(chm, chm.filename, overwrite = TRUE)

### remove the temp folder
unlink(temp.dir, recursive = T, force = T)
#*****************************************************************************************#




