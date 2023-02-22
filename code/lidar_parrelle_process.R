###########################################################################################
#
#        this script determines the environmental limiting factors for shrub fcovr
#
#    --- Last updated:  2020.11.13 By Daryl Yang <dediyang@bnl.gov>
###########################################################################################

#******************** close all devices and delete all variables *************************#
rm(list=ls(all=TRUE))   # clear workspace
graphics.off()          # close any open graphics
closeAllConnections()   # close any open connections to files
dlm <- .Platform$file.sep # <--- What is the platform specific delimiter?
#*****************************************************************************************#

#****************************** load required libraries **********************************#
### install and load required R packages
list.of.packages <- c("lidR", 'ggplot2', 'future', 'rgdal')  
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
out.dir <- "D:\\lidar_playwith\\output"
### Create output folders
if (! file.exists(out.dir)) dir.create(out.dir,recursive=TRUE)
### create a temporary folder
temp.dir <- paste0(out.dir, '/', 'temp')
if (! file.exists(temp.dir)) dir.create(temp.dir,recursive=TRUE)
#*****************************************************************************************#

#*************************************** load data ***************************************#
# load all laz file as catalog
data.dir <- 'D:\\lidar_playwith\\input'
ctg <- readLAScatalog(data.dir)

### set up paralell processing
future::plan(multisession)
set_lidr_threads(8)
get_lidr_threads()
opt_chunk_size(ctg) <- 50
opt_chunk_buffer(ctg) <- 0

### split laz file into tiles
opt_output_files(ctg) <- paste0(temp.dir, "\\retile_{XLEFT}_{YBOTTOM}")
newctg = catalog_retile(ctg)

ctg <- readLAScatalog(temp.dir)
las_check(ctg)

dens = 30
opt_output_files(ctg) <- paste0(temp.dir, "\\{*}_thinned")
thinned_ctg <- decimate_points(ctg, homogenize(dens, 1))

### classify point cloud
ws <- seq(3, 30, 6)
th <- seq(0.1, 0.3, length.out = length(ws))
opt_output_files(thinned_ctg) <- paste0(temp.dir, "\\{*}_classified")
classified_ctg <- classify_ground(thinned_ctg, algorithm = pmf(ws = ws, th = th))

### create dtm
opt_output_files(classified_ctg) <- opt_output_files(classified_ctg) <- paste0(temp.dir, "\\{*}_dtm")
opt_stop_early(classified_ctg) <- FALSE
dtm <- rasterize_terrain(classified_ctg, 1, tin())

### create canopy height laz file
opt_output_files(classified_ctg) <-  paste0(temp.dir, "\\{*}_norm")
#opt_merge(classified_ctg) <- TRUE
opt_stop_early(classified_ctg) <- FALSE
ctg_norm <- normalize_height(classified_ctg, dtm)
#opt_merge(ctg) <- TRUE

linear_stretch <- function(las) 
{
  las$Z <- las$Z*10
  return(las)
}
new_ctg <- catalog_map(ctg_norm, linear_stretch)

### identify shrubs
f <- function(x) {x * 0.1 + 3}
heights <- seq(0,5,1)
ws <- f(heights)
#opt_output_files(new_ctg) <- ""
#ttops <- locate_trees(new_ctg, lmf(f))
opt_output_files(new_ctg) <- "" #paste0(out.dir, "/{*}_treeid")
#opt_merge(new_ctg) <- TRUE
ttops <- locate_trees(new_ctg, lmf(f), uniqueness = "bitmerge")
#ttops$treeID <- 1:nrow(ttops)
ttops_merge <- do.call(rbind, ttops)
outname <- paste0(out.dir, '/', 'flight5_tree_id.shp')
#writeOGR(ttops_merge, outname, driver="ESRI Shapefile")
st_write(ttops_merge, outname, driver = "ESRI Shapefile")

### segment individual shrubs
opt_output_files(new_ctg) <- paste0(temp.dir, "/chm_{*}")
opt_merge(new_ctg) <- TRUE
opt_stop_early(new_ctg) <- FALSE
chm <- rasterize_canopy(new_ctg, 1, p2r(0.15), overwrite=TRUE)

opt_output_files(new_ctg) <- paste0(temp.dir, "/{*}_segmented")
algo <- dalponte2016(chm, ttops_merge)
opt_merge(new_ctg) <- TRUE
opt_stop_early(new_ctg) <- FALSE
ctg_segmented <- segment_trees(new_ctg, algo)
#outname <- paste0(data.dir, '/', 'tree_segment.shp')
#writeOGR(ctg_segmented, outname, driver = "ESRI Shapefile")
#algo <- dalponte2016(chm, ttops)
#ctg_segmented <- catalog_map(new_ctg, segment_trees(algo))

#ctg_segmented <- segment_trees(ctg_norm, algo)
outname <- paste0(out.dir, '/', 'flight5_tree_seg')
opt_output_files(ctg_segmented) <- outname
opt_chunk_buffer(ctg_segmented) <- 0
opt_chunk_size(ctg_segmented) <- 10000
singlefile_ctg <- catalog_retile(ctg_segmented)

temp <- crown_metrics(singlefile_ctg, func = .stdtreemetrics, geom = "convex")

unlink(temp.dir, recursive = T, force = T)





