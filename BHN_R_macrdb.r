#########################################
# Cost accumulation mapping in R
#
# Andy Nelson. Professor, Spatial Agriculture and Food Security
# Department of Natural Resources (Room 4-144)
# ITC - Faculty of Geo-Information Science and Earth Observation of the University of Twente
# PO Box 217, 7500 AE Enschede, the Netherlands

# Uses scripts or ideas by 
# Dan Weiss, Malaria Atlas Project, University of Oxford
# Jacob van Etten, Bioversity International

# Adaptations for Bhutan case by
# Rolf de By, GIP/ITC - Faculty of Geo-Information Science and Earth Observation 
# University of Twente
# Adaptations:
# -1- Operating system independence (win/macOS/linux)
# -2- Implementation on a 150m grid basis using national crs, as the data allows for this
# -2- Earlier move to metric crs in the script
# -3- Taken hard-coded data out and into external (csv) files.  
#     Such as: landcover_class_speeds.csv
# -4- Inclusion of a Tobler function for vehicle travel in sloped terrain
# -5- Repair of a which.min() induced error that assigned a zone id to NaN territory
# -6- Implementation of a tail end of the workflow, working with population census data
#     and housing information that allows fine-grain match of population with travel time
#     zones.
########################################
# call the required libraries (packages)
if (!require("raster"))
  install.packages("raster")
if (!require("gdistance"))
  install.packages("gdistance")
if (!require("rgdal"))
  install.packages("rgdal")

library(raster)
library(gdistance)
library(rgdal)

##########################################################################################
# A function we will need to test which OS we have:
get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

##########################################################################################
# You may need A LOT of RAM for this, depends on size of region and number of targets.
# Next command is only needed for Windows.  Other platforms allow dynamic mem assignation.
if (get_os()=='windows') memory.limit(999999)
if (get_os()=='windows') dirsep <- "\\\\" else dirsep <- "/"

filepath <- function(fpath){
  gsub("/",dirsep,fpath)
  }

# Be conservative on raster's memory usage; terrain() has caused crashes if we do not set
# chunksize.  This value has not been optimized and may work if bigger.  But all depends
# on your hardware.
rasterOptions(chunksize=1e+06)

##########################################################################################
# set the working directory
setwd(filepath("C:\\UT\\Work\\Workspace\\CSI\\accessibility_course_bhutan-master"))

##########################################################################################

# Create a template metric raster that we will use for various cookie cuttings
country_template <- raster()

# set the dimensions of the Bhutan area
dimensions <- extent(124000, 463000, 2953600, 3127600)

# extent it into the necessary dimensions
country_template <- setExtent(country_template, dimensions)

# define the resolution of the raster cells ( 150 meters)
res(country_template) <- 150

crs(country_template) <- '+init=EPSG:5266'

# give a value for each column of the raster
country_template[] <- 1:ncell(country_template)
country_template

##########################################################################################
# Start adding various data sets to the environment
# elevation raster
bhutan_dem_srtm <- raster(filepath("inputs/bhutan_dem_srtm.tif"))
bhutan_dem_srtm
plot(bhutan_dem_srtm, main="DTM")

# calculate the slope using the elevation in each point of the map
bhutan_dem_srtm_slope <- terrain(bhutan_dem_srtm, opt='slope', unit='degrees')
#plot(bhutan_dem_srtm_slope)

# write the raster with slope
writeRaster(bhutan_dem_srtm_slope, filename=filepath("processing/bhutan_srtm_slope.tif"),
            format="GTiff", overwrite=TRUE)

# landcover raster
landcover <- raster(filepath("inputs/landcover2010.tif"))
landcover
plot(landcover, main="land cover classes")

# target localities for accessibility; here we have chosen 23 hospitals
targets <- readOGR(filepath("inputs/bhu_facilities_point.shp"))
# transform targets also to the crs of the country raster template
targets <- spTransform(targets, crs(country_template))
targets
plot(targets, main="location of 23 hospitals")

# reproject landcover and elevation rasters to national Bhutan metric system and resample to 150m 
landcover_reprojected <- projectRaster(landcover, country_template, 
       crs = '+init=EPSG:5266', res = 150)

slope_reprojected <- projectRaster(bhutan_dem_srtm_slope, country_template, 
       crs = '+init=EPSG:5266', res = 150)

landcover
landcover_reprojected

bhutan_dem_srtm_slope
slope_reprojected

##########################################################################################
# This all is for walking the terrain, depending on landcover and slope.
# convert slope and landcover to speed rasters 
# Landcover -  reclassify the landcover classes into walking speeds (km/h)
#  values >= 0 and < 1 become 1, values >= 1 and < 2 become 3, etc.
#  do not use speed = 0 km/hr for rather obvious reasons
# See example data file, csv must be a 5-column (id,from,to,becomes,description) file,
# from which we immediately drop the description.
lcs = read.table("inputs/landcover_class_speeds.csv", 
                 colClasses=c("integer","integer","integer","numeric","character"), 
                 header=TRUE, sep=",", row.names=1)
lcs$description <- NULL
# view table for checking:
lcs

# apply the reclassification to get travel speeds per landcover and write to file
landcover_speed <- reclassify(landcover_reprojected,lcs,right=FALSE)
writeRaster(landcover_speed, filename=filepath("processing/landcover_speed.tif"), 
            format="GTiff", overwrite=TRUE)
plot(landcover_speed, main="speed in km per hour over land cover classes")

##########################################################################################
# Step 5

# slope for walking
# Tobler's walking speed is given by W = 6e^-3.5|tan(slope)+0.05|
# 6 km/h is maxspeed at a slight downward angle; 0.05 rad is that optimal angle for best horizontal speed;
# -3.5 is a factor of decay of maximally attainable speed as slope angle changes.
slp_walk <-  6 * exp(-0.4 * abs(tan(slope_reprojected*pi/180) + 0.05))
slp_walk
writeRaster(slp_walk, filename=filepath("processing/sloped_walk.tif"), 
            format="GTiff", overwrite=TRUE)

# divide by base speed of 5km/h to get speed factor and write that to a file
writeRaster(slp_walk/6.0, filename=filepath("processing/sloped_factor_walking.tif"), 
            format="GTiff", overwrite=TRUE)
# Slope adjusted walking speeds over landcover, write to file
terrain_walk_spd <- landcover_speed * slp_walk/6.0
writeRaster(terrain_walk_spd, filename=filepath("processing/terrain_walk_speed.tif"), 
            format="GTiff", overwrite=TRUE)
plot(terrain_walk_spd, main="slope adjusted walking speed")        

##########################################################################################
# Step 6

# slope for car driving
# Experimental Tobler's car speed is given by W = 50e^-2.4|tan(slope)+0.12|
# Original -2.4 gave far too much decay.
slp_car <-  50 * exp(-0.4 * abs(tan(slope_reprojected*pi/180) + 0.12))
writeRaster(slp_car, filename=filepath("processing/sloped_cardrive.tif"), 
            format="GTiff", overwrite=TRUE)
plot(slp_car, main="slope adjusted car speed")        

# read the road network shapefile; in this version, we use OSM
# The integer64 option is required or values will be read as strings and handled as
# numbered classes.
road_shp <- readOGR(filepath("inputs/osm_roads.shp"),integer64="allow.loss")
# transform to crs of the country raster template
road_shp <- spTransform(road_shp, crs(country_template))
road_shp
plot(road_shp,main="roads")

# Our footpaths current version are not good.
footpath_shp <- readOGR(filepath("inputs/footpaths.shp"))
# transform to the crs of the country raster template
footpath_shp <- spTransform(footpath_shp, crs(country_template))
footpath_shp
plot(footpath_shp,main="footpaths")

# create an empty raster for road speed
road_spd <- raster()
# define its dimensions
road_spd <- setExtent(road_spd, dimensions)
# and resolution
res(road_spd) <- res(country_template)
# and crs
crs(road_spd) <- crs(country_template)

# same steps for footpath
footpath_spd <- raster()
footpath_spd <- setExtent(footpath_spd, dimensions)
res(footpath_spd) <- res(country_template)
crs(footpath_spd) <- crs(country_template)

road_spd
footpath_spd

# assign raster cell values using the maxspeed field for roads
#road_spd <- rasterize(x=road_shp,y=road_spd,field="maxspeed",fun=max)
#writeRaster(road_spd, filename=filepath("processing/road_speed.tif"), 
#            format="GTiff", overwrite=TRUE)
road_spd <- raster(filepath("processing/road_speed.tif"))
road_spd

# assign raster cell values using the maxspeed field for footpaths
footpath_spd <- rasterize(x=footpath_shp,y=footpath_spd, field="maxspeed", fun=max)
writeRaster(footpath_spd, filename=filepath("processing/footpath_speed.tif"), 
            format="GTiff", overwrite=TRUE)
# unclear why this is needed, but above only gives cell values equal to 1, and should be 5
footpath_spd[footpath_spd==1] <- 5
footpath_spd

# other modes of transporation (train, ship, ...) could go here:

##########################################################################################

# if slope needs to be in play for road and footpaths speeds, here is where that happens
# we are dividing by the base speed for walking (5kmph) and road (50kmph)
sloped_footpath_spd <- footpath_spd * slp_walk / 5.0
sloped_road_spd <- road_spd * slp_car / 50.0
writeRaster(sloped_road_spd, filename=filepath("processing/sloped_road_speed.tif"), 
            format="GTiff", overwrite=TRUE)
plot(sloped_road_spd)                

##########################################################################################
# Step 7

# merging the various (road, footpath, ...) *network* speed rasters; ensure that you 
# prioritize properly
road_network_spd <- merge(sloped_road_spd, sloped_footpath_spd)
# write a raster with the road network (all travel modes) speed 
writeRaster(road_network_spd, filename=filepath("processing/road_network_speed.tif"), 
            format="GTiff", overwrite=TRUE)

# Merge the speed components and convert from travel speed to travel cost
# road takes priority over rail, rail takes priority over ship travel, it over walking 
# the terrain ...
merged_spd <- merge(road_network_spd,terrain_walk_spd)
writeRaster(merged_spd, filename=filepath("processing/merged_speed.tif"), 
            format="GTiff", overwrite=TRUE)
plot(merged_spd, main="merged speed rasters in kmph") 

##########################################################################################
# Step 8

# convert speed in km per hr to travel time in minutes per metre
# THIS IS THE FRICTION SURFACE
friction <- 1.0 / (merged_spd * 1000 / 60.0 )
writeRaster(friction, filename=filepath("processing/friction.tif"), 
            format="GTiff", overwrite=TRUE)
plot(friction, main="friction layer reporting time in minutes to travel one metre") 

##########################################################################################
# Step 9

# Make the graph for 8 directions using 1/mean(merge_mins) as the conductance value 
# between neighbouring cells 
T <- transition(friction, function(x) 1/mean(x), 8) 
# geo-corrected version of the graph to divide the conductance by the real distance
# between pixels.
# FALSE OBSERVATION: In Bhutan case, the next line is not needed as we are already having
# metric rasters at 150m resolution.
# CORRECT OBSERVATION: GC is also needed when directions = 8 or 16 (see gdistance documentation)
T.GC <- geoCorrection(T,type="c")

# accumulated cost calculation to the nearest target using the geo-corrected graph and 
# the target points. This will need the use of T.GC instead if produced above.
access_mins <- accCost(T.GC, targets)
# write the resulting raster showing time in minutes to the nearest target
writeRaster(access_mins, filename=filepath("outputs/access_mins.tif"), 
            format="GTiff", overwrite=TRUE)
plot(access_mins, main="travel time in minutes to nearest hospital")
plot(targets, add=TRUE)

##########################################################################################
# Step 10
# Compute cost allocation  to map which pixel is closest to which target
# derived from https://stat.ethz.ch/pipermail/r-sig-geo/2011-July/012208.html
# no function for this so...make a stach of accumulated cost rasters, one per target
# and then use the minimum pixel value through the stack to identify which pixel
# is closest to which target

# set up a stack with 2 layers, use existing data to fill it for now
accCost_stack <- stack(access_mins,access_mins)
# run accumulative cost for each target (targets from shapefile (3 columns, [ID X Y])
for(i in 1:length(targets)) {
    # grab the X and Y coord
    t  <- cbind(targets@coords[i,1],targets@coords[i,2])
    # add accumulated cost to the stack
    accCost_stack <- stack(accCost_stack, accCost(T, t))
}
# remove first two layers which have dud info
# THIS STACK CAN ALSO BE SAVED AS A MULTIBAND RASTER WITH ONE TRAVEL TIME LAYER PER TARGET
accCost_stack <- dropLayer(accCost_stack,1:2)

# make a new raster based on the layer ID (pixel values will be from 1 to n layers)
# Original code had just which.min(accCost_stack) *but this assigns the ID code of one target shed
# to the (no data) area outside of study space too*.  In current understanding, this is a flaw
# of which.min (or possibly deriving from a raster stack feed with not quite robust data).
# Hence, an extra step with a mask operator.
# Original code line
accCost_minID <- which.min(accCost_stack)

# and next, mask out area outside of study; we use reprojected_landcover but this could be any
# proper raster with appropriate nodata delineation.
accCost_minID_masked <- mask(accCost_minID,reprojected_landcover)

# make a new raster based on minimum value and compare to raster from section 4.5 (should be the same)
accCost_min <- min(accCost_stack)

# make a reclass table to assign target IDs to layer IDs.
# layer ID 1 becomes target ID 1 etc...
m <- c(1, 2, targets@data$gid[1])
for(i in 2:length(targets)) {
    m <- c(m,i, i+1, targets@data$gid[i])
}

rclcA <- matrix(m, ncol=3, byrow=TRUE)
# apply the reclassification to get allocation zone IDs - THIS IS THE COST ALLOCATION ZONE MAP
accCost_zones <- reclassify(accCost_minID_masked,rclcA,right=FALSE)

# write cost allocation to file
writeRaster(accCost_zones, filename="outputs/access_alloc.tif", format="GTiff", overwrite=TRUE)
# write minimum access to file and compare to output from 4.5 - should be identical
writeRaster(accCost_min, filename="outputs/access_minimum.tif", format="GTiff", overwrite=TRUE)

plot(accCost_zones, main="catchments around each hospital")
plot(targets, add=TRUE)

##########################################################################################
# Population data
# The code below is likely somewhat baroque and will make decent R coders chuckle.
# Additional data inputs required are:
# -a- a fairly complete set of building centroids: buildingcentroid.shp
# -b- a topologically correct set of census tracts: censustract.shp
# -c- a table of average nr of people per building, per censustract: 
#     censustract_paxperbuilding.csv
#    (this one should really be superfluous and piggy-back on censustract.shp)
# -d- a table that provides breakpoints for traveltime zonation: traveltime_zonebreaks.csv

##########################################################################################

# first get the building centroids
# we represent buildings by centroids (more trivial assignation to raster cells)
building_shp <- readOGR(filepath("inputs/buildingcentroid.shp"),integer64="allow.loss")

# towards a raster that holds building counts
# create an empty raster for buildings
buildings <- raster()
buildings <- setExtent(buildings, dimensions)
res(buildings) <- res(country_template)
crs(buildings) <- crs(country_template)

# assign count of buildings to raster
buildings <- rasterize(x=building_shp,y=buildings,field="gid",fun='count',background=0)
writeRaster(buildings, filename=filepath("outputs/building_counts.tif"), format="GTiff", 
            overwrite=TRUE)
            
# towards a raster that holds census tract id
# create an empty raster 
censustracts <- raster()
censustracts <- setExtent(censustracts, dimensions)
res(censustracts) <- res(country_template)
crs(censustracts) <- crs(country_template)
     
# obtain the census tracts:
censustracts_shp <- readOGR(filepath("inputs/censustract.shp"),integer64="allow.loss")
# in below rasterize(), just using " field='id' " I could not make to work, so instead:
# ids <- as.integer(as.matrix(censustracts_shp@data[["id"]]))
# Learned later that this problem is caused by reading the id field as a string not
# an integer; such is now remedied by includion of the integer64 option above.
# censustracts <- rasterize(x=censustracts_shp,y=censustracts,field=ids,fun='last')
censustracts <- rasterize(x=censustracts_shp,y=censustracts,field='id',fun='last')
writeRaster(censustracts, 
            filename=filepath("processing/censustract_ids.tif"), format="GTiff", 
            overwrite=TRUE)    

# paxperbuilding is census tract-specific; we have the correspondence in a csv file
# below line allows enough memspace for our 1000+ census tracts
options(max.print = 9999)

ppb = read.table("inputs/censustract_paxperbuilding.csv", 
                 colClasses=c("integer","integer","numeric"), 
                 header=TRUE, sep=",", row.names=1)

# reclassify censustracts to their ppb numeric value
# towards a raster that holds building counts
# create an empty raster for buildings
censustracts_housing <- raster()
censustracts_housing <- setExtent(censustracts_housing, dimensions)
res(censustracts_housing) <- res(country_template)
crs(censustracts_housing) <- crs(country_template)
censustracts_housing <- reclassify(censustracts,ppb,right=FALSE)
writeRaster(censustracts_housing, 
            filename=filepath("processing/censustract_housingfactor.tif"), format="GTiff", 
            overwrite=TRUE)     
            
# towards a raster that holds populationcounts
# create an empty raster 
population <- raster()
population <- setExtent(population, dimensions)
res(population) <- res(country_template)
crs(population) <- crs(country_template)   

population <- round(censustracts_housing*buildings,digits=0)
writeRaster(population, filename=filepath("outputs/population.tif"), format="GTiff", 
            overwrite=TRUE)
spplot(population,main='Population of Bhutan')
        
# Verify overall correctness of population assignation to raster cells
# For Bhutan, total population censused is 601,301.
# Observe that the below shows some discrepancy with the above number.
# This is caused by edged effects op popwogs, and above not being fully robust method.
# We compute #houses and thus ppb per popwog in the database on the basis of
# popwog vector polygons, however, above we have rasterized polygons for popwogs
# and those are not geometrically identical.  Since popwog boundaries often run along
# topographic features, such boundaries may have substantial numbers of houses ...
# See below for improvements
cellStats(population,stat='sum')


##########################################################################################
# Determine statistics and mapplot
# Read the breaks from traveltime_zones csv file:
tt_zones = read.table("inputs/traveltime_zonebreaks.csv", 
                 colClasses=c("integer","integer","character"), 
                 header=TRUE, sep=",", row.names=1)
breaks <- tt_zones$minutes
colors <- tt_zones$color
plotxnames <- 1:(length(breaks)-1)

# determine traveltime zones spatially by breaks given
travtime_zonation <- cut(accCost_min,breaks)
# save as raster and plot for pdf export
writeRaster(travtime_zonation, 
            filename=filepath("outputs/travtime_zonation.tif"), format="GTiff", 
            overwrite=TRUE)
spplot(travtime_zonation,main='Travel time to hospitals zonation')

# determine traveltime zonal statistics
zonestats <- zonal(population,travtime_zonation,fun='sum')

barplot(t(zonestats)[2,],main='Population at travel time to reach hospital',
        names.arg=plotxnames, xlab='Travel time [x 30min]',
        ylab='Population', col=colors)

spplot(accCost_min,main='Travel time to nearest hospital')

##########################################################################################
# Determine statistics and mapplot for one hospital shed, namely #23

popzone <- raster()
popzone <- setExtent(popzone, dimensions)
res(popzone)  <- res(country_template)
crs(popzone) <- crs(popzone)

popzone <- population
popzone <- mask(popzone,landcover_reprojected)
popzone[accCost_zones!=18] <- NA
popzone <- trim(popzone)
# create ttz to have same extent as popzone
ttz <- crop(travtime_zonation,popzone)
popzonestats <- zonal(popzone,ttz,fun='sum')

# Obtain ed popdata above for one hospital shed; now plto the map and the barplot

barplot(t(popzonestats)[2,],main='Hospital shed #18 population at travel time to reach hospital',
       names.arg=plotxnames, xlab='Travel time [x 30min]',
       ylab='Population', col=colors)
       
spplot(popzone, main='Population in hospital shed #18', xlab='Density for 150x150 m.')

# As well as the at-risk areas
poptravzone <- raster()
poptravzone <- setExtent(poptravzone, dimensions)
res(poptravzone)  <- res(country_template)
crs(poptravzone)  <- crs(country_template)
poptravzone <- travtime_zonation
poptravzone[accCost_zones!=18] <- NA
poptravzone <- trim(poptravzone)
spplot(poptravzone, main='Populations at risk in hospital shed #23')

##########################################################################################
# the barplots for all service areas:
# should still be made independent of constants such as 17

zero17 <- matrix(c(1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,0,9,0,10,0,11,0,12,0,13,0,
                   14,0,15,0,16,0,17,0), byrow=TRUE, nrow=17, ncol=2)
for(i in 1:length(targets)) {
  popzone <- raster()
  popzone <- setExtent(popzone, dimensions)
  res(popzone)  <- res(country_template)
  crs(popzone) <- crs(popzone)

  popzone <- population
  popzone <- mask(popzone,landcover_reprojected)
  popzone[accCost_zones!=i] <- NA
  popzone <- trim(popzone)
  # create ttz to have same extent as popzone
  ttz <- crop(travtime_zonation,popzone)
  # determine those stats per ttime zone
  popzonestats <- zonal(popzone,ttz,fun='sum')
  # pad results with zeroes to length 17 if needed
  if ( length(popzonestats[,1])>=17)
    { pzs <- popzonestats }
  else 
    { pzs <- rbind(popzonestats,zero17[(length(popzonestats[,1])+1):17,]) }
  barplot(t(pzs)[2,],
       main=paste('Hospital shed #',i,'population at travel time to reach hospital'),
       names.arg=plotxnames, xlab='Travel time [x 30min]',
       ylab='Population', col=colors)
}

##########################################################################################
## DONE
##########################################################################################
#  Validation;  Using travel times to the capital (proxy used is hospital there)
thimphu_layer <- accCost_stack@layers[[22]]
writeRaster(thimphu_layer, 
            filename=filepath("outputs/thimphu_travelcost.tif"), format="GTiff", 
            overwrite=TRUE)
            
##########################################################################################
#
#  POPULATION DISCREPANCY ANALYSIS
#

popzonestats <- matrix(c(0,0),byrow=TRUE, nrow=1, ncol=2)

for(i in censustracts_shp@data[["id"]]) { 
  popzone <- raster()
  popzone <- setExtent(popzone, dimensions)
  res(popzone)  <- res(country_template)
  crs(popzone) <- crs(popzone)

  popzone <- population
  popzone <- mask(popzone,landcover_reprojected)
  popzone[censustracts!=i] <- NA
  popzonestats <- rbind(popzonestats, c(i,cellStats(popzone,stat='sum')))
}
  
write.csv(popzonestats, "outputs/popzonestats.csv")

##########################################################################################
#
#  Same but now counting houses per census tract:
#

popzonestats <- matrix(c(0,0),byrow=TRUE, nrow=1, ncol=2)
buildings <- mask(buildings,censustracts)

for(i in censustracts_shp@data[["id"]]) { 
  popzone <- raster()
  popzone <- setExtent(popzone, dimensions)
  res(popzone)  <- res(country_template)
  crs(popzone) <- crs(popzone)

  popzone <- buildings
  popzone <- mask(popzone,censustracts)
  popzone[censustracts!=i] <- NA
  popzonestats <- rbind(popzonestats, c(i,cellStats(popzone,stat='sum')))
}
  
write.csv(popzonestats, "outputs/popzonebuildings.csv")
