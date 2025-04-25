rm(list=ls())
require(data.table)
require(raster)
require(sf)

#Set working directory
wd <- ""
setwd(wd)

#Set random number string
set.seed(1)

#Set Pacific area boundaries
Pacific <- st_bbox(c(xmin=-179,xmax=-117,ymin=32.5,ymax=61.5))

#Set Gulf of Mexico area boundaries
Gulf <- st_bbox(c(xmin=-99,xmax=-81.5,ymin=18,ymax=31))

#Get a list of all environmental map layers with the file extension .nc (NetCDF)
map_layers <- list.files(path="MapLayers",pattern = "\\.nc$")

#Convert all of the NetCDF map layers from Bio-Oracle into clipped rasters in tif format.
for(map_layer in map_layers){
  #Read in marine layer
  marine_layer <- brick(paste("MapLayers/",map_layer,sep=""))
  #Get base file name
  map_layer_name <- gsub("\\.nc$","", map_layer)
  #Convert from NetCDF to tif format.
  writeRaster(marine_layer, paste("MapLayers/",map_layer_name,".tif",sep=""), bylayer=FALSE,overwrite=TRUE)
  #Read in tiff formatted raster.
  marine_raster <- raster(paste("MapLayers/",map_layer_name,".tif",sep=""))
  #Crop the raster to the Pacific extent.
  pacific_raster <- crop(marine_raster,Pacific)
  #Export the Pacific cropped raster.
  writeRaster(pacific_raster,paste("MapLayers/Pacific_",map_layer_name,".tif",sep=""), bylayer=FALSE,overwrite=TRUE)
  #Crop the raster to the Gulf of Mexico extent.
  gulf_raster <- crop(marine_raster,Gulf)
  #Export the Pacific cropped raster.
  writeRaster(gulf_raster,paste("MapLayers/Gulf_",map_layer_name,".tif",sep=""), bylayer=FALSE,overwrite=TRUE)
}
