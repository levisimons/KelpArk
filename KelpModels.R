rm(list=ls())
require(data.table)
require(sf)
require(spThin)
require(raster)
require(stars)
require(terra)
require(dplyr)
require(plyr)
require(dismo)
require(randomForest)
require(DescTools)
require(ggplot2)
require(viridis)
require(virtualspecies)

#Set working directory
wd <- ""
setwd(wd)

#Set random number string
set.seed(1)

#Specific taxa found in the Pacific or Gulf of Mexico study area
Pacific_taxa <- c("mastocarpus","gracilaria_andersonii","hedophyllum_sessile","postelsia_palmaeformis","laminaria_setchellii","fucus_distichus","nereocystis_luetkeana","macrocystis_pyrifera","alaria_marginata","saccharina_latissima","pyropia")
Gulf_taxa <- c("sargassum","ulva","codium","gracilaria","eucheuma")
#Designate taxon
taxon <- "sargassum"

#Intake species occurrence data from GBIF into a data table.  Use the variable name taxon in this process.


#Convert species data to a spatial points object


#Spatially thin occurrence data so that no two points are closer than 9.26 km, the resolution of Bio-Oracle data.


#Get lists of all current and future environmental rasters in tif format.


#Build a raster stack of all environmental rasters.
#Rasters are generated using https://github.com/levisimons/KelpArk/blob/main/KelpRasters.R

#Update raster stack names for current and future environmental raster stacks so they match the environmental raster file names

#Filter collinear environmental variables using current environmental raster stack.
#Use a 0.75 Spearman correlation threshold, and 1000 background points, as per https://onlinelibrary.wiley.com/doi/pdf/10.1002/ece3.10901


#Build a raster stack of all current environmental rasters using just the filtered layers.

#Build a raster stack of all future environmental rasters using just the filtered layers.

#Update raster stack names for current and future environmental raster stacks so they match the environmental raster file names

#Extract current environmental raster values at species occurrence points


#Remove empty rows from the extracted environmental data


#Add presence column to extracted environmental data and set its value to 1.


#Set presence column to be a factor from being numeric for modeling.


#Set Pacific area boundaries
Pacific <- st_bbox(c(xmin=-179,xmax=-117,ymin=32.5,ymax=61.5))
#Set Gulf of Mexico area boundaries
Gulf <- st_bbox(c(xmin=-99,xmax=-81.5,ymin=18,ymax=31))
#Select study extent based on selected taxon.
if(taxon %in% Pacific_taxa){points_buffer <- st_as_sfc(st_bbox(Pacific))}
if(taxon %in% Gulf_taxa){points_buffer <- st_as_sfc(st_bbox(Gulf))}

#Determine the number of rows in the extracted environmental data

#Generate a set of random points, three times the size of the number of rows in the extracted environmental data, within points_buffer


#Convert single column coordinates to standard longitude/latitude columns
background_points <- sf::st_coordinates(background_points)

#Convert background points object to a data table
background_points <- as.data.table(background_points)

#Convert from longitude (X) and latitude (Y) columns to sf object.
background_points <- sf::st_as_sf(background_points,coords = c("X","Y"),  crs = 4326)

#Extract raster values at background points
background_extracted <- raster::extract(env_rasters, background_points)

#Update column names so the column names match the environmental raster file names
colnames(background_extracted) <- env_layers

#Remove empty rows
background_extracted <- as.data.frame(background_extracted[complete.cases(background_extracted),])

#Add presence column
background_extracted$presence <- 0

#Set presence variable to factor for modeling.
background_extracted$presence <- as.factor(background_extracted$presence)

#Set a predictable random number generator seed for reproducibility.
set.seed(1)

#Create an empty list to store prediction rasters.
raster_predict_list <- c()
#Create an empty list to store future prediction rasters.
future_raster_predict_list <- c()
#Create an empty list to store relative importance outputs.
importance_list <- c()
#Create an empty list to store accuracy outputs.
accuracy_list <- c()
#Create an empty list to store partial plot outputs.
partial_plot_list <- c()
j <- 1
i_max <- 50
for(i in 1:i_max){
  #Create a subset of the presence/background data with the following properties:
  #1. Composed of a randomly selected 80% of rows from scb_extracted.
  #2. Composed of rows randomly selected from background_extracted. The number of rows will also be 80% of rows found in scb_extracted.
  #3. Merged these two subsets together.
  subset_extracted <- rbind(env_extracted[sample(nrow(env_extracted),0.8*nrow(env_extracted)),],background_extracted[sample(nrow(background_extracted),0.8*nrow(env_extracted)),])
  
  #Run a random forest model over this data subset.
  rf1 <- suppressWarnings(tuneRF(x=subset_extracted[,!(colnames(subset_extracted) %in% "presence")],y=subset_extracted$presence,stepFactor=1,plot=FALSE,doBest=TRUE))
  
  #Make a prediction raster from the random forest model and store it in a list.
  raster_predict_list[[i]] <- dismo::predict(env_rasters,rf1,progress='text')
  #Plot predicted raster
  plot(raster_predict_list[[i]],col=viridis(10))
  
  #Make a future prediction raster from the random forest model and store it in a list.
  future_raster_predict_list[[i]] <- dismo::predict(future_env_rasters,rf1,progress='text')
  #Plot predicted raster
  plot(future_raster_predict_list[[i]],col=viridis(10))
  
  #Store relative importance of variable outputs as a temporary data frame.
  tmp <- as.data.frame(rf1$importance)
  #Set one column to store the variable names from the row names.
  tmp$VariableName <- rownames(tmp)
  #Store this importance data frame in the importance list.
  importance_list[[i]] <- tmp
  
  #Calculate the true skill statistic TSS to evaluate model accuracy.
  sensitivity <- rf1$confusion[[1]] / (rf1$confusion[[1]]+rf1$confusion[[2]])
  specificity <- rf1$confusion[[4]] / (rf1$confusion[[4]]+rf1$confusion[[3]])
  TSS <- sensitivity+specificity-1
  #Store TSS results
  accuracy_list[i] <- TSS
  
  #Loop through each environmental variable and store the partial response outputs in a temporary data frame.
  for(env_layer in env_retain){
    #Store partial plot chart data in a temporary data frame.
    tmp <- as.data.frame(partialPlot(rf1,subset_extracted[,!(colnames(subset_extracted) %in% "presence")],x.var=c(env_layer),plot=F))
    #Transform logistic probabilities to regular probabilities.
    tmp$y <- exp(tmp$y) / (1+exp(tmp$y))
    #Rename probability column
    colnames(tmp) <- c(env_layer,"Detection Probability")
    #Store partial plot data in a list of data frames.
    partial_plot_list[[j]] <- tmp
    j <- j+1
  }

  print(paste(i,Sys.time()))
}

#Stack the list of prediction rasters.
raster_predict <- brick(raster_predict_list)
#Sum the list of prediction rasters.
raster_predict <- calc(raster_predict, sum)
#Save raster output.
writeRaster(raster_predict,paste(taxon,"_prediction.tif",sep=""),overwrite=T)

# Convert full raster to data frame
raster_df <- as.data.frame(raster_predict, xy = TRUE)
names(raster_df)[3] <- "value"

# Extract non-zero points
raster_points <- subset(raster_df, value > 0.5*i_max)

#Plot prediction raster
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "lightblue", high = "lightblue", na.value = "grey") +
  guides(fill = "none") +
  geom_point(data = raster_points, aes(x = x, y = y, color = value), size = 1) +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = paste("Predicted occurrences of ",taxon," (2020)",sep=""),
       x="Longitude degrees East",y = "Latitude degrees North",
       color = paste("Predicted frequency\nout of ",i_max," models",sep=""))

#Get geographic range of predicted taxon occurrences
range(na.omit(raster_df[raster_df$value > 0.5*i_max,"x"]))
range(na.omit(raster_df[raster_df$value > 0.5*i_max,"y"]))
#Count the number of locations predicted to have suitable habitat.
nrow(raster_df[raster_df$value > 0.5*i_max,])

#Stack the list of future prediction rasters.
future_raster_predict <- brick(future_raster_predict_list)
#Sum the list of future prediction rasters.
future_raster_predict <- calc(future_raster_predict, sum)
#Save raster output.
writeRaster(raster_predict,paste(taxon,"_future_prediction.tif",sep=""),overwrite=T)

# Convert full raster to data frame
raster_df <- as.data.frame(future_raster_predict, xy = TRUE)
names(raster_df)[3] <- "value"

# Extract non-zero points
raster_points <- subset(raster_df, value > 0.5*i_max)

#Plot prediction raster
ggplot() +
  geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
  scale_fill_gradient(low = "lightblue", high = "lightblue", na.value = "grey") +
  guides(fill = "none") +
  geom_point(data = raster_points, aes(x = x, y = y, color = value), size = 1) +
  scale_color_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = paste("Predicted occurrences of ",taxon," (2060)",sep=""),
       x="Longitude degrees East",y = "Latitude degrees North",
       color = paste("Predicted frequency\nout of ",i_max," models",sep=""))

#Get geographic range of predicted taxon occurrences
range(na.omit(raster_df[raster_df$value > 0.5*i_max,"x"]))
range(na.omit(raster_df[raster_df$value > 0.5*i_max,"y"]))
#Count the number of locations predicted to have suitable habitat.
nrow(raster_df[raster_df$value > 0.5*i_max,])

#Calculate the mean TSS for the models
mean(accuracy_list)
sd(accuracy_list)

#Convert list of importance data frames to a single data frame.
importance_total <- rbind.fill(importance_list)
#Calculate the mean relative importance for each variable.
importance_total <- aggregate(x=importance_total$MeanDecreaseGini,by = list(importance_total$VariableName),FUN = mean)
#Rename columns.
colnames(importance_total) <- c("VariableName","Importance")
#Convert importance to rank importance.
importance_total$Importance <- rank(desc(importance_total$Importance))
#Save rank importance table.
write.table(importance_total,paste(taxon,"_rank_importance.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)

#Collapse partial plot outputs into single data frame.
partial_plots <- rbind.fill(partial_plot_list)
partial_plots <- as.data.frame(partial_plots)
write.table(partial_plots,paste(taxon,"_partial_plots.txt",sep=""),quote=FALSE,sep="\t",row.names = FALSE)
partial_plots <- read.table(paste(taxon,"_partial_plots.txt",sep=""), header=TRUE, sep="\t",as.is=T,skip=0,fill=TRUE,check.names=FALSE, encoding = "UTF-8")

#Plot partial dependence heat maps for continuous data.
k <- 7
ggplot(partial_plots, aes(x=!!sym(names(env_rasters[[k]])), y=`Detection Probability`) )+
  xlab(names(env_rasters[[k]]))+ylab("Detection\nProbability")+
  geom_bin2d(bins = 50)+
  scale_fill_continuous(type = "viridis",name=paste("Frequency\n(Out of ",i_max," models)",sep=""))+
  stat_smooth(aes(y = `Detection Probability`, fill=`Detection Probability`),method="auto",formula=y~x,color="violet",fill="red",n=0.1*sum(!is.na(partial_plots[,names(env_rasters[[k]])])))+
  theme_bw(base_size=25)
