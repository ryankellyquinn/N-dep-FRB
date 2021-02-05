
#Load libraries necessary to run script
#Define functions necessary to run script 

library(readr) #NOW
library(splitstackshape)
library(daymetr)
library(dplyr)
library(MALDIquant)
library(AICcmodavg)
library(tidyr)

se.error <- function(x)sd(x)/sqrt(length(x))


#Load datasets 

fs <- read_csv("data/field-sites.csv") #lat/lon, temp., etc. data for NEON field sites  
ndep <- read_csv("data/gbc20824-sup-0001-2018gb005990-ds01.csv") #simulated N deposition data from Ackerman et al. 
fred <- read_csv("data/FRED2_20180518.csv") #Root biomass data and data corresponding to collection location from FRED database (Fine Root Ecology Database- https://roots.ornl.gov/) 
neon.rb <- read_csv("data/NEON_rootsample.csv") #lat/lon, temp., etc. data for NEON field sites  
neon.cn <- read_csv("data/NEON_CN.csv") #lat/lon, temp., etc. data for NEON field sites  


#Combine neon.rb and neon.cn using the unique sample ID for each root sample. 
#Calculate the percentage of C and N in each root sample, and then the g of C and N in each sample. 

neon.rb$carbonpercent <- ifelse(is.na(match(paste(neon.rb$sampleID), paste(neon.cn$cnSampleID))),
                                0, neon.cn$carbonPercent)
neon.rb$nitrogenPercent <- ifelse(is.na(match(paste(neon.rb$sampleID), paste(neon.cn$cnSampleID))),
                                  0, neon.cn$nitrogenPercent)
neon.rb$cfrac <- neon.rb$carbonpercent/100
neon.rb$nfrac <- neon.rb$nitrogenPercent/100
neon.rb$rootc <- neon.rb$incrementRootBiomass*neon.rb$cfrac
neon.rb$rootn <- neon.rb$incrementRootBiomass*neon.rb$nfrac


#Split depthincrementID column (a column that contains the site, profile, and depth increment) into 3 columns

rb <- cSplit(neon.rb, "depthIncrementID", ".")
names(rb)[names(rb) == "depthIncrementID_1"] <- "site"
names(rb)[names(rb) == "depthIncrementID_2"] <- "profile"
names(rb)[names(rb) == "depthIncrementID_3"] <- "depth"


#Aggregate data (but n=1, no avging yet)

m.all <- aggregate(neon.rb$incrementRootBiomass, by = list(neon.rb$depthIncrementID, neon.rb$rootStatus, neon.rb$sizeCategory), "mean")
colnames(m.all) <- c("depthincrementID", "rootstatus", "sizecategory", "rootbiomass")
m.all.n <- aggregate(neon.rb$rootc, by = list(neon.rb$depthIncrementID, neon.rb$rootStatus, neon.rb$sizeCategory), "length")

colnames(m.all)


#Match nitrogen deposition date with NEON sites
#Currently must run 2x. FIX THIS! I think I need to add newlat/ newlon columns before running the loop 

#ndep.glob <- read.csv(folder)
ndep.tot.2014 <- ndep[,c("latitude", "longitude", "pixel_area_km2","tot_2014")]
ndep.global <- list()
for (n in 1:nrow(fs)){
  lat <- fs$Latitude[n]
  lon <- fs$Longitude[n]
  sorted_lat <- sort(ndep.tot.2014$latitude)
  fs$newlat[n] <- sorted_lat[MALDIquant::match.closest(lat, sorted_lat)][1]
  sorted_lon <- sort(ndep.tot.2014$longitude)
  fs$newlon[n] <- sorted_lon[MALDIquant::match.closest(lon, sorted_lon)][1]
  ndep.point <- ndep.tot.2014[which(ndep.tot.2014$latitude==fs$newlat[n] & ndep.tot.2014$longitude==fs$newlon[n]),]
  ndep.global[n] <- ndep.point$tot_2014/100 #convert per-km to per-hectare
}  
output <- unlist(ndep.global)
#names(output) <- fs$`Site ID`
newoutput <- as.data.frame(output)
newoutput$site <- fs$`Site ID` 


#Add MAT, MAP, and nitrogen deposition data stored in fs and newoutput (fs contains MAT, MAP, data specific to each NEON megapit site, newoutput contains simulated N deposition data specific to each NEON megapit site) to corresponding datapoint in rb (remember, rb contains the NEON root biomass and root cn data)   to rb.all
#Split columns where necessary

rb$MAP <- fs$`Mean Annual Precipitation`[match(rb$site, fs$`Site ID`)]
rb$MAT <- fs$`Mean Annual Temperature`[match(rb$site, fs$`Site ID`)]
rb$latitude <- fs$Latitude[match(rb$site, fs$`Site ID`)]
rb$longitude <- fs$Longitude[match(rb$site, fs$`Site ID`)]
rb$type <- fs$`class`[match(rb$site, fs$`Site ID`)]
rb <- cSplit(rb, "MAP", " ")
rb <- cSplit(rb, "MAT", "/")
rb <- cSplit(rb, "type", ",")
#rb <- cSplit(rb, "ty", " ")
rb$MAT_3 <- gsub( "C", "", as.character(rb$MAT_1))
rb$ndep <- newoutput$output[match(rb$site, newoutput$site)]
rb$domain <- fs$`Domain Name`[match(rb$site, newoutput$site)]


#Organize FRED rb data  

fred <- subset(fred, fred$sample_depth == 10) #only want to work with fine root biomass samples of the same volume during colletion, which is 10cm increments. Scaling is not sufficient here because it combines depths and mitigates potential differences between depth increments. 
fred$type_1 <- ifelse(fred$dominant_plant =="tree", "forest", 
                      ifelse(fred$dominant_plant == "graminoid", "grassland",
                             ifelse(fred$dominant_plant == "shrub", "shrub", 
                                    ifelse(fred$dominant_plant == "shrub/tree", "shrub", 
                                           "NA")))) #only interested in working with these dominant plant species, so make any other species NA. Rename the plant types to match NEON data 

fred <- subset(fred, !is.na(fred$type_1)) #remove where we don't know plant type 
fred <- subset(fred, !is.na(fred$Latitude)) #remove where we don't know Latitude
fred <- subset(fred, !is.na(fred$Longitude)) #remove where we don't know Longitude 

fred$MAP_1 <- fred$MAP #????????????WHY DID I DO THIS 
fred$MAT_3 <- fred$MAT #????????????WHY DID I DO THIS


#Add nitrogen deposition data to FRED database FRB data using lat/lon provided by FRED database (remember we deleted any data without a corresponding lat/lon in an earlier chunk of this spreadsheet)

#this loop was made by Zoey Werbin (thank you Zoey!)
#ndep.glob <- read.csv(folder)
ndep.tot.2014 <- ndep[,c("latitude", "longitude", "pixel_area_km2","tot_2014")]
ndep.global <- list()
for (n in 1:nrow(fred)){
  lat <- fred$Latitude[n]
  lon <- fred$Longitude[n]
  sorted_lat <- sort(ndep.tot.2014$latitude)
  fred$newlat[n] <- sorted_lat[MALDIquant::match.closest(lat, sorted_lat)][1]
  sorted_lon <- sort(ndep.tot.2014$longitude)
  fred$newlon[n] <- sorted_lon[MALDIquant::match.closest(lon, sorted_lon)][1]
  ndep.point <- ndep.tot.2014[which(ndep.tot.2014$latitude==fred$newlat[n] & ndep.tot.2014$longitude==fred$newlon[n]),]
  ndep.global[n] <- ndep.point$tot_2014/100 #convert per-km to per-hectare
}  
output <- unlist(ndep.global)
#names(output) <- fs$`Site ID`
newoutput <- as.data.frame(output)
#newoutput$site <- fs$`Site ID` 
fred$ndep <- newoutput[,1]
fred$incrementRootBiomass <- fred$rootbiomass


#Remove dead frb from NEON dataset, we are only interested in live frb 
#Create dataframes for each dataset with matching col names
#Combine the two  dataframes (thus combine datasets)
live.fine <- subset(rb, rb$rootStatus == "live" & rb$sizeCategory == "<=2mm")
rb.neon <- NULL
rb.neon <- data.frame("rb" = live.fine$incrementRootBiomass, 
                      "ndep" = live.fine$ndep, 
                      "map" = live.fine$MAP_1,
                      "mat" = live.fine$MAT_3,
                      "type" = live.fine$type_1,
                      "lat" = live.fine$latitude,
                      "lon" = live.fine$longitude,
                      "depth" = live.fine$depth) 

rb.fred <- NULL
rb.fred <- data.frame("rb" = fred$incrementRootBiomass, 
                      "map" = fred$MAP_1,
                      "mat" = fred$MAT, 
                      "type" = fred$type_1,
                      "lat" = fred$Latitude,
                      "lon" = fred$Longitude,
                      "ndep" = fred$ndep, 
                      "depth" = fred$mindepth)
rb.fred$source <- "fred"
rb.neon$source <- "neon"
df <- rbind(rb.fred, rb.neon)

