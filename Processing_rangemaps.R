install.packages('stringr', dependencies = T)
install.packages("plyr", dependencies = T)
install.packages('dplyr', dependencies = T)
install.packages('sf', dependencies = T)
install.packages('sp', dependencies = T)
install.packages('raster', dependencies = T)
install.packages('maptools', dependencies = T)
install.packages('spData', dependencies = T)

library(sf)
library(sp)
library(raster)
library(spData)
#library(spDataLarge)
library(stringr)
library(plyr)
library(dplyr)
library(maptools)
library(rgdal)

# main directory
#setwd('C:/Users/Quinten/Desktop/MASTER REPTILE DIVERSITY PROJECT/R_directory') Q Home
setwd('C:\\Users\\10534598\\DIRECT')  ### Q GIS LAB


rawdata = readOGR(file.path('rawdata/raw_rep/modeled_reptiles.shp'))
rawdata$Genus = word(rawdata$Binomial, 1) #adding a genus column to the raw data file, for my own sake

# using turtledata as a relatively small group to test formulas
turtledata = rawdata[rawdata$Group == 'turtle', ]
#turtledata[1:50,] %>% plot()

# this formula plots species of a certain group on a grid for visual appeal (DO NOT USE FOR GLOBAL GRID, since the grid fits on imput data)
RangeToGrid <- function(x, res,map = FALSE, ...) {
  rastertemp = raster(extent(x), resolution = res, crs = st_crs(x)$proj4string)
  rasterval = rasterize(x, rastertemp, field = 1, fun = 'count')
  rasterplot = plot(rasterval, ...)
  if(map) rasterplot = plot(world$geom, add= T)
  return(rasterplot)
}
#so this is the quick and dirty way to plot reptile species heatmaps
RangeToGrid(turtledata, 1, map = T)

### but the goal is to work with the raster cell polygons

land_grid <- readOGR(file.path('rawdata/shapefiles/360x114global_Land.shp'))
land_grid
#creating subsets to try the data on

Echis <- rawdata[rawdata$Genus == "Echis", ]
head(Echis)
Echis

vermiculatus <- rawdata[rawdata$Binomial == "Anolis vermiculatus",]
maura <- rawdata[rawdata$Binomial == "Natrix maura",]
Natrix <- rawdata[rawdata$Genus == "Natrix",]
writeOGR(vermiculatus, 'subsets', driver="ESRI Shapefile", layer = 'vermiculatus')
writeOGR(maura, 'subsets', driver="ESRI Shapefile", layer = 'maura')
writeOGR(Natrix, 'subsets', driver="ESRI Shapefile", layer = 'Natrix')
writeOGR(Echis, 'subsets', driver="ESRI Shapefile", layer = 'Echis')
writeOGR(turtledata, 'subsets', driver="ESRI Shapefile", layer = 'allturtles')
# the new files  eg. "vermiculatus.shp' contain range polygons of species or genera
# The function to overlap land_grid and other polygons or raster will be callen TDWGoverlap
TDWGoverlap <- function(x, dir, tdwg, type){ 
  #overlaps species range polygons or rasters with twdg 1 degree polygons
  #x = input file  
  #dir = directory where x is stored
  # tdwg = polygonfile representing tdwg degree polygons
  # type = "polygon" or " raster" depending on input data
  #returns dataframe of overlapped cells with cell ID
  # polygon inputs are range maps // raster imputs are rasterfiles with exact same extent & res as tdwg polygon  
  if(type == "polygon") { 
    polygon <- readOGR(file.path(dir,x))
    overlapList <- over(polygon, tdwg, returnList = TRUE)
    names(overlapList) <- 1:length(overlapList)
    overlapDF <- ldply(overlapList, .id = "polygon")
    overlapDF$SpecName <- gsub(x, pattern = ".shp ", replacement = "") 
  } 
  else {
    sp_raster <- raster(file.path(dir,x))
    sp_polys <- rasterToPolygons(sp_raster, fun = function(x){x>0})
    overlapList <- over(sp_polys, tdwg, returnList = TRUE)
    names(overlapList) <- 1:length(overlapList)
    overlapDF <- ldply(overlapList, .id = "gridcell") # note the grid cells are unique within species but not across species
    overlapDF$SpecName <- gsub(x, pattern = "\\.tif", replacement = "")  }  
  return(overlapDF)
} 
a <- TDWGoverlap(x='maura.shp', dir='subsets', tdwg=land_grid, type = "polygon")
b <- TDWGoverlap(x ="vermiculatus.shp", dir = 'subsets', tdwg = land_grid, type = "polygon")
c <- TDWGoverlap(x= "Natrix.shp", dir = 'subsets', tdwg = land_grid, type = "polygon")
d <- TDWGoverlap(x='Echis.shp', dir = 'subsets', tdwg = land_grid, type = 'polygon')
e <- TDWGoverlap(x='allturtles.shp', dir = 'subsets', land_grid = tdwg, type = 'polygon')

# these are data subsets to try functions on, 2 single species, 2 genera and the whole turtle dataset
# the next thing we need is to transform the outputdfs into usefull data
# I do this by sorting them by cell ID and sum of unique species
# so in this practice
quickCO <- function(x) {
 # quick count overlaps
 # x is a TDWGoverlap df, returns a dataframe of species count by HBWID cell
  A1 <- x[ , c("HBWID", "SpecName")]
  A1$SpecName <- 1
  names(A1)[2] <- "Spcount"
  A2 <- aggregate(.~ HBWID, data = A1, FUN = sum)
  return(A2)
}
# The quick Count Overlaps function turns a TDWIG overlap df into a df with number of species per HBWID 
# only usable in this script i guess
A <- quickCO(a)
B <- quickCO(b)
C <- quickCO(c)
D <- quickCO(d)
E <- quickCO(e)
hist(E$Spcount)
E
# now we test whether this worked, by comparing a range map with the overlap gridcells
#lets try 3 species : Anolis proboscis (pinnochio anole), Varanus salvator (water monitor)
# and the Echis carinatus (saw scale viper) nasty little thing
salvator <- rawdata[rawdata$Binomial == "Varanus salvator",]
plot(world$geom)
plot(salvator, add = T, col = "red")
## occurs in south east asia, sounds about right
#now we try to recreate this with gridcells
writeOGR(salvator, 'subsets', driver="ESRI Shapefile", layer = 'salvator')
salvagrid <- TDWGoverlap(x='salvator.shp', dir='subsets', tdwg=land_grid, type = "polygon")
occ_a <- land_grid[land_grid$HBWID %in% salvagrid$HBWID, ] #this line subsets the gridcells on the gridcells occurring in salvagrid
plot(world$geom)
plot(occ_a, add = T, col = "red")
## looks like the range to gricell procedure and the subsetting worked for the monitor
## lets also try with the other 2 species
proboscis <- rawdata[rawdata$Binomial == "Anolis proboscis",]
plot(world$geom)
plot(proboscis, add = T, col = "red")
plot(proboscis, add = F)
writeOGR(proboscis, 'subsets', driver="ESRI Shapefile", layer = 'proboscis')
pinocchiogrid <- TDWGoverlap(x='proboscis.shp', dir='subsets', tdwg=land_grid, type = "polygon")
occ_b <- land_grid[land_grid$HBWID %in% pinocchiogrid$HBWID, ] #this line subsets the gridcells on the gridcells occurring in salvagrid
plot(world$geom)
plot(occ_b, add = T, col = "green")
# this doesnt work yet for some reason, get to it alter
carinatus <- rawdata[rawdata$Binomial == "Echis carinatus",]
plot(world$geom)
plot(carinatus, add = T, col = "gold3")
writeOGR(carinatus, 'subsets', driver="ESRI Shapefile", layer = 'carinatus')
carinatusgrid <- TDWGoverlap(x='carinatus.shp', dir='subsets', tdwg=land_grid, type = "polygon")
occ_c <- land_grid[land_grid$HBWID %in% carinatusgrid$HBWID, ] #this line subsets the gridcells on the gridcells occurring in salvagrid
plot(world$geom)
plot(occ_c, add = T, col = "gold3")
# this works again
# the pinnochio anole grid doesnt work, the ragne polygon is a very small circle but it sould overlap with 
# a singe gridcell, so this needs to be fixed

# next up, lets see if we can do this for manys species at the same time, ENTRY OF THE TURTLES!
# the allturtles file is already created and ready for usen (we called it "e" earlier)
occ_turtles <- land_grid[land_grid$HBWID %in% e$HBWID, ]
plot(world$geom)
plot(turtledata, add = T, col = "blue") #this costs a lot of processing power, since its calculates many rangemaps
plot(world$geom)
plot(occ_turtles, add = T, col = "blue") # this is way faster because the many rangemaps are boiled down to gridcell polygons
###### these occ_ dataframes are fake, they are just subsetted gridcells and do not contain info, we need to include
###### sprichness into them and see if this can be plotted, this im trying with turtles as well
###### by trying to incorporate the occ_ frames with the quickCO data. 
E
occ_turtles
?merge
turtlesp <- merge(occ_turtles, E, by = "HBWID")
turtlesp
names(occ_turtles)
names(turtlesp)
hist(turtlesp$Spcount, 20)
dev.off()
plot(turtlesp, col = ifelse(turtlesp$Spcount == 1, "gray90",
             ifelse(turtlesp$Spcount <4 , "gray75", ifelse(turtlesp$Spcount <6, "gray60",
  ifelse(turtlesp$Spcount <9 , "gray45", ifelse(turtlesp$Spcount <13, "gray30", "gray15")) ))), border = 'transparent')
plot(world$geom, add = T)
#this seems to work, now i can color cells with cerain amounts of species or a range of richness values. 
# a way to test if this really works is making a global sp richness map like this and comparing it to the 
# maps produced by Roll et al.  (i think) 
# making a global richness map like this might take some computation time but lets go 
# sidenote, im not doing all reptiles but only snakes, turtles and lizards, so this is only to test the code.


#mokergrid <- TDWGoverlap(x='modeled_reptiles.shp', dir='rawdata/raw_rep', tdwg=tdwg, type = "polygon")
# ERROR: cannot allocate vector of size 2.0 GB     (ohboi)
# than we have to do the groups separately and compare them to some other roll et al figures. 
lizdata <- rawdata[rawdata$Group == "lizard", ]
writeOGR(lizdata, 'subsets', driver="ESRI Shapefile", layer = 'allliz')
overlapliz <- TDWGoverlap(x='allliz.shp', dir='subsets', tdwg=land_grid, type = "polygon")
lizgrids <- land_grid[land_grid$HBWID %in% overlapliz$HBWID, ]
lizrich <- quickCO(overlapliz)
lizmerge <- merge(lizgrids, lizrich, by = "HBWID")
# we can make a histogram of sprichness to see how to arrange the colors on the map
hist(lizrich$Spcount)
plot(world$geom, add = T)
plot(lizmerge, border = "transparent", col = ifelse(lizmerge$Spcount <3 , "gray95", ifelse(lizmerge$Spcount < 6, "gray80", 
     ifelse(lizmerge$Spcount < 12, "gray70", ifelse(lizmerge$Spcount < 20, "gray60",
      ifelse(lizmerge$Spcount < 26, "gray50", ifelse(lizmerge$Spcount < 35, "gray40", 
      ifelse(lizmerge$Spcount < 42, "gray30", ifelse(lizmerge$Spcount < 55, "gray20" ,"gray10")))))))))
## the color scale has to be changed a bit but this figure can now be compared to roll et al,
## along with the turtle figure and a snake figure if we please.
snakedata <- rawdata[rawdata$Group == "snake", ]
writeOGR(snakedata, 'subsets', driver="ESRI Shapefile", layer = 'allsnake')
overlapsnake <- TDWGoverlap(x='allsnake.shp', dir='subsets', tdwg=land_grid, type = "polygon")
snakegrids <- land_grid[land_grid$HBWID %in% overlapsnake$HBWID, ]
snakerich <- quickCO(overlapsnake)
snakemerge <- merge(snakegrids, snakerich, by = "HBWID")


# next, it is needed to create a full data.frame with all HBWID cells
# this complete data.frame must include spcount and all geodiveristy index values per HBWID cell
# so we must transform the Versteegh (2016) geodiversity index rasterdata into the HBWID data.frame
# thus the biggest challenge, working with the geodiveristy raster
# first we can easily make data.frames for each subclass
lizard.df <- as.data.frame(lizmerge)
turtle.df <- as.data.frame(turtlesp)
snake.df <- as.data.frame(snakemerge)
#bruteforcing the global dataframe
allgrid <- as.data.frame(tdwg)


lDF <- lizard.df[, c("HBWID", "Spcount")]
names(lDF) [2] <- "lizards" 
names(lDF)

sDF <- snake.df[, c("HBWID", "Spcount")]
names(sDF) [2] <- "snakes" 
names(sDF)

tDF <- turtle.df[, c("HBWID", "Spcount")]
names(tDF) [2] <- "turtles" 
names(tDF)
#as of now i dont know the real way to do this quickly

global.df <- merge(allgrid, lDF, by= "HBWID", all = T)
global.df <- merge(global.df, sDF, by = 'HBWID', all = T)
global.df <- merge(global.df, tDF, by = 'HBWID', all = T)
global.df[is.na(global.df)] <- 0 
global.df$reptileSP <- rowSums(x= cbind(global.df$turtles, global.df$snakes, global.df$lizards))
global.df
# now we have created a global data frame containing all gridcells and reptile richness values. 
# The raw range maps are now processed, the next thing I need to do is process geodiversity data!
