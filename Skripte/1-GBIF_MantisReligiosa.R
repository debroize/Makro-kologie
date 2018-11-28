#### GBIF und BIOCLIM Daten Vorbereitung zu SDM für ...
#### Mantis religiosa
rm(list= ls())

#### Start ####
start_time <- date()	# Startzeit speichern

### 0.1- Workspace setzen(Pfad in Zwischenablage) ####
 wdname <- "C:/Users/Denis/Documents/Makrooekologie/Workspace"
# wdname <- gsub( "\\\\",  "/",  readClipboard())
setwd(wdname); getwd(); rm(wdname)


### 0.2- Notwendige Ordnerstruktur in Working Directory ####
if("data" %in% list.files() == FALSE){dir.create("data/")}                    ## Daten Ordner
if("gbif" %in% list.files("data/") == FALSE){dir.create("data/gbif/")}        ## Unterordner für GBIF-Daten
if("bioclim" %in% list.files("data/") == FALSE){dir.create("data/bioclim/")}  ## Unterordner für GBIF-Daten


### 0.3- Pakete laden  ####
library(rgbif)                  ## Global Biodiversity Information Facility, Datenbank für Artvorkommen
library(raster)                 ## Rasterverarbeitung und Bioclim-Daten
library(dismo)                  ## Für MaxEnt-Modellierung
library(maptools)               ## Für Weltkarten-Polygone
library(colorRamps)
library(classInt)



###### 1- Datengenerierung #######

### 1.0- Parameter eintragen ####
speciesName <- 'Mantis religiosa'; species <-  gsub(" ","_", speciesName);    ## Art die beobachtet werden soll

  bioclimVarAll <- c('tmin', 'tmax', 'prec', 'bio')
  bioclimResAll <- c(0.5, 2.5, 5, 10)
bioclimVar <- bioclimVarAll[4]                                ## Variablen für Bioclim
bioclimRes <- bioclimResAll[3]                                ## Auflösung der Bioclim-Daten (0.5 nur mit Argumenten: lon, lat)
  
regNum <- 150; regnName <- "europe"            ## Region Nummer aus BioClim
regDel <- "Russia"         ## Name von Ländern die auseschlossen werden sollen

rm(bioclimVarAll, bioclimResAll)  

### 1.1- Rohdaten von GBIF, BioClim und Weltkarte laden, bzw.  downloaden und speichern ####
### 1.1.1- Daten von GBIF laden und speichern, bzw. laden ##
if(paste(species,".rds", sep= "") %in% list.files("data/gbif/")){       ## Überprüfen ob Daten schon vorliegen, falls JA -> direkt laden
 
   presences <- readRDS(paste("data/gbif/",species, ".rds", sep= ""))  

}else{
  
  key <- name_suggest(q= speciesName, rank='species')$key[1]
  n <- occ_count(taxonKey=key, georeferenced=TRUE); n                   ## Anzahl Datenpunkte überprüfen
  presences <- occ_search(taxonKey=key, limit=n)
  head(presences$data)
  
  presences <- as.data.frame(presences$data[, c("decimalLatitude", "decimalLongitude")])
  presences <- na.omit(presences)

  ## Speichern
  saveRDS(presences, paste("data/gbif/",species, ".rds", sep= ""))
}
rm(speciesName)

### 1.1.2- Bioclim Daten downoaden ##
  ## Parameter
  bioclimResFolder <- paste("wc", bioclimRes, "/", sep= "")
  bioclimNameZIP <- (paste(bioclimVar, bioclimRes, "bil.zip", sep= "_" ))
  
dlbioclim <- TRUE
if(bioclimNameZIP %in% list.files(paste("data/", bioclimResFolder, sep= ""))){dlbioclim <- FALSE}
bioclim <- getData('worldclim', download=dlbioclim, path="data/bioclim/", var= bioclimVar, res= bioclimRes) ## Einlesen bzw. Downloaden von worldclim Daten.  

  projection(bioclim)

rm(bioclimResFolder, bioclimNameZIP, dlbioclim)

### 1.1.3- Weltkarte laden ##
data(wrld_simpl)  
projection(wrld_simpl)



### 1.2- Rohdaten auf relevante Bereiche einschränken ####
### 1.2.1- Projektion der Daten gleichsetzen ##
projection(wrld_simpl) <- projection(bioclim)     ## Welkarte auf bioclim beziehen


### 1.2.2- Weltkarte auf Region beschränken ##
region <- wrld_simpl[wrld_simpl$REGION==regNum,]
  #region@data 
region <- region[region$NAME!=regDel,]

  #plot(region)
rm(regNum, regDel)

### 1.2.3- Punktdaten räumlich verorten ##
presences <- SpatialPoints(presences[,c("decimalLongitude", "decimalLatitude")], proj4string=CRS(projection(bioclim)))

  #plot(wrld_simpl)
  #points(decimalLatitude ~ decimalLongitude, data=presences, pch=16, cex=0.5, col="red")


### 1.2.4- Datenpunkte auf Region beschränken ##
  #presences_region <- presences[is.na(over(region, presences)),]  ### Funktioniert noch nicht (daten werden nict auf europa begrenzt)
presences_region <- presences[region]

  #plot(region)
  #points(presences_region, pch=16, cex=0.5, col="blue")

### 1.2.5- BioClim-Daten auf Region beschränken ##
bioclim_region <- crop(bioclim, region)
bioclim_region <- mask(bioclim_region, region)

rm(bioclim)  
### 1.2.6- Background-Datenpunkte ertellen ##
background <- randomPoints(bioclim_region, 50000)

  #plot(wrld_simpl)
  #points(background, pch=16, cex=0.5, col="blue")

### 1.2.7- Überprüfen ob Datenpunkte nur auf Region beschränkt sind
x11()
plot(wrld_simpl)
points(background, pch=16, cex=0.5, col="blue")
points(decimalLatitude ~ decimalLongitude, data=presences, pch=16, cex=0.5, col="red")
points(presences_region, pch=16, cex=0.5, col="green")

rm(presences)

### 1.3- Datensatz erstellen und abspeichern ####
rawData <- list(species, region, presences_region, background)
  ## Artdaten
  saveRDS(rawData, paste("data/gbif/", species, "_", regnName, "_dataset.rds", sep= ""))
  ## Umweltdaten
  saveRDS(bioclim_region, paste("data/bioclim/bioclim_", regnName, "_Var_", bioclimVar, "_Res_", bioclimRes, ".rds", sep = ""))


#### Ende ####
start_time; date()	## Start; und Endzeit abfragen ## Dauer: ca. 10 Minuten