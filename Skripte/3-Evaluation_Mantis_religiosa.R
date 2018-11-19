#### SDM auf Basis von GBIF und BIOCLIM f?r ...
#### Mantis religiosa
#### Modellbewertung
rm(list = ls())
### 0.1- Workspace setzen(Pfad in Zwischenablage) ####
 wdname <- "C:/Users/Denis/Documents/Makrooekologie/Workspace"
# wdname <- gsub( "\\\\",  "/",  readClipboard())
setwd(wdname); getwd(); rm(wdname)


### 0.2- Notwendige Ordnerstruktur in Working Directory ####
if("data" %in% list.files() == FALSE){dir.create("data/")}                      ## Daten Ordner
if("gbif" %in% list.files("data/") == FALSE){dir.create("data/gbif/")}          ## Unterordner f?r GBIF-Daten
if("bioclim" %in% list.files("data/") == FALSE){dir.create("data/bioclim/")}    ## Unterordner f?r GBIF-Daten
if("models" %in% list.files() == FALSE){dir.create("models/")}                  ## Modell Ordner
if("maxent" %in% list.files("models/") == FALSE){dir.create("models/maxent/")}  ## Unterordner f?r MaxEnt-Modelle


### 0.3- Pakete laden  ####
library(rgbif)                  ## Global Biodiversity Information Facility, Datenbank f?r Artvorkommen
library(raster)                 ## Rasterverarbeitung und Bioclim-Daten
library(dismo)
library(maptools)               ## Für Weltkarten-Polygone
library(colorRamps)
library(classInt)
library(rJava)
library(MaxentVariableSelection)
library(corrplot)
library(rgdal)
library(gbm)
library(hier.part)
library(maxent)


### 1- Daten einlesen + Parameter setzen ####
## 1.1- Datens?tze einlesen ##
## Artdaten
specDataReg <- readRDS("data/gbif/Mantis_religiosa_europe_dataset.rds")

## Umweltdaten
enviData <- readRDS("data/bioclim/bioclim_europe_Var_bio_Res_5.rds")

## 1.2- Daten extrahieren ##
species <- specDataReg[[1]]
region <- specDataReg[[2]]
presences_region <- specDataReg[[3]]
background <- specDataReg[[4]]

rm(specDataReg)

## 1.3- Umweltvariablen Bioclim ####
##BIO1 = Annual Mean Temperature, ##BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)), ##BIO3 = Isothermality (BIO2/BIO7) (* 100), BIO4 = Temperature Seasonality (standard deviation *100)
##BIO5 = Max Temperature of Warmest Month, ##BIO6 = Min Temperature of Coldest Month, ##BIO7 = Temperature Annual Range (BIO5-BIO6)
##BIO8 = Mean Temperature of Wettest Quarter, ##BIO9 = Mean Temperature of Driest Quarter, ##BIO10 = Mean Temperature of Warmest Quarter
##BIO11 = Mean Temperature of Coldest Quarter, ##BIO12 = Annual Precipitation, ##BIO13 = Precipitation of Wettest Month
##BIO14 = Precipitation of Driest Month, ##BIO15 = Precipitation Seasonality (Coefficient of Variation), ##BIO16 = Precipitation of Wettest Quarter
##BIO17 = Precipitation of Driest Quarter,##BIO18 = Precipitation of Warmest Quarter, ##BIO19 = Precipitation of Coldest Quarter





## 1.4- Modelle einlesen ####
maxMods <- readRDS("models/maxent/Mantis_religiosa_europe_MaxEntModels.rds")
me <- unlist(maxMods)

meS <- maxMods[[1]]
meA <- maxMods[[2]]
meM <- maxMods[[3]]

#### 2- Bewertung ####
### 2.1- Contribution
x11()
par(mfrow=c(3,5))
for(i in 1:length(me)){
  plot(me[[i]], main = paste("AUC", me[[i]]@results[5]))
}




