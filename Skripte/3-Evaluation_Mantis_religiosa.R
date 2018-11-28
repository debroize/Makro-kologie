#### SDM auf Basis von GBIF und BIOCLIM f?r ...
#### Mantis religiosa
#### Modellbewertung
rm(list = ls())
#28.11.
#### Start ####
start_time <- date()	# Startzeit speichern

### 0.1- Workspace setzen(Pfad in Zwischenablage) ####
 wdname <- "C:/Users/Denis/Documents/Makrooekologie/Workspace"
# wdname <- gsub( "\\\\",  "/",  readClipboard())
setwd(wdname); getwd(); rm(wdname)


### 0.2- Notwendige Ordnerstruktur in Working Directory ####
if("data" %in% list.files() == FALSE){dir.create("data/")}                      ## Daten Ordner
if("gbif" %in% list.files("data/") == FALSE){dir.create("data/gbif/")}          ## Unterordner für GBIF-Daten
if("bioclim" %in% list.files("data/") == FALSE){dir.create("data/bioclim/")}    ## Unterordner für Bioclim-Daten
if("figures" %in% list.files() == FALSE){dir.create("figures/")}                ## Ordner für Grafiken
if("models" %in% list.files() == FALSE){dir.create("models/")}                  ## Modell Ordner
if("maxent" %in% list.files("models/") == FALSE){dir.create("models/maxent/")}  ## Unterordner für MaxEnt-Modelle


### 0.3- Pakete laden  ####
library(rgbif)                  ## Global Biodiversity Information Facility, Datenbank f?r Artvorkommen
library(raster)                 ## Rasterverarbeitung und Bioclim-Daten
library(dismo)                  ## Für MaxEnt-Modellierung
library(maptools)               ## Für Weltkarten-Polygone
library(colorRamps)             
library(classInt)
library(rJava)                  ## Java Implementierung
library(MaxentVariableSelection)
library(corrplot)               
library(rgdal)
library(gbm)
library(hier.part)



### 1- Daten einlesen  ####
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



#### 2- Contribution ####
#x11()
png(paste("figures/Cont_", species, "_", "europe", ".png", sep= ""), width = 1200, height = 800, res = 120)
par(mfrow=c(3,5))
for(i in 1:length(me)){
  plot(me[[i]], main = paste("AUC", me[[i]]@results[5]))
}
dev.off()



### 3- Evaluation ####
## 3.1- Auswählen welche Modelle betrachtet werden sollen ####
meAll <- me[[length(me)]]

meSing <- me[[1]]
  for(i in 1:5){    ## Wählt besten einzelnen Parameter aus
   if(me[[i]]@results[5] > meSing@results[5]){meSing <- me[[i]]}
  }  
meSing

meImp <- me[[10]]

me1   <- me[[11]]


## 3.2- Evaluieren ####
eAll <- evaluate(presences_region, background, meSing, enviData)
eSing <- evaluate(presences_region, background, meAll, enviData)
eImp <- evaluate(presences_region, background, meImp, enviData)

e1 <- evaluate(presences_region, background, me1, enviData) ## Gutes Modell mit nur 3 Predictoren (untereinander kaum korreliert)

eAll
eSing
eImp
e1

### 4- Betrachten der Evaluation ####
## 4.1- Thresholds ##
thrAll <- threshold(eAll)
thrSing <- threshold(eSing)
thrImp <- threshold(eImp)
thr1 <- threshold(e1)


## 4.2- ROC-Kurven und Density-Plots ##
x11()
par(mfrow=c(1,2)) 
plot(eAll, "ROC", sub= "eAll")
density(eAll)

x11()
par(mfrow=c(1,2)) 
plot(eSing, "ROC", sub= "eSing")
density(eSing)

x11()
par(mfrow=c(1,2))
plot(eImp, "ROC", sub= "eImp")
density(eImp)

x11()
par(mfrow=c(1,2))
plot(e1, "ROC", sub= "e1")
density(e1)



#### 5- Response Kurven ####
### 5.1- 1D Response Kurven###
x11()
response(meAll)
response(meImp)
response(me1)



### 5.2- 2D Response Kurven (in Bearbeitung)###
np <- 30
newdata <- expand.grid(bio10=seq(145, 200, len=np), bio18=seq(0, 240, len=np))
newdata$pred <- predict(me1, newdata)

## 3.4.1 Use threshold to show distribution
newdata$pred[newdata$pred<thr$sensitivity] <- NA

## 3.4.2- Create classes of site suitability
cInt <- classIntervals((newdata$pred))

xdiff <-diff(unique(newdata$bio10))[1]
ydiff <-diff(unique(newdata$bio18))[1]

mypalette <- colorRampPalette(c("lightgreen", "darkgreen"))
newdata$colors <- findColours(cInt, mypalette(length(cInt$brks)))

par(mfrow=c(1,1), mar=c(5,5,1,1))
symbols(x=newdata$bio10, y=newdata$bio18, rectangles=matrix(rep(c(xdiff, ydiff), nrow(newdata)), ncol=2, byrow=T), bg=newdata$colors, fg="white", inches=F, xlab="Temperature of warmest quarter (°dC)", ylab="Precipitation of warmest quarter (mm)")
contour(x=unique(newdata$bio10), y=unique(newdata$bio18), z=matrix(newdata$pred, nrow=np), add=T, levels=unique(round(cInt$brks,1)), labcex = 1.3)
mtext(species, side=3, line=-1.3, font=3)
mtext(paste0("AUC = " , round(e@auc, 2), " "), side=1, line=-2.3, adj=1)
mtext(paste0("Pearson r = " , round(e@cor, 2), " "), side=1, line=-1.3, adj=1)



#### Ende ####
start_time; date()	## Start; und Endzeit abfragen ## Dauer: ca. 1 Min