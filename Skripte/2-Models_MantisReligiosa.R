#### SDM auf Basis von GBIF und BIOCLIM für ...
#### Mantis religiosa
#### Modellerstellung
 rm(list = ls())

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

### 1- Daten einlesen + Parameter setzen ####
## 1.1- Datensätze einlesen ####
## Artdaten
specDataReg <- readRDS("data/gbif/Mantis_religiosa_europe_dataset.rds")

## Umweltdaten
enviData <- readRDS("data/bioclim/bioclim_europe_Var_bio_Res_5.rds")

## 1.2- Daten extrahieren ####
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



#### 2- Modellerstellung ####

## 2.1- Mit allen Variablen ####
meAll <- maxent(enviData, presences_region@coords, background)

# me1Jack <- maxent(enviData, presences_region@coords, background, args= "-J") ### mit Jacknife, Dauer ca. 10min


## 2.2- Variablen selektieren ####

# 2.2.1- Automatisiert ####
## Umweltvariablen nach Contribution sortieren ##
ContAll <- meAll@results[7:(6+length(names(enviData))),1]
ContAll <- sort(ContAll, decreasing = TRUE)

## Umweltvariablen auswählen deren Bedeutung >= dem 75% Quartil sind ##
ContAll75Quant <- ContAll[ContAll>quantile(ContAll)[4]]

ContImp <- as.character(strsplit(names(ContAll75Quant), ".contribution"))

## Datensatz mit Bedeutenden Umweltvariablen erstellen ##
# Alle
enviImpAll <- enviData[[ContImp]]

# Einzelne
enviImpNamesSing <- c(1:length(ContImp))
for(i in 1:length(ContImp)){
  enviImpNamesSing[i] <- paste("enviImpS", i, sep= "")
  
  assign(enviImpNamesSing[i], enviData[[ContImp[[i]]]])
}

# Ansteigend
enviImpNamesAsc <- c(1:length(ContImp))
j <- 1
for(i in 1:length(ContImp)){
  enviImpNamesAsc[i] <- paste("enviImpA", 1, ":", i, sep= "")
  assign(enviImpNamesAsc[i], enviData[[ContImp[1:i]]])
  j<- j+1
}

rm(ContAll, ContAll75Quant)
# Automatisiert - Ende ###

# 2.2.2- Händisch ####
envi_1 <- enviData[[c("bio3", "bio15", "bio2")]]      ## bio4 und bio6 entfert: Korrelieren beide mit bio3,
envi_2 <- enviData[[c(ContImp, "bio18")]]             ## bio18 hinzugenommen weil es bei einzelnem fehlen den größten negativeffekt zeigt



### 2.3- Modelle aus selektierten Variablen erstellen ####
# 2.3.1- Automatisiert ####
meImpAll <- maxent(enviImpAll, presences_region@coords, background)

# Einzelne
meImpS1 <- maxent(brick(get(enviImpNamesSing[1])), presences_region@coords, background)  ## bio3
meImpS2 <- maxent(brick(get(enviImpNamesSing[2])), presences_region@coords, background)  ## bio15
meImpS3 <- maxent(brick(get(enviImpNamesSing[3])), presences_region@coords, background)  ## bio4
meImpS4 <- maxent(brick(get(enviImpNamesSing[4])), presences_region@coords, background)  ## bio2
meImpS5 <- maxent(brick(get(enviImpNamesSing[5])), presences_region@coords, background)  ## bio6

meImpS <- list(meImpS1, meImpS2, meImpS3, meImpS4, meImpS5)

# Aufsteigende
meImpA1 <- maxent(brick(get(enviImpNamesAsc[1])), presences_region@coords, background)  ## bis bio3 == meImp1
meImpA2 <- maxent(get(enviImpNamesAsc[2]), presences_region@coords, background)  ## bis bio15
meImpA3 <- maxent(get(enviImpNamesAsc[3]), presences_region@coords, background)  ## bis bio4
meImpA4 <- maxent(get(enviImpNamesAsc[4]), presences_region@coords, background)  ## bis bio2
meImpA5 <- maxent(get(enviImpNamesAsc[5]), presences_region@coords, background)  ## bis bio6 == meImpAll

meImpA <- list(meImpA1, meImpA2, meImpA3, meImpA4, meImpA5)

# 2.3.2- Händische ####
# Predictoren wählen ##
envi_1 <- enviData[[c("bio3", "bio15", "bio2")]]      ## bio4 und bio6 entfert: Korrelieren beide mit bio3,
envi_2 <- enviData[[c(ContImp, "bio18")]]             ## bio18 hinzugenommen weil es bei einzelnem fehlen den größten negativeffekt zeigt

# Modelle erstellen
me_1 <- maxent(envi_1, presences_region@coords, background)
me_2 <- maxent(envi_2, presences_region@coords, background)

meM <- list(me_1, me_2)

rm(envi_1, envi_2)

#### 3- MaxEnt-Modelle abspeichern ####

maxMods <- list(meImpS, meImpA, meM, meAll)

saveRDS(maxMods, paste("models/maxent/", species, "_", "europe", "_MaxEntModels.rds", sep= ""))


#### 4- Korrelationen der Parameter ####
### 4.1- Daten an Presence- und Absence-Punkten auswählen ###
presvals <- extract(enviData, as.data.frame(presences_region))   
absvals <- extract(enviData, background)
head(presvals)
head(absvals)


### 4.2- Datensatz der Predictoren erstellen ###
pb <- c(rep(1, nrow(presvals)), rep(0, nrow(absvals)))
sdmdata <- data.frame(cbind(pb, rbind(presvals, absvals)))
sdmdata <- na.omit(sdmdata)
head(sdmdata)                   ## head: gibt die ersten 6 Zeilen wieder
tail(sdmdata)                   ## tail: gibt die letzten 6 Zeilen wieder  
na.omit(sdmdata)

### 4.3- Korrelation zwischen verwendeteten Prediktoren ###
# Scatterplot
png(paste("data/correlations/scatter", species, "europe.png", sep= "_"), pointsize = 16)
pairs(sdmdata[,c(ContImp)], cex=0.1, fig=TRUE)
dev.off()
# Contribution
png(paste("data/correlations/corplot", species, "europe.png", sep= "_"), pointsize = 16)
corrplot(cor(sdmdata[,c(ContImp)]), type= "lower", diag=FALSE)
dev.off()


#### Ende ####
start_time; date()	## Start; und Endzeit abfragen ## Dauer: ca. 5 Minuten
























### 5- Climate Change Projektion ####
cc <- getData('CMIP5', var="bio", res=5, rcp=85, model='HD', year=70, download=TRUE, path="data")
cc <- crop(cc, enviData)
cc <- mask(cc, enviData)
names(cc) <- names(enviData)

pred_cc <- predict(me, cc)

distr_cc <- pred_cc
distr_cc[distr_cc[] < thr$sensitivity] <- NA

### Create figure for distribution for current and climate change projection ####
x11(width=12)
pdf("figures/maps.pdf", width=12, pointsize = 16)
par(mfrow=c(1,2))
plot(distr, col=mypalette(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(region, add=T)
mtext("Current climate ", 3,-2.2)
mtext(species, 3,-1.2, font=3)
mtext("WGS84 ", 1,-1.2, adj=1)
plot(distr_cc, col=mypalette(length(cInt$brks)), breaks=cInt$brks, legend=F, xlab="Longitude", ylab="Latitude")
plot(region, add=T)
mtext("Climate scenario RCP85 HD ", 3,-2.2)
mtext(species, 3,-1.2, font=3)
mtext("WGS84 ", 1,-1.2, adj=1)
dev.off()
