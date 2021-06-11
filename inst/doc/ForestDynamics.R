## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(medfate)

## -----------------------------------------------------------------------------
spar = defaultSoilParams(2)

## -----------------------------------------------------------------------------
examplesoil = soil(spar)

## -----------------------------------------------------------------------------
data(exampleforestMED)
exampleforestMED

## -----------------------------------------------------------------------------
data(examplemeteo)
head(examplemeteo)

## -----------------------------------------------------------------------------
data("SpParamsMED")

## -----------------------------------------------------------------------------
control = defaultControl("Granier")

## -----------------------------------------------------------------------------
meteo2001 = examplemeteo
meteo2002 = examplemeteo
meteo2002$Precipitation = meteo2002$Precipitation/2
meteo2003 = examplemeteo
row.names(meteo2002) = seq(as.Date("2002-01-01"), 
                           as.Date("2002-12-31"), by="day")
row.names(meteo2003) = seq(as.Date("2003-01-01"), 
                           as.Date("2003-12-31"), by="day")
meteo_01_03 = rbind(meteo2001, meteo2002, meteo2003)

## -----------------------------------------------------------------------------
fd<-fordyn(exampleforestMED, examplesoil, SpParamsMED, meteo_01_03, control, 
           latitude = 41.82592, elevation = 100)

## -----------------------------------------------------------------------------
fd$StandSummary

## -----------------------------------------------------------------------------
fd$SpeciesSummary

## -----------------------------------------------------------------------------
fd$TreeTable

## -----------------------------------------------------------------------------
fd$DeadTreeTable

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(fd$GrowthResults[[1]], "LeafArea", bySpecies = T)
plot(fd$GrowthResults[[2]], "LeafArea", bySpecies = T)
plot(fd$GrowthResults[[3]], "LeafArea", bySpecies = T)

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(fd$GrowthResults[[1]], "SapwoodArea", bySpecies = T)
plot(fd$GrowthResults[[2]], "SapwoodArea", bySpecies = T)
plot(fd$GrowthResults[[3]], "SapwoodArea", bySpecies = T)

