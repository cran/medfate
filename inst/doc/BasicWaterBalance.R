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
x = forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control)

## -----------------------------------------------------------------------------
names(x)

## -----------------------------------------------------------------------------
d = 100
sd1<-spwb_day(x, rownames(examplemeteo)[d],  
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, 
             slope= 0, aspect = 0, prec = examplemeteo$Precipitation[d])

## -----------------------------------------------------------------------------
names(sd1)

## -----------------------------------------------------------------------------
sd1

## -----------------------------------------------------------------------------
x$soil$W

## -----------------------------------------------------------------------------
resetInputs(x)
x$soil$W

## -----------------------------------------------------------------------------
S = spwb(x, examplemeteo, latitude = 41.82592, elevation = 100)

## -----------------------------------------------------------------------------
class(S)

## -----------------------------------------------------------------------------
names(S)

## -----------------------------------------------------------------------------
names(S$spwbInput)

## -----------------------------------------------------------------------------
head(S$WaterBalance)

## -----------------------------------------------------------------------------
head(S$Plants$PlantPsi)

## ---- fig=TRUE, fig.align="center", fig.width=5, fig.height = 3.5-------------
plot(S, type = "PET_Precipitation")
plot(S, type = "Snow")
plot(S, type = "Export")

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 3.5-------------
plot(S, type = "Evapotranspiration")

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 4---------------
plot(S, type="SoilTheta")
plot(S, type="SoilRWC")
plot(S, type="SoilPsi")
plot(S, type="SoilVol")

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 4---------------
plot(S, type="Transpiration")
plot(S, type="GrossPhotosynthesis")
plot(S, type="PlantPsi")
plot(S, type="PlantStress")

## -----------------------------------------------------------------------------
summary(S, freq="months",FUN=mean, output="Soil")

## -----------------------------------------------------------------------------
summary(S, freq="months",FUN=mean, output="PlantStress")

## -----------------------------------------------------------------------------
head(summary(S, freq="day", output="PlantStress", bySpecies = TRUE))

## -----------------------------------------------------------------------------
summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)

## -----------------------------------------------------------------------------
spwb_stress(S, index = "NDD", freq = "years", draw=FALSE)
spwb_stress(S, index = "MDS", freq = "years", draw=FALSE)

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 4---------------
spwb_stress(S, index = "WSI", freq = "months", draw=TRUE)

## -----------------------------------------------------------------------------
spwb_waterUseEfficiency(S, type = "Stand Ag/E", freq = "months", draw=FALSE)

