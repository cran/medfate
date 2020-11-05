## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(medfate)

## -----------------------------------------------------------------------------
spar = defaultSoilParams(2)
print(spar)

## -----------------------------------------------------------------------------
examplesoil = soil(spar)
class(examplesoil)

## -----------------------------------------------------------------------------
names(examplesoil)

## -----------------------------------------------------------------------------
examplesoil$W

## -----------------------------------------------------------------------------
print(examplesoil, model = "SX")

## -----------------------------------------------------------------------------
print(examplesoil, model="VG")

## ---- fig = TRUE, fig.width= 5, fig.height=3, fig.align= 'center', echo=TRUE----
soil_retentionCurvePlot(examplesoil, model="both")

## -----------------------------------------------------------------------------
data("SpParamsMED")

## -----------------------------------------------------------------------------
names(SpParamsMED)

## -----------------------------------------------------------------------------
data(exampleforestMED)
exampleforestMED

## -----------------------------------------------------------------------------
above = forest2aboveground(exampleforestMED, SpParamsMED)
above

## ---- fig = TRUE, fig.width= 4, fig.height=3, fig.align= 'center', echo=TRUE----
vprofile_leafAreaDensity(above, byCohorts = F)

## ---- fig = TRUE, fig.width= 5, fig.height=3, fig.align= 'center', echo=TRUE----
vprofile_leafAreaDensity(above, byCohorts = T)

## -----------------------------------------------------------------------------
Z50 = c(exampleforestMED$treeData$Z50, exampleforestMED$shrubData$Z50)
Z95 = c(exampleforestMED$treeData$Z95, exampleforestMED$shrubData$Z95)

## ---- fig = TRUE, fig.width= 5, fig.height=3, fig.align= 'center', echo=TRUE----
vprofile_rootDistribution(exampleforestMED, SpParamsMED)

## -----------------------------------------------------------------------------
data(examplemeteo)
head(examplemeteo)

## -----------------------------------------------------------------------------
control = defaultControl()
control

## -----------------------------------------------------------------------------
x = spwbInput(above, Z50, Z95, examplesoil, SpParamsMED, control)

## -----------------------------------------------------------------------------
names(x)

## -----------------------------------------------------------------------------
x = forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control)

## -----------------------------------------------------------------------------
d = 100
sd1<-spwb_day(x, examplesoil, rownames(examplemeteo)[d],  
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
examplesoil$W

## -----------------------------------------------------------------------------
resetInputs(x, examplesoil)
examplesoil$W

## -----------------------------------------------------------------------------
S = spwb(x, examplesoil, examplemeteo, latitude = 41.82592, elevation = 100)

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
plot(S, type="PlantTranspiration")
plot(S, type="PlantPhotosynthesis")
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
spwb_waterUseEfficiency(S, type = "Stand An/E", freq = "months", draw=FALSE)

