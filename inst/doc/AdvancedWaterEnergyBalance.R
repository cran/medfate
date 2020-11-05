## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(medfate)

## -----------------------------------------------------------------------------
examplesoil = soil(defaultSoilParams(2))

## -----------------------------------------------------------------------------
examplesoil$Temp

## -----------------------------------------------------------------------------
data("SpParamsMED")

## -----------------------------------------------------------------------------
names(SpParamsMED)

## -----------------------------------------------------------------------------
data(exampleforestMED)
exampleforestMED

## -----------------------------------------------------------------------------
data(examplemeteo)
head(examplemeteo)

## -----------------------------------------------------------------------------
control = defaultControl()

## -----------------------------------------------------------------------------
control$transpirationMode = "Sperry"
control$soilFunctions = "VG"

## -----------------------------------------------------------------------------
x = forest2spwbInput(exampleforestMED, examplesoil, SpParamsMED, control)

## -----------------------------------------------------------------------------
x$cohorts

## -----------------------------------------------------------------------------
x$canopy

## -----------------------------------------------------------------------------
x$above

## -----------------------------------------------------------------------------
x$below

## -----------------------------------------------------------------------------
x$belowLayers

## -----------------------------------------------------------------------------
x$paramsAnatomy

## -----------------------------------------------------------------------------
x$paramsTranspiration

## -----------------------------------------------------------------------------
x$paramsWaterStorage

## ---- fig=TRUE, fig.align="center", fig.width=5, fig.height = 3.5-------------
hydraulics_vulnerabilityCurvePlot(x, type="leaf")
hydraulics_vulnerabilityCurvePlot(x, type="stem")
hydraulics_vulnerabilityCurvePlot(x, type="root")
hydraulics_vulnerabilityCurvePlot(x, examplesoil, type="rhizo")

## ---- fig=TRUE, fig.align="center", fig.width=5, fig.height = 3.5-------------
hydraulics_supplyFunctionPlot(x, examplesoil, type="E")
hydraulics_supplyFunctionPlot(x, examplesoil, type="ERhizo")
hydraulics_supplyFunctionPlot(x, examplesoil, type="dEdP")
hydraulics_supplyFunctionPlot(x, examplesoil, type="StemPsi")

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 3.5-------------
d = 100
transp_stomatalRegulationPlot(x, examplesoil, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="E")
transp_stomatalRegulationPlot(x, examplesoil, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="An")
transp_stomatalRegulationPlot(x, examplesoil, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="Gw")
transp_stomatalRegulationPlot(x, examplesoil, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="T")
transp_stomatalRegulationPlot(x, examplesoil, examplemeteo, day = d, timestep=12,
                              latitude = 41.82592, elevation = 100, type="VPD")

## ---- fig=TRUE, fig.align="center", fig.width=5, fig.height = 3.5-------------
moisture_pressureVolumeCurvePlot(x, segment="leaf", fraction="symplastic")
moisture_pressureVolumeCurvePlot(x, segment="leaf", fraction="apoplastic")
moisture_pressureVolumeCurvePlot(x, segment="stem", fraction="symplastic")
moisture_pressureVolumeCurvePlot(x, segment="stem", fraction="apoplastic")

## -----------------------------------------------------------------------------
sd1<-spwb_day(x, examplesoil, rownames(examplemeteo)[d],  
             examplemeteo$MinTemperature[d], examplemeteo$MaxTemperature[d], 
             examplemeteo$MinRelativeHumidity[d], examplemeteo$MaxRelativeHumidity[d], 
             examplemeteo$Radiation[d], examplemeteo$WindSpeed[d], 
             latitude = 41.82592, elevation = 100, 
             slope= 0, aspect = 0, prec = examplemeteo$Precipitation[d])

## -----------------------------------------------------------------------------
names(sd1)

## -----------------------------------------------------------------------------
sd1$WaterBalance

## -----------------------------------------------------------------------------
sd1$Soil

## -----------------------------------------------------------------------------
names(sd1$EnergyBalance)

## ---- fig=TRUE, fig.align="center", fig.width=5, fig.height = 4---------------
plot(sd1, type = "Temperature")
plot(sd1, type = "CanopyEnergyBalance")
plot(sd1, type = "SoilEnergyBalance")

## -----------------------------------------------------------------------------
sd1$Plants

## ---- fig=TRUE, fig.align="center", fig.width=5, fig.height = 4---------------
plot(sd1, type = "PlantTranspiration", bySpecies = T)
plot(sd1, type = "TranspirationPerLeaf", bySpecies = T)
plot(sd1, type = "NetPhotosynthesis", bySpecies = T)
plot(sd1, type = "LeafPsiAverage", bySpecies = T)

## -----------------------------------------------------------------------------
sd1$SunlitLeaves
sd1$ShadeLeaves

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 4---------------
plot(sd1, type = "LeafTemperature", bySpecies=TRUE)

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 4---------------
plot(sd1, type = "LeafGrossPhotosynthesis", bySpecies=TRUE)

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 4---------------
plot(sd1, type = "LeafNetPhotosynthesis", bySpecies=TRUE)

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 4---------------
plot(sd1, type = "LeafStomatalConductance", bySpecies=TRUE)

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 4---------------
plot(sd1, type = "LeafVPD", bySpecies=TRUE)

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 4---------------
plot(sd1, type = "LeafPsi", bySpecies=TRUE)

## ---- fig=TRUE, fig.align="center", fig.width=7.5, fig.height = 4-------------
plot(sd1, type = "LeafCi", bySpecies=TRUE)

## ---- fig=TRUE, fig.align="center", fig.width=7.5, fig.height = 4-------------
plot(sd1, type = "LeafIntrinsicWUE", bySpecies=TRUE)

## -----------------------------------------------------------------------------
examplesoil$W

## -----------------------------------------------------------------------------
examplesoil$Temp

## -----------------------------------------------------------------------------
x$canopy

## -----------------------------------------------------------------------------
resetInputs(x, examplesoil)
examplesoil$W
examplesoil$Temp
x$canopy

## -----------------------------------------------------------------------------
S = spwb(x, examplesoil, examplemeteo[110:170,], latitude = 41.82592, elevation = 100)

## -----------------------------------------------------------------------------
names(S)

## -----------------------------------------------------------------------------
names(S$spwbInput)

## -----------------------------------------------------------------------------
head(S$WaterBalance)

## -----------------------------------------------------------------------------
head(S$Plants$LeafPsiMin)

## ---- fig=TRUE, fig.align="center", fig.width=7, fig.height = 3.5-------------
plot(S, type="Evapotranspiration")
plot(S, type="SoilPsi", bySpecies = TRUE)
plot(S, type="PlantTranspiration", bySpecies = TRUE)
plot(S, type="LeafPsiMin", bySpecies = TRUE)

## -----------------------------------------------------------------------------
summary(S, freq="months",FUN=mean, output="Soil")

## -----------------------------------------------------------------------------
summary(S, freq="months",FUN=mean, output="PlantStress")

## -----------------------------------------------------------------------------
head(summary(S, freq="day", output="PlantStress", bySpecies = TRUE))

## -----------------------------------------------------------------------------
summary(S, freq="month", FUN = mean, output="PlantStress", bySpecies = TRUE)

