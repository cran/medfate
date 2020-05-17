## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(medfate)

## -----------------------------------------------------------------------------
spar = defaultSoilParams(2)
print(spar)

## -----------------------------------------------------------------------------
examplesoil = soil(spar, VG_PTF = "Toth")
names(examplesoil)

## -----------------------------------------------------------------------------
examplesoil$W

## -----------------------------------------------------------------------------
examplesoil$Temp

## -----------------------------------------------------------------------------
print(examplesoil, model = "SX")

## -----------------------------------------------------------------------------
print(examplesoil, model="VG")

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

## -----------------------------------------------------------------------------
Z = c(exampleforestMED$treeData$Z95, exampleforestMED$shrubData$Z)
Z
below = forest2belowground(exampleforestMED, examplesoil, SpParamsMED)
below

## -----------------------------------------------------------------------------
root_ldrDistribution(exampleforestMED$treeData$Z50[1], 
                     exampleforestMED$treeData$Z95[1],
                     examplesoil$dVec)

## -----------------------------------------------------------------------------
data(examplemeteo)
head(examplemeteo)

## -----------------------------------------------------------------------------
control = defaultControl()

## -----------------------------------------------------------------------------
names(control)

## -----------------------------------------------------------------------------
x = growthInput(above, Z, below, examplesoil, SpParamsMED, control)

## -----------------------------------------------------------------------------
x = forest2growthInput(exampleforestMED, examplesoil, SpParamsMED, control)

## -----------------------------------------------------------------------------
names(x)

## -----------------------------------------------------------------------------
x$cohorts

## -----------------------------------------------------------------------------
x$above

## -----------------------------------------------------------------------------
x$paramsTransp

## -----------------------------------------------------------------------------
x$below

## -----------------------------------------------------------------------------
x$paramsGrowth
x$paramsAllometries

## -----------------------------------------------------------------------------
G1<-growth(x, examplesoil, examplemeteo, latitude = 41.82592, elevation = 100)


## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(G1, "GrossPhotosynthesis", bySpecies = T)

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(G1, "MaintenanceRespiration", bySpecies = T)
plot(G1, "GrowthRespiration", bySpecies = T)

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(G1, "CarbonBalance", bySpecies = T)

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(G1, "SugarLeaf", bySpecies = T)
plot(G1, "StarchLeaf", bySpecies = T)

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(G1, "SugarSapwood", bySpecies = T)
plot(G1, "StarchSapwood", bySpecies = T)

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(G1, "LeafArea", bySpecies = T)
plot(G1, "SapwoodArea", bySpecies = T)

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(G1, "LAgrowth", bySpecies = T)
plot(G1, "SAgrowth", bySpecies = T)

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
G1$StandStructures

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
G1$StandStructures

