## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(medfate)

## ------------------------------------------------------------------------
spar = defaultSoilParams(2)
print(spar)

## ------------------------------------------------------------------------
examplesoil = soil(spar, VG_PTF = "Toth")
names(examplesoil)

## ------------------------------------------------------------------------
examplesoil$W

## ------------------------------------------------------------------------
examplesoil$Temp

## ------------------------------------------------------------------------
print(examplesoil, model = "SX")

## ------------------------------------------------------------------------
print(examplesoil, model="VG")

## ------------------------------------------------------------------------
data("SpParamsMED")

## ------------------------------------------------------------------------
names(SpParamsMED)

## ------------------------------------------------------------------------
data(exampleforest)
exampleforest

## ------------------------------------------------------------------------
above = forest2aboveground(exampleforest, SpParamsMED)
above

## ------------------------------------------------------------------------
Z = c(exampleforest$treeData$Z95, exampleforest$shrubData$Z)
Z
below = forest2belowground(exampleforest, examplesoil, SpParamsMED)
below

## ------------------------------------------------------------------------
root_ldrDistribution(exampleforest$treeData$Z50[1], 
                     exampleforest$treeData$Z95[1],
                     examplesoil$dVec)

## ------------------------------------------------------------------------
data(examplemeteo)
head(examplemeteo)

## ------------------------------------------------------------------------
control = defaultControl()

## ------------------------------------------------------------------------
names(control)

## ------------------------------------------------------------------------
control$storagePool = "two"

## ------------------------------------------------------------------------
x = growthInput(above, Z, below, examplesoil, SpParamsMED, control)

## ------------------------------------------------------------------------
x = forest2growthInput(exampleforest, examplesoil, SpParamsMED, control)

## ------------------------------------------------------------------------
names(x)

## ------------------------------------------------------------------------
x$cohorts

## ------------------------------------------------------------------------
x$canopy

## ------------------------------------------------------------------------
x$above

## ------------------------------------------------------------------------
x$paramsTransp

## ------------------------------------------------------------------------
x$below

## ------------------------------------------------------------------------
x$paramsGrowth
x$paramsAllometries

## ------------------------------------------------------------------------
G1<-growth(x, examplesoil, examplemeteo, elevation = 100)


## ---- fig=TRUE, echo=FALSE, fig.width=7, fig.height=4--------------------
plot(G1, "PlantPhotosynthesisPerLeaf")
plot(G1, "PlantTranspirationPerLeaf")

## ---- fig=TRUE, echo=FALSE, fig.width=7, fig.height=4--------------------
plot(G1, "PlantRespiration")
plot(G1, "PlantRespirationPerLeaf")
plot(G1, "PlantRespirationPerIndividual")
plot(G1, "PlantCBalance")
plot(G1, "PlantCBalancePerLeaf")
plot(G1, "PlantCBalancePerIndividual")

## ---- fig=TRUE, echo=FALSE, fig.width=7, fig.height=4--------------------
plot(G1, "PlantSAgrowth")
plot(G1, "PlantRelativeSAgrowth")

## ---- fig=TRUE, echo=FALSE, fig.width=7, fig.height=4--------------------
plot(G1, "PlantLAIlive")
plot(G1, "PlantLAIdead")

## ---- fig=TRUE, echo=FALSE, fig.width=7.5, fig.height=3.5----------------
plot(G1, "PlantCstorageFast")
plot(G1, "PlantCstorageSlow")

