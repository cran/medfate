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
x$internalCarbon

## -----------------------------------------------------------------------------
G1<-growth(x, examplemeteo, latitude = 41.82592, elevation = 100)


## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(G1, "GrossPhotosynthesis", bySpecies = T)

## ---- fig=TRUE, echo=TRUE, fig.width=7, fig.height=3--------------------------
plot(G1, "MaintenanceRespiration", bySpecies = T)
plot(G1, "GrowthCosts", bySpecies = T)

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

## ---- eval = FALSE------------------------------------------------------------
#  shinyplot(G1)

