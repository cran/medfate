## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE, include=FALSE-----------------------------------------------
library(medfate)

## -----------------------------------------------------------------------------
data("SpParamsMED")

## -----------------------------------------------------------------------------
names(SpParamsMED)

## -----------------------------------------------------------------------------
data("SpParamsDefinition")
knitr::kable(SpParamsDefinition[,-2])

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
spar = defaultSoilParams(2)
print(spar)

## -----------------------------------------------------------------------------
examplesoil = soil(spar)
class(soil)

## -----------------------------------------------------------------------------
names(examplesoil)

## -----------------------------------------------------------------------------
print(examplesoil, model = "SX")

## -----------------------------------------------------------------------------
examplesoil$W

## -----------------------------------------------------------------------------
examplesoil$Temp

## -----------------------------------------------------------------------------
print(examplesoil, model="VG")

## ---- fig = TRUE, fig.width= 5, fig.height=3, fig.align= 'center', echo=TRUE----
soil_retentionCurvePlot(examplesoil, model="both")

## -----------------------------------------------------------------------------
data(examplemeteo)
head(examplemeteo)

## -----------------------------------------------------------------------------
control = defaultControl()
names(control)

