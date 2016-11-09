### R code from vignette source 'FireBehaviour.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: FireBehaviour.Rnw:23-27
###################################################
options(width=67)
library(medfate)
data("SpParamsMED")
data(exampleforest)


###################################################
### code chunk number 2: FireBehaviour.Rnw:179-186
###################################################
vprofile.FuelBulkDensity(exampleforest, SpParamsMED, z = seq(0,1000, by=1))
fs = fuel.stratification(exampleforest,SpParamsMED, bulkDensityThreshold = 0.04)
abline(v=0.04, col="red", lty=2)
abline(h=fs$surfaceLayerBaseHeight, col="gray", lty=3, lwd=1.5)
abline(h=fs$surfaceLayerTopHeight, col="gray", lty=3, lwd=1.5)
abline(h=fs$canopyBaseHeight, col="gray", lty=3, lwd=1.5)
abline(h=fs$canopyTopHeight, col="gray", lty=3, lwd=1.5)


