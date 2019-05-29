## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=TRUE, include=FALSE-------------------------------------------
library(medfate)

## ------------------------------------------------------------------------
data(exampleforest)
exampleforest

## ------------------------------------------------------------------------
spar = defaultSoilParams(2)
print(spar)

## ------------------------------------------------------------------------
examplesoil = soil(spar)
names(examplesoil)

