### R code from vignette source 'RootSystems.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: RootSystems.Rnw:21-23
###################################################
options(width=67)
library(medfate)


###################################################
### code chunk number 2: RootSystems.Rnw:37-42
###################################################
P = root.ldrProfile(c(100,500), c(200,1500), rep(50,31))
matplot(x=t(P), y=seq(0,1500, by=50), type="l", xlim=c(0,0.5), ylim=c(1500,0), ylab = "Depth (mm)", 
        xlab="Root density")
legend("bottomright", legend=c("Z50 = 100, Z95 = 200", 
       "Z50 = 200, Z95 = 1500"), col=c("black","red"), lty=1:2, bty="n")


