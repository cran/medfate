### R code from vignette source 'Hydraulics.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: Hydraulics.Rnw:21-23
###################################################
options(width=67)
library(medfate)


###################################################
### code chunk number 2: Hydraulics.Rnw:42-45
###################################################
kxmax = 8.0
c = 2
d = 3


###################################################
### code chunk number 3: Hydraulics.Rnw:49-53
###################################################
par(mar=c(4,4,1,1))
psiVec = seq(-0.1, -7.0, by =-0.01)
k = unlist(lapply(psiVec, hydraulics.xylemConductance, kxmax, c, d))
plot(-psiVec, k, type="l",ylab="Xylem conductance", xlab="Canopy sap pressure (-MPa)", lwd=1.5)


###################################################
### code chunk number 4: Hydraulics.Rnw:62-77
###################################################
#Maximum soil conductance
#Textural parameters (1 MPa = 0.00009804139432 cm)
#Sandy clay loam 
krmax1 = 80000000000
alpha1 = 602 #602 MPa-1 = 0.059 cm-1   
n1 = 1.48
#Silt loam
krmax2 = 4000000000
alpha2 = 203 #203 MPa-1 = 0.020 cm-1
n2 = 1.41
#Silty clay
krmax3 = 400000000
alpha3 = 50.99887 #50.99887 MPa-1 = 0.005 cm-1
n3 = 1.09
textures = c("Sandy clay loam","Silt loam", "Silty clay")


###################################################
### code chunk number 5: Hydraulics.Rnw:81-90
###################################################
par(mar=c(4,4,1,1))
psiVec = seq(-0.01, -7.0, by =-0.01)
k1 = unlist(lapply(psiVec, hydraulics.vanGenuchtenConductance, krmax1, n1, alpha1))
k2 = unlist(lapply(psiVec, hydraulics.vanGenuchtenConductance, krmax2, n2, alpha2))
k3 = unlist(lapply(psiVec, hydraulics.vanGenuchtenConductance, krmax3, n3, alpha3))
plot(-psiVec, k1, type="l", col="black", ylab="Rhizosphere conductance", xlab = "Soil hydraulic pressure (-MPa)",  xlim=c(0, max(-psiVec)), ylim=c(0, 100))
lines(-psiVec, k2, lty=2, lwd=1.5)
lines(-psiVec, k3, lty=3, lwd=1.5)
legend("topright", bty="n", legend=textures, lty=1:3, cex=0.8, lwd=1.5)


###################################################
### code chunk number 6: Hydraulics.Rnw:103-123
###################################################
par(mar=c(4,4,1,1))
psiSoil1 = -0.1
psiSoil2 = -0.5
psiSoil3 = -2.0
psiSoil4 = -3.0
psiVec = seq(-0.1, -7.0, by =-0.01)
E1 = unlist(lapply(psiVec, hydraulics.EXylem, psiSoil1, kxmax, c, d))
E2 = unlist(lapply(psiVec, hydraulics.EXylem, psiSoil2, kxmax, c, d))
E3 = unlist(lapply(psiVec, hydraulics.EXylem, psiSoil3, kxmax, c, d))
E4 = unlist(lapply(psiVec, hydraulics.EXylem, psiSoil4, kxmax, c, d))
# Ec1 = Ecrit(psiSoil1, kxmax, c, d)
# Ec2 = Ecrit(psiSoil2, kxmax, c, d)
# Ec3 = Ecrit(psiSoil3, kxmax, c, d)
# Ec4 = Ecrit(psiSoil4, kxmax, c, d)
Psic = hydraulics.psiCrit(c,d)
plot(-psiVec, E1, type="l", ylab="Xylem flow rate", xlab="Canopy sap pressure (-MPa)", ylim=c(0,20), lwd=1.5)
lines(-psiVec, E2, lty=1, lwd=1.5)
lines(-psiVec, E3, lty=1, lwd=1.5)
lines(-psiVec, E4, lty=1, lwd=1.5)
# points(rep(-Psic,4), c(Ec1,Ec2, Ec3,Ec4))


###################################################
### code chunk number 7: Hydraulics.Rnw:128-159
###################################################
par(mar=c(4,4,1,1))
Evec = seq(0, 10, by=0.1)
psiSoil1 = -0.01
psiSoil2 = -1.0
psiSoil3 = -2.0
psiSoil4 = -3.0
psiVG11 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil1, krmax1, n1, alpha1))
psiVG12 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil2, krmax1, n1, alpha1))
psiVG13 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil3, krmax1, n1, alpha1))
psiVG14 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil4, krmax1, n1, alpha1))
psiVG21 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil1, krmax2, n2, alpha2))
psiVG22 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil2, krmax2, n2, alpha2))
psiVG23 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil3, krmax2, n2, alpha2))
psiVG24 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil4, krmax2, n2, alpha2))
psiVG31 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil1, krmax3, n3, alpha3))
psiVG32 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil2, krmax3, n3, alpha3))
psiVG33 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil3, krmax3, n3, alpha3))
psiVG34 = unlist(lapply(Evec, hydraulics.E2psiVanGenuchten, psiSoil4, krmax3, n3, alpha3))
plot(-psiVG11, Evec, type="l", col="black", ylab="Rhizosphere flow rate", xlab = "Root sap pressure (-MPa)",  xlim=c(0, 8))
lines(-psiVG21, Evec, col="black", lty=2, lwd=1.5)
lines(-psiVG31, Evec, col="black", lty=3, lwd=1.5)
lines(-psiVG12, Evec, col="black", lty=1, lwd=1.5)
lines(-psiVG22, Evec, col="black", lty=2, lwd=1.5)
lines(-psiVG32, Evec, col="black", lty=3, lwd=1.5)
lines(-psiVG13, Evec, col="black", lty=1, lwd=1.5)
lines(-psiVG23, Evec, col="black", lty=2, lwd=1.5)
lines(-psiVG33, Evec, col="black", lty=3, lwd=1.5)
lines(-psiVG14, Evec, col="black", lty=1, lwd=1.5)
lines(-psiVG24, Evec, col="black", lty=2, lwd=1.5)
lines(-psiVG34, Evec, col="black", lty=3, lwd=1.5)
legend("topright", bty="n", legend=textures, lty=1:3, lwd=1.5, cex=0.8)


###################################################
### code chunk number 8: Hydraulics.Rnw:167-198
###################################################
par(mar=c(4,4,1,1))
dE =0.1
kxmax = 8.0
c = 2
d = 3
psi2E11 = hydraulics.supplyFunction(20, psiSoil1, krmax1, kxmax, n1, alpha1, c,d, dE)
psi2E12 = hydraulics.supplyFunction(20, psiSoil2, krmax1, kxmax, n1, alpha1, c,d, dE)
psi2E13 = hydraulics.supplyFunction(20, psiSoil3, krmax1, kxmax, n1, alpha1, c,d, dE)
psi2E14 = hydraulics.supplyFunction(20, psiSoil4, krmax1, kxmax, n1, alpha1, c,d, dE)
psi2E21 = hydraulics.supplyFunction(20, psiSoil1, krmax2, kxmax, n2, alpha2, c,d, dE)
psi2E22 = hydraulics.supplyFunction(20, psiSoil2, krmax2, kxmax, n2, alpha2, c,d, dE)
psi2E23 = hydraulics.supplyFunction(20, psiSoil3, krmax2, kxmax, n2, alpha2, c,d, dE)
psi2E24 = hydraulics.supplyFunction(20, psiSoil4, krmax2, kxmax, n2, alpha2, c,d, dE)
psi2E31 = hydraulics.supplyFunction(20, psiSoil1, krmax3, kxmax, n3, alpha3, c,d, dE)
psi2E32 = hydraulics.supplyFunction(20, psiSoil2, krmax3, kxmax, n3, alpha3, c,d, dE)
psi2E33 = hydraulics.supplyFunction(20, psiSoil3, krmax3, kxmax, n3, alpha3, c,d, dE)
psi2E34 = hydraulics.supplyFunction(20, psiSoil4, krmax3, kxmax, n3, alpha3, c,d, dE)

plot(-psi2E11$PsiPlant, psi2E11$FittedE, type="l", col="black", ylab="Flow rate", xlab = "Canopy sap pressure (-MPa)",  xlim=c(0, 7), ylim=c(0,20), lwd=1)
lines(-psi2E12$PsiPlant, psi2E12$FittedE, lty=1, lwd=1)
lines(-psi2E13$PsiPlant, psi2E13$FittedE, lty=1, lwd=1)
lines(-psi2E14$PsiPlant, psi2E14$FittedE, lty=1, lwd=1)
lines(-psi2E21$PsiPlant, psi2E21$FittedE, lty=2, lwd=1)
lines(-psi2E22$PsiPlant, psi2E22$FittedE, lty=2, lwd=1)
lines(-psi2E23$PsiPlant, psi2E23$FittedE, lty=2, lwd=1)
lines(-psi2E24$PsiPlant, psi2E24$FittedE, lty=2, lwd=1)
lines(-psi2E31$PsiPlant, psi2E31$FittedE, lty=3, lwd=1)
lines(-psi2E32$PsiPlant, psi2E32$FittedE, lty=3, lwd=1)
lines(-psi2E33$PsiPlant, psi2E33$FittedE, lty=3, lwd=1)
lines(-psi2E34$PsiPlant, psi2E34$FittedE, lty=3, lwd=1)
legend("topright", bty="n", legend=textures, lwd=1, lty=1:3, cex=0.8)


###################################################
### code chunk number 9: Hydraulics.Rnw:203-217
###################################################
par(mar=c(4,4,1,1))
plot(-psi2E11$PsiPlant, psi2E11$dEdP, type="l", col="black", ylab="dEdP", xlab = "Canopy sap pressure (-MPa)",  xlim=c(0, 7), ylim=c(0,10), lwd=1)
lines(-psi2E12$PsiPlant, psi2E12$dEdP, lty=1, lwd=1)
lines(-psi2E13$PsiPlant, psi2E13$dEdP, lty=1, lwd=1)
lines(-psi2E14$PsiPlant, psi2E14$dEdP, lty=1, lwd=1)
lines(-psi2E21$PsiPlant, psi2E21$dEdP, lty=2, lwd=1)
lines(-psi2E22$PsiPlant, psi2E22$dEdP, lty=2, lwd=1)
lines(-psi2E23$PsiPlant, psi2E23$dEdP, lty=2, lwd=1)
lines(-psi2E24$PsiPlant, psi2E24$dEdP, lty=2, lwd=1)
lines(-psi2E31$PsiPlant, psi2E31$dEdP, lty=3, lwd=1)
lines(-psi2E32$PsiPlant, psi2E32$dEdP, lty=3, lwd=1)
lines(-psi2E33$PsiPlant, psi2E33$dEdP, lty=3, lwd=1)
lines(-psi2E34$PsiPlant, psi2E34$dEdP, lty=3, lwd=1)
legend("topright", bty="n", legend=textures, lwd=1, lty=1:3, cex=0.8)


###################################################
### code chunk number 10: Hydraulics.Rnw:230-279
###################################################
par(mar=c(4,4,1,1))
plot(-psiVec, unlist(lapply(psiVec, hydraulics.EXylem, -0.25, kxmax, c, d)), type="l", ylab="Transpiration rate", xlab="Canopy sap pressure (-MPa)", xlim = c(0,5))
for(psi in seq(-0.5, -4.5, by=-0.25)) {
  lines(-psiVec, unlist(lapply(psiVec, hydraulics.EXylem, psi, kxmax, c, d)))
}
Ein1 = 4
Ein2 = 8
Ein3 = 15
psiVec = seq(-0.25, -4.5, by =-0.25)
LU1 = rep(NA, length(psiVec))
LR1 = rep(NA, length(psiVec))
PU1 = rep(NA, length(psiVec))
PR1 = rep(NA, length(psiVec))
G1 = rep(NA, length(psiVec))
LU2 = rep(NA, length(psiVec))
LR2 = rep(NA, length(psiVec))
PU2 = rep(NA, length(psiVec))
PR2 = rep(NA, length(psiVec))
G2 = rep(NA, length(psiVec))
LU3 = rep(NA, length(psiVec))
LR3 = rep(NA, length(psiVec))
PU3 = rep(NA, length(psiVec))
PR3 = rep(NA, length(psiVec))
G3 = rep(NA, length(psiVec))
for(i in 1:length(psiVec)) {
  rp1 = hydraulics.regulatedPsiXylem(Ein1, psiVec[i], kxmax, c,d)
  PU1[i] = rp1[1]
  PR1[i] = rp1[2]
  LU1[i] = rp1[3]
  LR1[i] = rp1[4]
  G1[i] = rp1[6]
  rp2 = hydraulics.regulatedPsiXylem(Ein2, psiVec[i], kxmax, c,d)
  PU2[i] = rp2[1]
  PR2[i] = rp2[2]
  LU2[i] = rp2[3]
  LR2[i] = rp2[4]
  G2[i] = rp2[6]
  rp3 = hydraulics.regulatedPsiXylem(Ein3, psiVec[i], kxmax, c,d)
  PU3[i] = rp3[1]
  PR3[i] = rp3[2]
  LU3[i] = rp3[3]
  LR3[i] = rp3[4]
  G3[i] = rp3[6]
}
lines(-PR1, LR1, col="gray", lwd=1.5, lty=1)
lines(-PR2, LR2, col="gray", lwd=1.5, lty=2)
lines(-PR3, LR3, col="gray", lwd=1.5, lty=3)
legend("topleft", col=rep("gray",3), lty=1:3, lwd=1.5, 
       legend=c(paste("E\'=",Ein1,sep=""), paste("E\'=",Ein2,sep=""), paste("E\'=",Ein3,sep="")), bty="n", cex=0.8)


###################################################
### code chunk number 11: Hydraulics.Rnw:284-337
###################################################
 par(mar=c(4,4,3,1))
 Ein1 = 4
 Ein2 = 8
 Ein3 = 15
 c = 2
 d = 3
 dE = 0.1
 psiVec = seq(-0.01, -2.0, by =-0.03)
 LU1 = rep(NA, length(psiVec))
 LR1 = rep(NA, length(psiVec))
 PU1 = rep(NA, length(psiVec))
 PR1 = rep(NA, length(psiVec))
 G1 = rep(NA, length(psiVec))
 LU2 = rep(NA, length(psiVec))
 LR2 = rep(NA, length(psiVec))
 PU2 = rep(NA, length(psiVec))
 PR2 = rep(NA, length(psiVec))
 G2 = rep(NA, length(psiVec))
 LU3 = rep(NA, length(psiVec))
 LR3 = rep(NA, length(psiVec))
 PU3 = rep(NA, length(psiVec))
 PR3 = rep(NA, length(psiVec))
 G3 = rep(NA, length(psiVec))
 for(i in 1:length(psiVec)) {
   rp1 = hydraulics.regulatedPsiTwoElements(Ein1, psiVec[i], krmax1, kxmax, n1, alpha1, c, d,dE)
   PU1[i] = rp1[1]
   PR1[i] = rp1[2]
   LU1[i] = rp1[3]
   LR1[i] = rp1[4]
   G1[i] = LR1[i]/LU1[i]
   rp2 = hydraulics.regulatedPsiTwoElements(Ein2, psiVec[i], krmax1, kxmax, n1, alpha1, c,d,dE)
   PU2[i] = rp2[1]
   PR2[i] = rp2[2]
   LU2[i] = rp2[3]
   LR2[i] = rp2[4]
   G2[i] = LR2[i]/LU2[i]
   rp3 = hydraulics.regulatedPsiTwoElements(Ein3, psiVec[i], krmax1, kxmax, n1, alpha1, c,d,dE)
   PU3[i] = rp3[1]
   PR3[i] = rp3[2]
   LU3[i] = rp3[3]
   LR3[i] = rp3[4]
   G3[i] = LR3[i]/LU3[i]
 }
 plot(-PR1, LR1, type="l", col="black", ylab="Transpiration rate", xlab = "Canopy sap pressure (-MPa)",  main = textures[1],  xlim=c(0, 2.0), ylim=c(0,10), lwd=1.5, lty=1)
 lines(-PR2, LR2, col="gray", lwd=1.5, lty=1)
 lines(-PR3, LR3, col="red", lwd=1.5, lty=1)
 for(i in 1:length(psiVec)) {
   pp = hydraulics.supplyFunction(15, psiVec[i], krmax1, kxmax, n1, alpha1, c,d, dE = dE)
   lines(-pp$PsiPlant, pp$FittedE, col="black", lwd=1, lty=1)
 }
  legend("topleft", col=c("black","gray","red"), lty=1, lwd=1.5, 
        legend=c(paste("E\'=",Ein1,sep=""), paste("E\'=",Ein2,sep=""), paste("E\'=",Ein3,sep="")), 
        bty="n", cex=0.8)


###################################################
### code chunk number 12: Hydraulics.Rnw:341-349
###################################################
  par(mar=c(4,4,3,1))
 plot(-psiVec, -(PR1-psiVec), type="l", asp= 1, col="black", ylab="Regulated pressure drop (-MPa)", xlab = "Soil pressure (-MPa)",  main = textures[1],  xlim=c(0, 2.0), ylim=c(0,2), lwd=1.5, lty=1)
 lines(-psiVec, -(PR2-psiVec), col="gray", lwd=1.5, lty=1)
 lines(-psiVec, -(PR3-psiVec), col="red", lwd=1.5, lty=1)
 for(a in seq(0, 5, by=0.25)) abline(a=a, b=-1, col="gray")
  legend("topleft", col=c("black","gray","red"), lty=1, lwd=1.5, 
        legend=c(paste("E\'=",Ein1,sep=""), paste("E\'=",Ein2,sep=""), paste("E\'=",Ein3,sep="")), 
        bty="n", cex=0.8)


###################################################
### code chunk number 13: Hydraulics.Rnw:364-418
###################################################
 par(mar=c(4,4,3,1))
 Ein1 = 4
 Ein2 = 8
 Ein3 = 15
 c = 2
 d = 3
 dE = 0.1
 psiVec = seq(-0.25, -4.5, by =-0.25)
 LU1 = rep(NA, length(psiVec))
 LR1 = rep(NA, length(psiVec))
 PU1 = rep(NA, length(psiVec))
 PR1 = rep(NA, length(psiVec))
 G1 = rep(NA, length(psiVec))
 LU2 = rep(NA, length(psiVec))
 LR2 = rep(NA, length(psiVec))
 PU2 = rep(NA, length(psiVec))
 PR2 = rep(NA, length(psiVec))
 G2 = rep(NA, length(psiVec))
 LU3 = rep(NA, length(psiVec))
 LR3 = rep(NA, length(psiVec))
 PU3 = rep(NA, length(psiVec))
 PR3 = rep(NA, length(psiVec))
 G3 = rep(NA, length(psiVec))
 for(i in 1:length(psiVec)) {
   rp1 = hydraulics.regulatedPsiTwoElements(Ein1, psiVec[i], krmax2, kxmax, n2, alpha2, c, d, dE)
   PU1[i] = rp1[1]
   PR1[i] = rp1[2]
   LU1[i] = rp1[3]
   LR1[i] = rp1[4]
   G1[i] = LR1[i]/LU1[i]
   rp2 = hydraulics.regulatedPsiTwoElements(Ein2, psiVec[i], krmax2, kxmax, n2, alpha2, c,d, dE)
   PU2[i] = rp2[1]
   PR2[i] = rp2[2]
   LU2[i] = rp2[3]
   LR2[i] = rp2[4]
   G2[i] = LR2[i]/LU2[i]
   rp3 = hydraulics.regulatedPsiTwoElements(Ein3, psiVec[i], krmax2, kxmax, n2, alpha2, c,d, dE)
   PU3[i] = rp3[1]
   PR3[i] = rp3[2]
   LU3[i] = rp3[3]
   LR3[i] = rp3[4]
   G3[i] = LR3[i]/LU3[i]
 }
 plot(-PR1, LR1, type="l", col="black", ylab="Transpiration rate", xlab = "Canopy sap pressure (-MPa)",   main = textures[2], xlim=c(0, 5), ylim=c(0,15), lwd=1.5, lty=2)
 lines(-PR2, LR2, col="gray", lwd=1.5, lty=2)
 lines(-PR3, LR3, col="red", lwd=1.5, lty=2)
 for(i in 1:length(psiVec)) {
   pp = hydraulics.supplyFunction(15, psiVec[i], krmax2, kxmax, n2, alpha2, c,d, dE = dE)
   lines(-pp$PsiPlant, pp$FittedE, col="black", lwd=1, lty=1)
 }
 
 legend("topleft", col=c("black","gray","red"), lty=2, lwd=1.5, 
        legend=c(paste("E\'=",Ein1,sep=""), paste("E\'=",Ein2,sep=""), paste("E\'=",Ein3,sep="")), 
        bty="n", cex=0.8)


###################################################
### code chunk number 14: Hydraulics.Rnw:422-430
###################################################
  par(mar=c(4,4,3,1))
 plot(-psiVec, -(PR1-psiVec), type="l", asp= 1, col="black", ylab="Regulated pressure drop (-MPa)", xlab = "Soil pressure (-MPa)",  main = textures[2],  xlim=c(0, 2.0), ylim=c(0,2), lwd=1.5, lty=2)
 lines(-psiVec, -(PR2-psiVec), col="gray", lwd=1.5, lty=2)
 lines(-psiVec, -(PR3-psiVec), col="red", lwd=1.5, lty=2)
 for(a in seq(0, 5, by=0.25)) abline(a=a, b=-1, col="gray")
  legend("topleft", col=c("black","gray","red"), lty=2, lwd=1.5, 
        legend=c(paste("E\'=",Ein1,sep=""), paste("E\'=",Ein2,sep=""), paste("E\'=",Ein3,sep="")), 
        bty="n", cex=0.8)


###################################################
### code chunk number 15: Hydraulics.Rnw:446-493
###################################################
 par(mar=c(4,4,3,1))
 psiVec = seq(-0.25, -4.5, by =-0.25)
LU1 = rep(NA, length(psiVec))
 LR1 = rep(NA, length(psiVec))
 PU1 = rep(NA, length(psiVec))
 PR1 = rep(NA, length(psiVec))
 G1 = rep(NA, length(psiVec))
 LU2 = rep(NA, length(psiVec))
 LR2 = rep(NA, length(psiVec))
 PU2 = rep(NA, length(psiVec))
 PR2 = rep(NA, length(psiVec))
 G2 = rep(NA, length(psiVec))
 LU3 = rep(NA, length(psiVec))
 LR3 = rep(NA, length(psiVec))
 PU3 = rep(NA, length(psiVec))
 PR3 = rep(NA, length(psiVec))
 G3 = rep(NA, length(psiVec))
 for(i in 1:length(psiVec)) {
   rp1 = hydraulics.regulatedPsiTwoElements(Ein1, psiVec[i], krmax3, kxmax, n3, alpha3, c, d)
   PU1[i] = rp1[1]
   PR1[i] = rp1[2]
   LU1[i] = rp1[3]
   LR1[i] = rp1[4]
   G1[i] = LR1[i]/LU1[i]
   rp2 = hydraulics.regulatedPsiTwoElements(Ein2, psiVec[i], krmax3, kxmax, n3, alpha3, c,d)
   PU2[i] = rp2[1]
   PR2[i] = rp2[2]
   LU2[i] = rp2[3]
   LR2[i] = rp2[4]
   G2[i] = LR2[i]/LU2[i]
   rp3 = hydraulics.regulatedPsiTwoElements(Ein3, psiVec[i], krmax3, kxmax, n3, alpha3, c,d)
   PU3[i] = rp3[1]
   PR3[i] = rp3[2]
   LU3[i] = rp3[3]
   LR3[i] = rp3[4]
   G3[i] = LR3[i]/LU3[i]
 }
 plot(-PR1, LR1, type="l", col="black", ylab="Transpiration rate", xlab = "Canopy sap pressure (-MPa)",  main = textures[3], xlim=c(0, 5), ylim=c(0,15), lwd=1.5, lty=3)
 lines(-PR2, LR2, col="gray", lwd=1.5, lty=3)
 lines(-PR3, LR3, col="red", lwd=1.5, lty=3)
 for(i in 1:length(psiVec)) {
   pp = hydraulics.supplyFunction(15, psiVec[i], krmax3, kxmax, n3, alpha3, c,d, dE = 0.1)
   lines(-pp$PsiPlant, pp$E, col="black", lwd=1, lty=1)
 }
 legend("topleft", col=c("black","gray","red"), lty=3, lwd=1.5, 
        legend=c(paste("E=",Ein1,sep=""), paste("E=",Ein2,sep=""), paste("E=",Ein3,sep="")), 
        bty="n", cex=0.8)


###################################################
### code chunk number 16: Hydraulics.Rnw:497-505
###################################################
  par(mar=c(4,4,3,1))
 plot(-psiVec, -(PR1-psiVec), type="l", asp= 1, col="black", ylab="Regulated pressure drop (-MPa)", xlab = "Soil pressure (-MPa)",  main = textures[3],  xlim=c(0, 2.0), ylim=c(0,2), lwd=1.5, lty=3)
 lines(-psiVec, -(PR2-psiVec), col="gray", lwd=1.5, lty=3)
 lines(-psiVec, -(PR3-psiVec), col="red", lwd=1.5, lty=3)
 for(a in seq(0, 5, by=0.25)) abline(a=a, b=-1, col="gray")
  legend("topleft", col=c("black","gray","red"), lty=3, lwd=1.5, 
        legend=c(paste("E\'=",Ein1,sep=""), paste("E\'=",Ein2,sep=""), paste("E\'=",Ein3,sep="")), 
        bty="n", cex=0.8)


###################################################
### code chunk number 17: Hydraulics.Rnw:523-539
###################################################
par(mar=c(4,4,1,1))
psiVec = seq(-0.1, -7.0, by =-0.01)
k = unlist(lapply(psiVec, hydraulics.xylemConductance, kxmax, c, d))
plot(-psiVec, k, type="l",ylab="Hydraulic conductance", xlab="Pressure (-MPa)", lwd=1.5)

krmax21 = 16000000000
krmax22 = 8000000000
krmax23 = 4000000000
k1 = unlist(lapply(psiVec, hydraulics.vanGenuchtenConductance, krmax21, n2, alpha2))
k2 = unlist(lapply(psiVec, hydraulics.vanGenuchtenConductance, krmax22, n2, alpha2))
k3 = unlist(lapply(psiVec, hydraulics.vanGenuchtenConductance, krmax23, n2, alpha2))
lines(-psiVec, k1, lty=3, lwd=1.5, col="black")
lines(-psiVec, k2, lty=3, lwd=1.5, col="gray")
lines(-psiVec, k3, lty=3, lwd=1.5, col="red")
 legend("topright", col=c("black","black","gray","red"), lty=c(1,3,3,3), lwd=1.5, legend=c("xylem curve", paste("krmax=",krmax21,sep=""), paste("krmax=",krmax22,sep=""), paste("krmax=",krmax23,sep="")), 
        bty="n", cex=0.8)


###################################################
### code chunk number 18: Hydraulics.Rnw:544-599
###################################################
par(mar=c(4,4,1,1))
Ein2 = 4
psiVec = seq(-0.01, -5.0, by =-0.01)
LU1 = rep(NA, length(psiVec))
LR1 = rep(NA, length(psiVec))
PU1 = rep(NA, length(psiVec))
PR1 = rep(NA, length(psiVec))
G1 = rep(NA, length(psiVec))
 LU2 = rep(NA, length(psiVec))
 LR2 = rep(NA, length(psiVec))
 PU2 = rep(NA, length(psiVec))
 PR2 = rep(NA, length(psiVec))
 G2 = rep(NA, length(psiVec))
 LU3 = rep(NA, length(psiVec))
 LR3 = rep(NA, length(psiVec))
 PU3 = rep(NA, length(psiVec))
 PR3 = rep(NA, length(psiVec))
 G3 = rep(NA, length(psiVec))
 for(i in 1:length(psiVec)) {
   rp1 =  hydraulics.regulatedPsiXylem(Ein2, psiVec[i], kxmax, c,d)
   PU1[i] = rp1[1]
   PR1[i] = rp1[2]
   LU1[i] = rp1[3]
   LR1[i] = rp1[4]
   G1[i] = LR1[i]/LU1[i]
}
plot(-PR1, LR1, type="l", col="black", ylab="Transpiration rate", xlab = "Canopy sap pressure (-MPa)", xlim=c(0, 5), ylim=c(0,4), lwd=1.5, lty=1)

 for(i in 1:length(psiVec)) {
   rp1 = hydraulics.regulatedPsiTwoElements(Ein2, psiVec[i], krmax21, kxmax, n2, alpha2, c, d, dE = 0.1)
   PU1[i] = rp1[1]
   PR1[i] = rp1[2]
   LU1[i] = rp1[3]
   LR1[i] = rp1[4]
   G1[i] = LR1[i]/LU1[i]
   rp2 = hydraulics.regulatedPsiTwoElements(Ein2, psiVec[i], krmax22, kxmax, n2, alpha2, c,d,dE = 0.1)
   PU2[i] = rp2[1]
   PR2[i] = rp2[2]
   LU2[i] = rp2[3]
   LR2[i] = rp2[4]
   G2[i] = LR2[i]/LU2[i]
   rp3 = hydraulics.regulatedPsiTwoElements(Ein2, psiVec[i], krmax23, kxmax, n2, alpha2, c,d,dE = 0.1)
   PU3[i] = rp3[1]
   PR3[i] = rp3[2]
   LU3[i] = rp3[3]
   LR3[i] = rp3[4]
   G3[i] = LR3[i]/LU3[i]
 }
 lines(-PR1, LR1, col="black", lwd=1.5, lty=3)
 lines(-PR2, LR2, col="gray", lwd=1.5, lty=3)
 lines(-PR3, LR3, col="red", lwd=1.5, lty=3)
 legend("topright", col=c("black", "black","gray","red"), lty=c(1,3,3,3), lwd=1.5, 
        legend=c("xylem only", paste("krmax=",krmax1,sep=""), paste("krmax=",krmax2,sep=""), paste("krmax=",krmax3,sep="")), 
        bty="n", cex=0.8)
 


###################################################
### code chunk number 19: Hydraulics.Rnw:608-656
###################################################
data(examplemeteo)
par(mar=c(4,4,1,1))
psiSoil = -1.5
E = examplemeteo$PET
LU1 = rep(NA, length(E))
LR1 = rep(NA, length(E))
PU1 = rep(NA, length(E))
PR1 = rep(NA, length(E))
G1 = rep(NA, length(E))
 LU2 = rep(NA, length(E))
 LR2 = rep(NA, length(E))
 PU2 = rep(NA, length(E))
 PR2 = rep(NA, length(E))
 G2 = rep(NA, length(E))
 LU3 = rep(NA, length(E))
 LR3 = rep(NA, length(E))
 PU3 = rep(NA, length(E))
 PR3 = rep(NA, length(E))
 G3 = rep(NA, length(E))
 for(i in 1:length(E)) {
   rp1 = hydraulics.regulatedPsiTwoElements(E[i], psiSoil, krmax1, kxmax, n1, alpha1, c, d,dE)
   PU1[i] = rp1[1]
   PR1[i] = rp1[2]
   LU1[i] = rp1[3]
   LR1[i] = rp1[4]
   G1[i] = LR1[i]/E[i]
   if(E[i]==0) G1[i] = NA
   rp2 = hydraulics.regulatedPsiTwoElements(E[i], psiSoil, krmax2, kxmax, n2, alpha2, c,d,dE)
   PU2[i] = rp2[1]
   PR2[i] = rp2[2]
   LU2[i] = rp2[3]
   LR2[i] = rp2[4]
   G2[i] = LR2[i]/E[i]
   if(E[i]==0) G2[i] = NA
   rp3 = hydraulics.regulatedPsiTwoElements(E[i], psiSoil, krmax3, kxmax, n3, alpha3, c,d,dE)
   PU3[i] = rp3[1]
   PR3[i] = rp3[2]
   LU3[i] = rp3[3]
   LR3[i] = rp3[4]
   G3[i] = LR3[i]/E[i]
   if(E[i]==0) G3[i] = NA
 }

plot(LR1, type="l", col="black", ylab="Transpiration rate", xlab = "DOY", ylim=c(0,8), lwd=1, lty=1)
lines(examplemeteo$PET, lty=1, col="gray", lwd=1)
lines(LR2, lty=2, col="black", lwd=1)
lines(LR3, lty=3, col="black", lwd=1)
legend("topleft", col=c("gray","black","black","black"), legend=c("PET", textures), lty=c(1,1,2,3), bty="n", cex=0.8)


###################################################
### code chunk number 20: Hydraulics.Rnw:661-666
###################################################
par(mar=c(4,4,1,1))
plot(PU1, type="l", col="black", ylab="Unregulated Psi_canopy", ylim=c(-10,0), xlab = "DOY", lwd=1, lty=1)
lines(PU2, col="black", lty=2)
lines(PU3, col="black", lty=3)
legend("bottomleft", col="black", legend=textures, lty=c(1,2,3), bty="n", cex=0.8)


###################################################
### code chunk number 21: Hydraulics.Rnw:670-675
###################################################
par(mar=c(4,4,1,1))
plot(PR1, type="l", col="black", ylab="Regulated Psi_canopy", ylim=c(-2.5,-1.5), xlab = "DOY", lwd=1, lty=1)
lines(PR2, col="black", lty=2)
lines(PR3, col="black", lty=3)
legend("bottomleft", col="black", legend=textures, lty=c(1,2,3), bty="n", cex=0.8)


###################################################
### code chunk number 22: Hydraulics.Rnw:680-685
###################################################
par(mar=c(4,4,1,1))
plot(G1, type="l", col="black", ylab="Relative G", ylim=c(0,1), xlab = "DOY", lwd=1, lty=1)
lines(G2, col="black", lty=2)
lines(G3, col="black", lty=3)
legend("bottomleft", col="black", legend=textures, lty=c(1,2,3), bty="n", cex=0.8)


###################################################
### code chunk number 23: Hydraulics.Rnw:691-695
###################################################
par(mar=c(4,4,1,1))
E = examplemeteo$PET
P = examplemeteo$Precipitation
barplot(P, col="black", ylab="Precipitation", xlab = "DOY")


###################################################
### code chunk number 24: Hydraulics.Rnw:700-774
###################################################
par(mar=c(4,4,1,1))
sand1 = 65
sand2 = 30
sand3 = 10
clay1 = 30
clay2 = 10
clay3 = 50
FC1 = soil.psi2theta(clay1, sand1, -33)
FC2 = soil.psi2theta(clay2, sand2, -33)
FC3 = soil.psi2theta(clay3, sand3, -33)
V1 = FC1*1000
V2 = FC2*1000
V3 = FC3*1000
psiSoil1 = rep(NA, length(E))
psiSoil2 = rep(NA, length(E))
psiSoil3 = rep(NA, length(E))
W1 = rep(1, length(E))
W2 = rep(1, length(E))
W3 = rep(1, length(E))
LU1 = rep(NA, length(E))
LR1 = rep(NA, length(E))
PU1 = rep(NA, length(E))
PR1 = rep(NA, length(E))
G1 = rep(NA, length(E))
 LU2 = rep(NA, length(E))
 LR2 = rep(NA, length(E))
 PU2 = rep(NA, length(E))
 PR2 = rep(NA, length(E))
 G2 = rep(NA, length(E))
 LU3 = rep(NA, length(E))
 LR3 = rep(NA, length(E))
 PU3 = rep(NA, length(E))
 PR3 = rep(NA, length(E))
 G3 = rep(NA, length(E))
 for(i in 1:length(E)) {
   W1[i] = min((W1[i]*V1+P[i])/V1,1.0)
   W2[i] = min((W2[i]*V2+P[i])/V2,1.0)
   W3[i] = min((W3[i]*V3+P[i])/V3,1.0)
   psiSoil1[i] = soil.theta2psi(clay1,sand1,W1[i]*FC1)/1000
   psiSoil2[i] = soil.theta2psi(clay2,sand2,W2[i]*FC2)/1000
   psiSoil3[i] = soil.theta2psi(clay3,sand3,W3[i]*FC3)/1000
   rp1 = hydraulics.regulatedPsiTwoElements(E[i], psiSoil1[i], krmax1, kxmax, n1, alpha1, c, d,dE)
   PU1[i] = rp1[1]
   PR1[i] = rp1[2]
   LU1[i] = rp1[3]
   LR1[i] = rp1[4]
   G1[i] = LR1[i]/E[i]
   if(E[i]==0) G1[i] = NA
   rp2 = hydraulics.regulatedPsiTwoElements(E[i], psiSoil2[i], krmax2, kxmax, n2, alpha2, c,d,dE)
   PU2[i] = rp2[1]
   PR2[i] = rp2[2]
   LU2[i] = rp2[3]
   LR2[i] = rp2[4]
   G2[i] = LR2[i]/E[i]
   if(E[i]==0) G2[i] = NA
   rp3 = hydraulics.regulatedPsiTwoElements(E[i], psiSoil3[i], krmax3, kxmax, n3, alpha3, c,d,dE)
   PU3[i] = rp3[1]
   PR3[i] = rp3[2]
   LU3[i] = rp3[3]
   LR3[i] = rp3[4]
   G3[i] = LR3[i]/E[i]
   if(E[i]==0) G3[i] = NA
   if(i<length(E)) {
    W1[i+1] = min((W1[i]*V1-LR1[i])/V1,1.0)
    W2[i+1] = min((W2[i]*V2-LR2[i])/V2,1.0)
    W3[i+1] = min((W3[i]*V3-LR3[i])/V3,1.0)
   }
 }

plot(LR1, type="l", col="black", ylab="Transpiration rate", xlab = "DOY", ylim=c(0,8), lwd=1, lty=1)
lines(examplemeteo$PET, lty=1, col="gray", lwd=1)
lines(LR2, lty=2, col="black", lwd=1)
lines(LR3, lty=3, col="black", lwd=1)
legend("topleft", col=c("gray","black","black","black"), legend=c("PET", textures), lty=c(1,1,2,3), bty="n", cex=0.8)


###################################################
### code chunk number 25: Hydraulics.Rnw:779-784
###################################################
par(mar=c(4,4,1,1))
plot(W1, type="l", col="black", ylab="W", ylim=c(0,1), xlab = "DOY", lwd=1, lty=1)
lines(W2, col="black", lty=2)
lines(W3, col="black", lty=3)
legend("bottomleft", col="black", legend=textures, lty=c(1,2,3), bty="n", cex=0.8)


###################################################
### code chunk number 26: Hydraulics.Rnw:790-795
###################################################
par(mar=c(4,4,1,1))
plot(psiSoil1, type="l", col="black", ylab="Psi_soil", ylim=c(-10,0), xlab = "DOY", lwd=1, lty=1)
lines(psiSoil2, col="black", lty=2)
lines(psiSoil3, col="black", lty=3)
legend("bottomleft", col="black", legend=textures, lty=c(1,2,3), bty="n", cex=0.8)


###################################################
### code chunk number 27: Hydraulics.Rnw:799-804
###################################################
par(mar=c(4,4,1,1))
plot(PR1, type="l", col="black", ylab="Regulated Psi_canopy", ylim=c(-4,-0), xlab = "DOY", lwd=1, lty=1)
lines(PR2, col="black", lty=2)
lines(PR3, col="black", lty=3)
legend("bottomleft", col="black", legend=textures, lty=c(1,2,3), bty="n", cex=0.8)


###################################################
### code chunk number 28: Hydraulics.Rnw:809-814
###################################################
par(mar=c(4,4,1,1))
plot(G1, type="l", col="black", ylab="Relative G", ylim=c(0,1), xlab = "DOY", lwd=1, lty=1)
lines(G2, col="black", lty=2)
lines(G3, col="black", lty=3)
legend("bottomleft", col="black", legend=textures, lty=c(1,2,3), bty="n", cex=0.8)


