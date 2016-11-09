### R code from vignette source 'SoilWaterBalance.Rnw'
### Encoding: ISO8859-1

###################################################
### code chunk number 1: SoilWaterBalance.Rnw:21-23
###################################################
options(width=67)
library(medfate)


###################################################
### code chunk number 2: SoilWaterBalance.Rnw:189-230
###################################################
par(mfrow=c(2,2), mar=c(5,5,5,1))
throughfallMatrixGash<-function(P = seq(1,50, by=1), Cm = seq(1,5, by=1), 
                                ER = 0.08,p=0.8) {
  m2<-P-interception(P,Cm[1],p,ER=ER)
  for(i in 2:length(Cm)) {
    m2<-rbind(m2,P-interception(P,Cm[i],p,ER=ER))
  }
  colnames(m2)<-P
  rownames(m2)<-Cm
  return(m2)
}

Cm = c(0.5,seq(1,4, by=1))
P = seq(1,50, by=1)

m2 = throughfallMatrixGash(P=P, p=0.2, Cm=Cm,ER = 0.05)
rt = sweep(m2,2,P,"/")*100
matplot(t(rt), type="l", axes=TRUE, ylab="Relative throughfall (%)", 
        xlab="Gross rainfall (mm)", xlim=c(0,length(P)), 
        lty=1:length(Cm), col="black", ylim=c(0,100))
title(main="p = 0.2 E/R = 0.05")
m2 = throughfallMatrixGash(P=P, p=0.8, Cm=Cm,ER = 0.05)
rt = sweep(m2,2,P,"/")*100
matplot(t(rt), type="l", axes=TRUE, ylab="Relative throughfall (%)", 
        xlab="Gross rainfall (mm)", xlim=c(0,length(P)), 
        lty=1:length(Cm), col="black", ylim=c(0,100))
title(main="p = 0.8 E/R = 0.05")
m2 = throughfallMatrixGash(P=P, p=0.2, Cm=Cm,ER = 0.2)
rt = sweep(m2,2,P,"/")*100
matplot(t(rt), type="l", axes=TRUE, ylab="Relative throughfall (%)", 
        xlab="Gross rainfall (mm)", xlim=c(0,length(P)), 
        lty=1:length(Cm), col="black", ylim=c(0,100))
title(main="p = 0.2 E/R = 0.2")
m2 = throughfallMatrixGash(P=P, p=0.8, Cm=Cm,ER = 0.2)
rt = sweep(m2,2,P,"/")*100
matplot(t(rt), type="l", axes=TRUE, ylab="Relative throughfall (%)", 
        xlab="Gross rainfall (mm)", xlim=c(0,length(P)), 
        lty=1:length(Cm), col="black", ylim=c(0,100))
title(main="p = 0.8 E/R = 0.2")

legend("bottomright",lty=1:length(Cm), legend=paste("S =",Cm), bty="n")


###################################################
### code chunk number 3: SoilWaterBalance.Rnw:239-287
###################################################
par(mfrow=c(1,2), mar=c(5,5,5,1))

SoilDepth = c(200,400,800,1200,1500)

#TOPSOIL LAYERS
d1 = pmin(SoilDepth, 300) #<300
#SUBSOIL LAYERS
d2 = pmax(0, pmin(SoilDepth-300,1200)) #300-1500 mm
#ROCK LAYER
d3 = 4000-(d1+d2) #From SoilDepth down to 4.0 m

TS_clay = 15
TS_sand = 25
SS_clay = 15
SS_sand = 25
RL_clay = 15
RL_sand = 25
TS_gravel = 20
SS_gravel = 40
RL_gravel = 95

Theta_FC1=soil.psi2theta(TS_clay, TS_sand, -33) #in m3/m3
Theta_FC2=soil.psi2theta(SS_clay, SS_sand, -33) #in m3/m3
Theta_FC3=soil.psi2theta(RL_clay, RL_sand, -33) #in m3/m3
pcTS_gravel = 1-(TS_gravel/100)
pcSS_gravel = 1-(SS_gravel/100)
pcRL_gravel = 1-(RL_gravel/100)
MaxVol1 = (d1*Theta_FC1*pcTS_gravel)
MaxVol2 = (d2*Theta_FC2*pcSS_gravel)
MaxVol3 = (d3*Theta_FC3*pcRL_gravel)
V = MaxVol1+MaxVol2+MaxVol3

par(mar=c(5,5,1,1), mfrow=c(1,2))
NP = seq(0,60, by=1)
plot(NP,soil.infiltration(NP, V[1]), type="l", xlim=c(0,60), ylim=c(0,60), 
     ylab="Infiltration (mm)", xlab="Net rainfall (mm)", frame=FALSE)
lines(NP,soil.infiltration(NP, V[2]), lty=2)
lines(NP,soil.infiltration(NP, V[3]), lty=3)
lines(NP,soil.infiltration(NP, V[4]), lty=4)
lines(NP,soil.infiltration(NP, V[5]), lty=5)
legend("topleft", bty="n", lty=1:5, legend=c(paste("d =", SoilDepth, "Vsoil =",round(V),"mm")))
plot(NP,NP-soil.infiltration(NP, V[1]), type="l", xlim=c(0,60), ylim=c(0,60), 
     ylab="Runoff (mm)", xlab="Net rainfall (mm)", frame=FALSE)
lines(NP,NP-soil.infiltration(NP, V[2]), lty=2)
lines(NP,NP-soil.infiltration(NP, V[3]), lty=3)
lines(NP,NP-soil.infiltration(NP, V[4]), lty=4)
lines(NP,NP-soil.infiltration(NP, V[5]), lty=5)
legend("topleft", bty="n", lty=1:5, legend=c(paste("d =", SoilDepth,"Vsoil =",round(V),"mm")))


###################################################
### code chunk number 4: SoilWaterBalance.Rnw:310-314
###################################################
par(mar=c(4,4,1,1))
LAIc = seq(0,10, by=0.01)
TmaxPET = -0.006*(LAIc^2) + 0.134*LAIc + 0.036
plot(LAIc, TmaxPET, type="l", ylab="PTcanopy/PET", xlab="LAIstand", ylim=c(0,1))


###################################################
### code chunk number 5: SoilWaterBalance.Rnw:375-385
###################################################
par(mar=c(4,4,1,1))
x = seq(-10000, 0, by=10)
plot(x/1000,exp(-0.6931472*(x/-2000)^3)*100, type="l", ylab="K (relative conductance, %)",
     xlim=c(-10.0,0), ylim=c(0,100),xlab="Y (soil potential, MPa)", frame=FALSE)
lines(x/1000, exp(-0.6931472*(x/-3000)^3)*100, lty=2)
lines(x/1000, exp(-0.6931472*(x/-4000)^3)*100, lty=3)
lines(x/1000, exp(-0.6931472*(x/-5000)^3)*100, lty=4)
legend("topleft", lty=1:4, col=c(rep("black",4)), 
       legend = paste("Ysp = ", c(-2.0,-3.0,-4.0, -5.0), "MPa"), bty="n")
abline(h=50, col="gray", lwd=2)


###################################################
### code chunk number 6: SoilWaterBalance.Rnw:468-559
###################################################
TS_clay=10
TS_silt=65
TS_sand=25
TS_gravel=40
SS_clay=10
SS_silt=65
SS_sand = 25
SS_gravel=40
TS_macro=0.25
TS_micro = 0.75
SS_macro=0.10
SS_micro=0.90
#Rock layer is like subsoil but with 95% of rocks
RL_clay = SS_clay
RL_sand = SS_sand
RL_macro = SS_macro
RL_micro = SS_micro
RL_gravel = 95


RunEvaporation<-function(Gsoil, Ksoil, d1,d2,d3, numDays = 15){
  PET = 100 #Not limited by PET
  Lground = 1
  
  Theta_FC1=soil.psi2theta(TS_clay, TS_sand, -33) #in m3/m3
  Theta_FC2=soil.psi2theta(SS_clay, SS_sand, -33) #in m3/m3
  Theta_FC3=soil.psi2theta(RL_clay, RL_sand, -33) #in m3/m3
  pcTS_gravel = 1-(TS_gravel/100)
  pcSS_gravel = 1-(SS_gravel/100)
  pcRL_gravel = 1-(RL_gravel/100)
  MaxVol1 = (d1*Theta_FC1*pcTS_gravel)
  MaxVol2 = (d2*Theta_FC2*pcSS_gravel)
  MaxVol3 = (d3*Theta_FC3*pcRL_gravel)
  Ssoil = MaxVol1 + MaxVol2 + MaxVol3

  W1=rep(0, numDays)
  W2=rep(0, numDays)
  W3=rep(0, numDays)
  W1[1] = 1
  W2[1] = 1
  W3[1] = 1
  Esoil = rep(NA,numDays)
  EsoilCum = rep(NA,numDays)
  t = rep(NA, numDays)
  for(i in 1:numDays){
    #Evaporation from bare soil
    Esoil[i] = soil.evaporation(DEF=(MaxVol1*(1 - W1[i])), PETs = PET*Lground, Gsoil = Gsoil)
    if(i==1) EsoilCum[i] = Esoil[i]
    else EsoilCum[i] = EsoilCum[i-1]+Esoil[i]
    #Exponential decay to divide bare soil evaporation among layers
    Esoil1 = Esoil[i]*(1-exp(-Ksoil*d1))
    Esoil2 = Esoil[i]*(exp(-Ksoil*d1)-exp(-Ksoil*(d1+d2)))
    Esoil3 = Esoil[i]*(exp(-Ksoil*(d1+d2)))
    if(i<numDays){
      W1[i+1] = max(W1[i]-(Esoil1)/MaxVol1,0)
      W2[i+1] = max(min(W2[i]-(Esoil2)/MaxVol2,1),0)
      W3[i+1] = max(min(W3[i]-(Esoil3)/MaxVol3,1),0)
    }  
  }
  return(list(Esoil = Esoil, EsoilCum = EsoilCum))  
}

E11=RunEvaporation(Gsoil=1, Ksoil = 0.05, d1=300, d2=1200, d3= 2500)
E12=RunEvaporation(Gsoil=2, Ksoil = 0.05, d1=300, d2=1200, d3= 2500)
E13=RunEvaporation(Gsoil=3, Ksoil = 0.05, d1=300, d2=1200, d3= 2500)
E21=RunEvaporation(Gsoil=1, Ksoil = 0.005, d1=300, d2=1200, d3= 2500)
E22=RunEvaporation(Gsoil=2, Ksoil = 0.005, d1=300, d2=1200, d3= 2500)
E23=RunEvaporation(Gsoil=3, Ksoil = 0.005, d1=300, d2=1200, d3= 2500)


par(mar=c(4,4,1,1))
plot(x=1:length(E11$EsoilCum), y=E11$EsoilCum, ylim=c(0,15), ylab="Cummulative soil evaporation (mm)", xlab="day", type="l", axes=FALSE)
axis(1, at=1:length(E11$EsoilCum), cex.axis=0.7)
axis(2)
points(x=1:length(E11$EsoilCum), y=E11$EsoilCum, pch=1)
lines(x=1:length(E12$EsoilCum), y=E12$EsoilCum, lty=2)
points(x=1:length(E12$EsoilCum), y=E12$EsoilCum, pch=1)
lines(x=1:length(E13$EsoilCum), y=E13$EsoilCum, lty=3)
points(x=1:length(E13$EsoilCum), y=E13$EsoilCum, pch=1)
lines(x=1:length(E21$EsoilCum), y=E21$EsoilCum, lty=1)
points(x=1:length(E21$EsoilCum), y=E21$EsoilCum, pch=2)
lines(x=1:length(E22$EsoilCum), y=E22$EsoilCum, lty=2)
points(x=1:length(E22$EsoilCum), y=E22$EsoilCum, pch=2)
lines(x=1:length(E23$EsoilCum), y=E23$EsoilCum, lty=3)
points(x=1:length(E23$EsoilCum), y=E23$EsoilCum, pch=2)
legend("topleft", lty=rep(1:3,2), pch=c(1,1,1,2,2,2), legend=c("Gsoil = 1 Ksoil = 0.05", 
                                    "Gsoil = 2 Ksoil = 0.05", 
                                    "Gsoil = 3 Ksoil = 0.05",
                                    "Gsoil = 1 Ksoil = 0.005", 
                                    "Gsoil = 2 Ksoil = 0.005", 
                                    "Gsoil = 3 Ksoil = 0.005"), cex = 0.7, bty="n")


