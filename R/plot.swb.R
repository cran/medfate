plot.swb<-function(x, yearAxis=FALSE, type="PET_Precipitation", xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,...) {
  dates = as.Date(rownames(x$DailyBalance))
  DailyBalance = x$DailyBalance
  SoilWaterBalance = x$SoilWaterBalance
  nlayers = 3
  TYPES = c("PET_Precipitation","PET_NetPrec","ET","Psi","Theta","Vol","PlantStress","Export")
  type = match.arg(type,TYPES)  
  numDays = length(dates)
  numYears = round(numDays/365)
  firstYear=as.numeric(format(dates[1],"%Y"))
  
  plotAxes<-function(){
    if(!yearAxis) axis(1)
    else {
      axis(1, at = (0:numYears)*365, labels=FALSE)
      axis(1, at = -182+365*(1:numYears), tick = FALSE, line=FALSE, labels=firstYear:(firstYear+numYears-1))
    }
    axis(2)    
  }
  mnp = max(DailyBalance$Precipitation)
  if(is.null(xlab)) xlab = ifelse(yearAxis,"Year", "Day")  
  if(type=="PET_Precipitation") {
    if(is.null(ylab)) ylab = "mm water"
    if(!is.null(xlim)) span = xlim[1]:xlim[2]
    else span = 1:numDays
    if(is.null(ylim)) ylim = c(0,mnp)
    barplot(DailyBalance$Precipitation[span], ylim=ylim, col="black",space=0, ylab=ylab, 
            xlab=xlab, axes=FALSE)
    plotAxes()
    lines(1:length(span), DailyBalance$PET[span], col="gray")    
    legend("topleft", bty="n", col=c("black","gray"),lty=c(1,1), lwd=2,
           legend=c("Precipitation","PET"))
    
  } else if(type=="PET_NetPrec") {
    if(is.null(ylab)) ylab = "mm water"    
    if(!is.null(xlim)) span = xlim[1]:xlim[2]
    else span = 1:numDays
    if(is.null(ylim)) ylim = c(0,mnp)
    barplot(DailyBalance$NetPrec[span], ylim=ylim, col="black",space=0, ylab=ylab, 
            xlab=xlab, axes=FALSE)
    plotAxes()
    lines(1:length(span), DailyBalance$PET[span], col="gray")    
    legend("topleft", bty="n", col=c("black","gray"),lty=c(1,1), lwd=2,
           legend=c("NetPrec","PET"))        
  } else if(type=="ET") {
    if(is.null(ylab)) ylab = "mm water"
    if(is.null(ylim)) ylim = c(0,max(DailyBalance$Etot))
    plot(1:numDays, DailyBalance$Etot, ylim=ylim, type="l", ylab=ylab, 
         xlab=xlab, xlim=xlim,frame=FALSE, col="black", axes=FALSE, lwd=2)
    plotAxes()
    lines(1:numDays, DailyBalance$Eplanttot, col="gray", lty=2, lwd=1.5)
    lines(1:numDays, DailyBalance$Esoil, col="black", lty=3, lwd=1.5)
    legend("topleft", bty="n", col=c("black","gray","black"),lty=c(1,2,3), lwd=c(2,1.5,1.5),
           legend=c("Total evaporation","Plant ET","Bare soil E"))
    
  } else if(type=="Psi") {
    if(is.null(ylab)) ylab = "Psi (kPa)"    
    if(is.null(ylim)) ylim =c(-8000,0)
    plot(1:numDays, SoilWaterBalance$psi1, lwd=1.5,
         ylim=ylim, type="l", ylab=ylab, xlab=xlab, xlim=xlim,
         frame=FALSE, axes=FALSE)
    plotAxes()
    lines(1:numDays, SoilWaterBalance$psi2,lwd=1.5,lty=2)
    lines(1:numDays, SoilWaterBalance$psi3, lwd=1.5,lty=3)
    legend("bottomleft", bty="n", lty=c(1,2,3,4,5), col="black", lwd=1.5,
           legend=paste("Layer", 1:nlayers))
    
  } else if(type=="Theta") {
    W1 = SoilWaterBalance$W1
    W2 = SoilWaterBalance$W2
    W3 = SoilWaterBalance$W3
    if(is.null(ylab)) ylab = "% field capacity"
    if(is.null(ylim)) ylim = c(0,100)
    plot(1:numDays, W1*100, ylim= ylim,lwd=1.5, type="l", ylab=ylab, xlim=xlim,
         xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    lines(1:numDays, W2*100, lwd=1.5,lty=2)
    lines(1:numDays, W3*100, lwd=1.5,lty=3)
    legend("bottomleft", bty="n", col="black",lty=1:nlayers, lwd=1.5,
           legend=paste("Layer", 1:nlayers))
  } else if(type=="Vol") {
    ML1 = SoilWaterBalance$ML1
    ML2 = SoilWaterBalance$ML2
    ML3 = SoilWaterBalance$ML3
    MLTot = SoilWaterBalance$MLTot
    if(is.null(ylim)) ylim =c(0,max(MLTot)*1.3)
    if(is.null(ylab)) ylab = "mm soil water"
    plot(1:numDays, MLTot, ylim=ylim, lwd=2, type="l",xlim=xlim,
         ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()
    lines(1:numDays, ML1, lwd=1.5,lty=1)
    lines(1:numDays, ML2, lwd=1.5,lty=2)
    lines(1:numDays, ML3, lwd=1.5,lty=3)
    legend("topleft", bty="n", col="black",lty=c(1,1:nlayers), lwd=c(2,rep(1.5,5)),
           legend=c("Total",paste("Layer", 1:nlayers)))
    
  } else if(type=="PlantStress") {
    if(is.null(ylab)) ylab = "Drought stress"
    if(is.null(ylim)) ylim = c(0,1)
    matplot(1:numDays, x$PlantStress, ylim = ylim, lwd=1, type="l", xlim=xlim,
            ylab=ylab, xlab=xlab, frame=FALSE, axes=FALSE)
    plotAxes()        
  } else if(type=="Export") {
    if(is.null(ylab)) ylab = "mm water"    
    mnp = max(DailyBalance$DeepDrainage+DailyBalance$Runoff)    
    if(is.null(ylim)) ylim = c(0,mnp)
    plot(1:numDays, DailyBalance$DeepDrainage+DailyBalance$Runoff, ylim=ylim, col="black", type="l", 
         ylab=ylab, xlab=xlab, xlim=xlim,
         frame=FALSE, axes=FALSE)
    lines(1:numDays, DailyBalance$DeepDrainage, col="blue")
    lines(1:numDays, DailyBalance$Runoff, col="red")
    plotAxes()
    legend("topright", bty="n", col=c("black","blue","red"),lty=c(1,1,1), lwd=c(1.5,1,1),
           legend=c("DD+R","Deep drainage (DD)","Runoff (R)"))        
  }
}
