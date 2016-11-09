swbpoints<-function(y, SpParams, meteo, dates = NULL, 
                    control = defaultControl(), summaryFunction=NULL, args=NULL) {
  
  #Check input
  if(!inherits(meteo,"data.frame") && 
     !inherits(meteo,"SpatialPointsMeteorology") && 
     !inherits(meteo,"SpatialPointsDataFrame")) 
    stop("'meteo' has to be of class 'data.frame', 'SpatialPointsMeteorology' or 'SpatialPointsDataFrame'.")
  
  #Get spatial object properties
  sp = SpatialPoints(y@coords, y@proj4string)
  IDs = names(y@forestlist)
  if(inherits(meteo,"SpatialPoints")) {
    if(sum(y@coords == meteo@coords)!=2*nrow(y@coords)) stop("Coordinates of 'y' and 'meteo' must be the same.")
  }
  # .checkSpatialInput(y, SpParams, dates, Temperature, Rainfall, Radiation)
  
  longlat = spTransform(y,CRS("+proj=longlat"))
  latitude = longlat@coords[,2]
  
  reslist = vector("list",length(IDs))
  names(reslist) = IDs
  for(i in 1:length(IDs)) {
    id = IDs[i]
    if(!is.na(id)) {
      if(control$verbose) cat(paste("\nProcessing ",id," at (", y@coords[i,1],", ", y@coords[i,2],")", sep=""))
      else cat(".")
      if(control$verbose) cat(" - Meteorology")
      if(inherits(meteo,"data.frame")) met = meteo
      else if(inherits(meteo,"SpatialPointsMeteorology")) {
        met = meteo@data[[i]]
      } else {
        f = paste(meteo@data$dir[i], meteo@data$filename[i],sep="/")
        if(!file.exists(f)) stop(paste("Meteorology file '", f,"' does not exist!", sep=""))
        met = readmeteorologypoint(f)
      }
      if(!is.null(dates)) met = met[as.character(dates),] #subset dates
      yid = y@forestlist[[i]] 
      soil = y@soillist[[i]]
      if(control$verbose) cat(" - Soil water balance")
      swbInput = forest2swbInput(yid, SpParams, soil$dVec, petMode = control$petMode, hydraulicMode = control$hydraulicMode)
      S<-swb(swbInput, soil, meteo=met, petMode = control$petMode, hydraulicMode = control$hydraulicMode, 
             elevation = y@data$elevation[i], latitude = latitude[i],
            verbose = FALSE)      
      if(!is.null(summaryFunction)){
        #Fill argument list for  summary function
        argList = list(object=S)
        if(!is.null(args)) {
          for(j in 1:length(args)) argList[names(args)[j]]=args[j]
        }
        reslist[[i]] = do.call(summaryFunction, args=argList)
      } else {
        reslist[[i]] = S
      }
    } else {
      if(control$verbose) cat("\nEmpty cell at (", y@coords[i,1],", ", y@coords[i,2],")", sep="")
      else cat(".")
    }
  }
  l = list(sp = sp, result = reslist)
  class(l) = c("swbpoints","list")
  return(l)
}