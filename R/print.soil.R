print.soil<-function(x,...) {
  #Depth
  cat(paste("Soil depth (mm):", round(x$SoilDepth, digits=0),
            "  Rock layer depth (mm):", round(x$RockLayerDepth, digits=0),"\n"))
  #Soil parameters related to texture
  if(x$dVec[1]>0) cat(paste("\nTopsoil \n    clay (%):", round(x$clay[1]),"sand (%):", round(x$sand[1]),
                            "[", x$usda_Type[1],"]\n    Rock fragment content (%):", round(x$rfc[1]),"Macroporosity (%):", round(x$macro[1]*100),  
                            "\n    Theta FC (%):", round(100*x$Theta_FC[1]),"Vol FC (mm):", round(x$Water_FC[1]),"\n"))
  if(x$dVec[2]>0) cat(paste("\nSubsoil\n    clay (%):", round(x$clay[2]),"sand (%):", round(x$sand[2]),
                            "[", x$usda_Type[2], "]\n    Rock fragment content (%):", round(x$rfc[2]),"Macroporosity (%):", round(x$macro[2]*100),
                            "\n    Theta FC (%):", round(100*x$Theta_FC[2]),"Vol FC (mm):", round(x$Water_FC[2]),"\n"))
  if(x$dVec[3]>0) cat(paste("\nRock layer\n    clay (%):", round(x$clay[3]),"sand (%):", round(x$sand[3]),
                            "[", x$usda_Type[3], "]\n    Rock fragment content (%):", round(x$rfc[3]),"Macroporosity (%):", round(x$macro[3]*100),
                            "\n    Theta FC (%):", round(100*x$Theta_FC[3]),"Vol FC (mm):", round(x$Water_FC[3]),"\n"))
  cat(paste("\nTotal soil water holding capacity (mm): ", round(sum(x$Water_FC), digits=0),"\n"))  
  cat(paste("\nCurrent Vol1 (mm):", round(x$W[1]*x$Water_FC[1]),"Vol2 (mm):", round(x$W[2]*x$Water_FC[2]),"Vol3 (mm):", round(x$W[3]*x$Water_FC[3]),"\n"))
}