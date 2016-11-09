#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export("soil.theta2psi")]]
double theta2psi(double clay, double sand, double theta) {
  double A = 100 * exp(-4.396 - (0.0715*clay)-(0.0004880*pow(sand,2)) - (0.00004285*pow(sand,2)*clay));
  double B = -3.140 - (0.00222*pow(clay,2)) - (0.00003484*pow(sand,2)*(clay));
  double psi = -1*A*pow(theta,B);
  if(psi < -40000) psi = -40000;
  if(theta==0.0) psi = -40000;
  return(psi);
}
// [[Rcpp::export("soil.psi2theta")]]
double psi2theta(double clay, double sand, double psi) {
  double A = 100 * exp(-4.396 - (0.0715*clay)-(0.0004880*pow(sand,2)) - (0.00004285*pow(sand,2)*clay));
  double B = -3.140 - (0.00222*pow(clay,2)) - (0.00003484*pow(sand,2)*(clay));
  return(pow(std::abs(psi)/A, 1.0/B));
}

// [[Rcpp::export("soil.USDAType")]]
String soilUSDAType(double clay, double sand) {
  double silt = 100 - clay - sand;
  if((silt+1.5*clay)<15) return("Sand");
  else if(((silt+1.5*clay)>=15) & ((silt + 2.0*clay)<30)) return("Loamy sand");
  else if(((clay>=7) & (clay<20) & (sand>52) & ((silt + 2.0*clay)>=30)) | ((clay < 7) & (silt < 50) & ((silt + 2.0*clay)>=30))) return("Sandy loam");
  else if(((clay>=7) & (clay<27)) & ((silt>=28) & (silt<50)) & (sand<=52)) return("Loam");
  else if(((silt>=50) & ((clay>=12) & (clay<27))) | ((silt>=50) & (silt<80) & (clay <12))) return("Silt loam");
  else if((silt>=80) & (clay<12)) return("Silt");
  else if(((clay>=20) & (clay<35)) & (silt<28) & (sand>45)) return("Sandy clay loam");
  else if(((clay>=27) & (clay<40)) & ((sand>20) & (sand<=45))) return("Clay loam");
  else if(((clay>=27) & (clay<40)) & (sand<=20)) return("Silty clay loam");
  else if((clay>=35) & (sand>45)) return("Sandy clay");
  else if((clay>=40) & (silt>=40)) return("Silty clay");
  else if((clay>=40) & (sand<=45) &(silt<40)) return("Clay");
  return("Unknown");
}




/* 
 * Leij, F.J., Alves, W.J., Genuchten, M.T. Van, Williams, J.R., 1996. The UNSODA Unsaturated Soil Hydraulic Database User’s Manual Version 1.0.
 * Textural parameters (1 MPa = 0.00009804139432 cm)
 */
// [[Rcpp::export("soil.vanGenuchtenParams")]]
NumericVector vanGenuchtenParams(String soilType) {
  NumericVector vg(3,NA_REAL);
  if(soilType=="Sand") {vg[0]=1478.967; vg[1]=2.68;                 vg[2] = 80000000000;}
  else if(soilType=="Loamy sand") {vg[0]=1264.772; vg[1]=2.28;      vg[2] = 80000000000;}
  else if(soilType=="Sandy loam") {vg[0]=764.983; vg[1]=1.89;       vg[2] = 80000000000;}
  else if(soilType=="Loam") {vg[0]=367.1918; vg[1]=1.56;            vg[2] = 80000000000;}
  else if(soilType=="Silt") {vg[0]=163.1964; vg[1]=1.37;            vg[2] = 80000000000;}
  else if(soilType=="Silt loam") {vg[0]=203.9955; vg[1]=1.41;       vg[2] = 400000000;}
  else if(soilType=="Sandy clay loam") {vg[0]=601.7866; vg[1]=1.48; vg[2] = 80000000000;}
  else if(soilType=="Clay loam") {vg[0]=193.7957; vg[1]=1.31;       vg[2] = 400000000;}
  else if(soilType=="Silty clay loam") {vg[0]=101.9977; vg[1]=1.23; vg[2] = 400000000;}
  else if(soilType=="Sandy clay") {vg[0]=275.3939; vg[1]=1.23;      vg[2] = 40000000;}
  else if(soilType=="Silty clay") {vg[0]=50.99887; vg[1]=1.09;      vg[2] = 40000000;}
  else if(soilType=="Clay") {vg[0]=81.59819; vg[1]=1.09;            vg[2] = 40000000;}
  return(vg);
}

// [[Rcpp::export("soil")]]
List soil(List SoilParams, NumericVector W = NumericVector::create(1.0,1.0,1.0)) {
    double SoilDepth = SoilParams["SoilDepth"];
  SoilDepth = std::min(SoilDepth,4000.0);
  double RockLayerDepth = SoilParams["RockLayerDepth"];
  RockLayerDepth = std::max(RockLayerDepth, SoilDepth);
  //TOPSOIL LAYER
  double d1 = std::min(SoilDepth, 300.0); //0-300 mm
  //SUBSOIL LAYER
  double d2 = std::max(0.0, std::min(SoilDepth-300.0,3700.0)); //300-4000 mm
  //ROCK LAYER
  double d3 = std::max(0.0,RockLayerDepth - (d1+d2));  //From SoilDepth down to RockLayerDepth
  
  NumericVector dVec = NumericVector::create(d1,d2,d3);
  
    //Soil parameters related to texture
  NumericVector clay = NumericVector::create(0.0,0.0,0.0);
  NumericVector sand = NumericVector::create(0.0,0.0,0.0);
  NumericVector macro = NumericVector::create(0.0,0.0,0.0);
  NumericVector rfc = NumericVector::create(0.0,0.0,0.0);

  clay[0] = SoilParams["TS_clay"];
  sand[0] = SoilParams["TS_sand"];
  macro[0] = SoilParams["TS_macro"];
  rfc[0] = SoilParams["TS_rfc"];
  if(d2 > 0.0) {
    clay[1] = SoilParams["SS_clay"];
    sand[1] = SoilParams["SS_sand"];
    macro[1] = SoilParams["SS_macro"];
    rfc[1] = SoilParams["SS_rfc"];
  }
  if(d3 > 0.0) {
    clay[2] = SoilParams["RL_clay"];
    sand[2] = SoilParams["RL_sand"];
    macro[2] = SoilParams["RL_macro"];
    rfc[2] = SoilParams["RL_rfc"];
  }


  NumericVector SoilPRFC = NumericVector::create(1.0-(rfc[0]/100.0),
                            1.0-(rfc[1]/100.0), 1.0-(rfc[2]/100.0));
  NumericVector Theta_FC = NumericVector::create(psi2theta(clay[0], sand[0], -33.0),
                    psi2theta(clay[1], sand[1], -33.0),psi2theta(clay[2], sand[2], -33.0));
  NumericVector Water_FC = NumericVector::create(d1*Theta_FC[0]*SoilPRFC[0],
                    d2*Theta_FC[1]*SoilPRFC[1], d3*Theta_FC[2]*SoilPRFC[2]);
  
  CharacterVector usda_Type = CharacterVector::create(soilUSDAType(clay[0],sand[0]),
                                                      soilUSDAType(clay[1],sand[1]),
                                                      soilUSDAType(clay[2],sand[2]));
  
  NumericVector psi = NumericVector::create(theta2psi(clay[0], sand[0], W[0]*Theta_FC[0]),
                        theta2psi(clay[1], sand[1], W[1]*Theta_FC[1]),
                        theta2psi(clay[2], sand[2], W[2]*Theta_FC[2]));    

  NumericVector VG1 = vanGenuchtenParams(usda_Type[0]);
  NumericVector VG2 = vanGenuchtenParams(usda_Type[1]);
  NumericVector VG3 = vanGenuchtenParams(usda_Type[2]);
  NumericVector VG_alpha = NumericVector::create(VG1[0],VG2[0],VG3[0]);
  NumericVector VG_n = NumericVector::create(VG1[1],VG2[1],VG3[1]);
  NumericVector VG_ksmax = NumericVector::create(VG1[2],VG2[2],VG3[2]);
  double Ssoil = Water_FC[0];
  List l = List::create(_["SoilDepth"] = SoilDepth, _["RockLayerDepth"] = RockLayerDepth,
                      _["W"] = W, _["psi"] = psi, _["Ksoil"] = SoilParams["Ksoil"], _["Gsoil"] = SoilParams["Gsoil"],
                      _["Ssoil"] = Ssoil,
                      _["dVec"] = dVec,
                      _["sand"] = sand, _["clay"] = clay,
                      _["usda_Type"] = usda_Type,
                      _["VG_alpha"] = VG_alpha,_["VG_n"] = VG_n,_["VG_ksmax"] = VG_ksmax, 
                      _["macro"] = macro, _["rfc"] = rfc,
                      _["Theta_FC"] = Theta_FC, _["Water_FC"] = Water_FC);
  l.attr("class") = CharacterVector::create("soil","list");
  return(l);
}
/*** R
usdaV = numeric(5052)
sandV = numeric(5052)
siltV = numeric(5052)
clayV = numeric(5052)
alphaV = numeric(5052)
nV = numeric(5052)
cnt=0
for(clay in 0:100) for(sand in (clay+1):100) {
  cnt = cnt+1
  usdaV[cnt] = soil.soilUSDAType(clay,sand)
  clayV[cnt] = clay
  sandV[cnt] = sand
  siltV[cnt] = 100- clay-sand
  p = soil.vanGenuchtenParams(usdaV[cnt])
  alphaV[cnt] = p[1]
  nV[cnt] = p[2]
}
table(usdaV)
t = data.frame(clayV, siltV, sandV, usdaV, alphaV, nV)
print(head(t))  
t[usdaV=="Unknown",]
sum(is.na(t$alphaV))
*/
