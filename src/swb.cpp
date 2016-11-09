#include <numeric>
#include "lightextinction.h"
#include "hydraulics.h"
#include "PenmanMonteith.h"
#include "soil.h"
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(".er")]]
NumericVector er(IntegerVector DOY, double ERconv=0.05, double ERsyn = 0.2){
  int nDays = DOY.size();
  NumericVector ER=rep(0.0,nDays);
  for(int i=0;i<nDays;i++){
    if((DOY[i]<=120)|(DOY[i]>=335)) {
      ER[i] = ERsyn;
    } else {
      ER[i] = ERconv;
    }
  }
  return(ER);
  
}
// [[Rcpp::export(".gdd")]]
NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0){
  int nDays = Temp.size();
  NumericVector GDD(nDays);
  double cum = 0.0;
  for(int i=0;i<nDays;i++){
    if((Temp[i]-Tbase < 0.0) & (DOY[i]>180)) {
      cum = 0.0;
    } else {
      if(Temp[i]-Tbase>0.0) cum = cum + (Temp[i]-Tbase);
    }
    GDD[i] = cum;
    if(DOY[i] >= 365) cum = 0.0;
  }
  return(GDD);
}

// [[Rcpp::export("swb.SoilEvaporation")]]
double soilevaporation(double DEF,double PETs, double Gsoil){
  double t = pow(DEF/Gsoil, 2.0);
  double Esoil = 0.0;
  Esoil = std::min(Gsoil*(sqrt(t+1)-sqrt(t)), PETs);
  return(Esoil);
}
// [[Rcpp::export(".infiltrationDay")]]
double infiltrationDay(double NetPrec, double Ssoil) {
  double I = 0;
  if(NetPrec>0.2*Ssoil) {
    I = NetPrec-(pow(NetPrec-0.2*Ssoil,2.0)/(NetPrec+0.8*Ssoil));
  } else {
    I = NetPrec;
  }
  return(I);
}

/**
 * Whole-plant conductance function
 */
NumericVector Psi2K(double psi, NumericVector psi_extr, int ws) {
  int n = psi_extr.size();
  NumericVector k(n);
  for(int i=0; i<n; i++) {
    k[i] = exp(-0.6931472*pow((std::abs(psi)/std::abs(psi_extr[i])),ws));
  }
  return k;
}
/**
 * Inverse of the whole-plant conductance function. Used to obtain the 'average' soil water
 * potential perceived by each plant cohort.
 */
NumericVector K2Psi(NumericVector K, NumericVector psi_extr, int ws) {
  int n = psi_extr.size();
  NumericVector psi(n);
  for(int i=0; i<n; i++) {
    psi[i] = psi_extr[i]*pow(log(K[i])/(-0.6931472),1.0/ws);
    if(psi[i]>0.0) psi[i] = -psi[i]; //Usually psi_extr is a positive number
  }
  return psi;
}

// [[Rcpp::export(".interceptionGashDay")]]
double interceptionGashDay(double Precipitation, double Cm, double p, double ER=0.05) {
    double I = 0.0;
    double PG = (-Cm/(ER*(1.0-p)))*log(1.0-ER); //Precipitation need to saturate the canopy
    if(Cm==0.0 || p==1.0) PG = 0.0; //Avoid NAs
    if(Precipitation>PG) {
      I = (1-p)*PG + (1-p)*ER*(Precipitation-PG);
    } else {
      I = (1-p)*Precipitation;
    }
    return(I);
}


// Soil water balance with PET as input
// [[Rcpp::export(".swbDay1")]]
List swbDay1(DataFrame x, List soil, double gdd, double pet, double rain, double er, double runon=0.0, String hydraulicMode = "Simple", bool verbose = false) {
  NumericVector LAI = Rcpp::as<Rcpp::NumericVector>(x["LAI"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(x["H"]);
  NumericVector V1 = Rcpp::as<Rcpp::NumericVector>(x["V.1"]);
  NumericVector V2 = Rcpp::as<Rcpp::NumericVector>(x["V.2"]);
  NumericVector V3 = Rcpp::as<Rcpp::NumericVector>(x["V.3"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(x["Sgdd"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(x["k"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(x["g"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(x["CR"]);
  NumericVector transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);

  NumericVector W = soil["W"];
  NumericVector psi = soil["psi"];
  NumericVector dVec = soil["dVec"];
  NumericVector Theta_FC = soil["Theta_FC"];
  NumericVector Water_FC = soil["Water_FC"];
  NumericVector macro = soil["macro"];
  NumericVector clay = soil["clay"];
  NumericVector sand = soil["sand"];
  double d1 = dVec[0], d2 = dVec[1], d3 = dVec[2];
  
  int numCohorts = LAI.size();
  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  NumericVector LAIphe(numCohorts);
  NumericVector Phe = pmin(pmax(gdd/Sgdd,0.0),1.0);
  double s = 0.0, LAIcell = 0.0, Cm = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(Sgdd[c]==0.0) Phe[c]=1.0;
    LAIphe[c] = LAI[c]*Phe[c]; //LAI modified by phenology
    s += (kPAR[c]*LAIphe[c]);
    LAIcell += LAIphe[c];
    Cm += LAIphe[c]*gRainIntercept[c];
  }
  NumericVector CohASWRF = cohortAbsorbedSWRFraction(LAIphe,  H, CR, kPAR);
  NumericVector CohPAR = parcohortC(H, LAIphe, kPAR, CR)/100.0;
  double LgroundPAR = exp((-1)*s);
  double LgroundSWR = 1.0 - std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  
  //Hydrologic input
  double NetPrec = 0.0, Infiltration= 0.0, Runoff= 0.0, DeepDrainage= 0.0;
  if(rain>0.0) {
    //Interception
    NetPrec = rain - interceptionGashDay(rain,Cm,LgroundPAR,er);
    //Net Runoff and infiltration
    Infiltration = infiltrationDay(NetPrec+runon, soil["Ssoil"]);
    Runoff = (NetPrec+runon) - Infiltration;
    //Input of the first soil layer is infiltration
    double VI = Infiltration;
    double Wn;
    //Update topsoil layer
    if(VI>0) {
      Wn = W[0]*Water_FC[0] + VI*(1.0-macro[0]); //Update water volume
      VI = VI*macro[0] + std::max(Wn - Water_FC[0],0.0); //Add the excess to the infiltrating water (saturated flow)
      W[0] = std::min(std::max(0.0,Wn/Water_FC[0]),1.0); //Update theta
    }
    //Update subsoil layer
    if((d2>0) & (VI>0)){
      Wn = W[1]*Water_FC[1] + VI*(1.0-macro[1]); //Update water volume
      VI = VI*macro[1] + std::max(Wn - Water_FC[1], 0.0); //Add the excess to the percolating water (saturated flow)
      W[1] = std::min(std::max(0.0,Wn/Water_FC[1]),1.0); //Update theta
    }
    // Update rocky layer
    if((d3>0) & (VI>0)){
      Wn = W[2]*Water_FC[2] + VI*(1.0-macro[2]);
      DeepDrainage = VI*macro[2]+std::max(Wn - Water_FC[2],0.0); //Set the excess as deep drainage
      W[2] = std::min(std::max(0.0,Wn/Water_FC[2]),1.0); //Update theta
    } else if(VI>0) {
      DeepDrainage = VI;
    }
  }

  psi[0] = theta2psi(clay[0], sand[0], W[0]*Theta_FC[0]);
  psi[1] = theta2psi(clay[1], sand[1], W[1]*Theta_FC[1]);
  psi[2] = theta2psi(clay[2], sand[2], W[2]*Theta_FC[2]);


  //Proportion of transpiration that absorbed by each plant cohort (old version)
  // NumericVector PP = CohLight*LAIphe;
  // NumericVector f = PP/std::accumulate(PP.begin(),PP.end(),0.0); 
  // if(LAIcell==0.0) f = rep(0.0,numCohorts); //Avoids NaN values

  //Apply fractions to potential evapotranspiration
  //Maximum canopy transpiration
  //    Tmax = PET[i]*(-0.006*pow(LAIcell[i],2.0)+0.134*LAIcell[i]+0.036); //From Granier (1999)
  double Tmax = pet*(-0.006*pow(LAIcell,2.0)+0.134*LAIcell); //From Granier (1999)
  double PETsoil = pet*LgroundSWR;

  //Fraction of Tmax attributed to each plant cohort
  double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  NumericVector TmaxCoh(numCohorts);
  if(pabs>0.0) TmaxCoh = Tmax*(CohASWRF/pabs);
  
  //Actual plant transpiration
  NumericVector K1(numCohorts), K2(numCohorts), K3(numCohorts);
  NumericVector Eplant1Coh(numCohorts),Eplant2Coh(numCohorts),Eplant3Coh(numCohorts);
  NumericVector PlantPsi(numCohorts, NA_REAL);
  
  if(hydraulicMode =="Simple") {
    int WeibullShape=3;
    NumericVector Psi_Extract = Rcpp::as<Rcpp::NumericVector>(x["Psi_Extract"]);
    K1 = Psi2K(psi[0], Psi_Extract, WeibullShape);
    K2 = Psi2K(psi[1], Psi_Extract, WeibullShape);
    K3 = Psi2K(psi[2], Psi_Extract, WeibullShape);
    Eplant1Coh = pmax(TmaxCoh*K1*V1,0.0);
    Eplant2Coh = pmax(TmaxCoh*K2*V2,0.0);
    Eplant3Coh = pmax(TmaxCoh*K3*V3,0.0);
  } else if(hydraulicMode=="Sperry") {
    NumericVector VG_ksmax = soil["VG_ksmax"];
    NumericVector VG_n = soil["VG_n"];
    NumericVector VG_alpha = soil["VG_alpha"];
    NumericVector VC_kxmax = x["VC_kxmax"];
    NumericVector VC_c = x["VC_c"];
    NumericVector VC_d = x["VC_d"];
    for(int c=0;c<numCohorts;c++){
      double TmaxperLAI = TmaxCoh[c]/LAIphe[c];
      if(LAIphe[c]<=0.0) TmaxperLAI = 0.0;
      if(TmaxperLAI>0.0) { 
        NumericVector sperry1 = regulatedPsiTwoElements(TmaxperLAI, psi[0]/1000.0, VG_ksmax[0], VC_kxmax[c], 
                                                        VG_n[0], VG_alpha[0], VC_c[c], VC_d[c], 0.1, -10.0);
        NumericVector sperry2 = regulatedPsiTwoElements(TmaxperLAI, psi[1]/1000.0, VG_ksmax[1], VC_kxmax[c],
                                                        VG_n[1], VG_alpha[1], VC_c[c], VC_d[c], 0.1, -10.0);
        NumericVector sperry3 = regulatedPsiTwoElements(TmaxperLAI, psi[2]/1000.0, VG_ksmax[2], VC_kxmax[c],
                                                        VG_n[2], VG_alpha[2], VC_c[c], VC_d[c], 0.1, -10.0);
        
        //Actual transpirated water from each soil layer depending on the amount of roots
        Eplant1Coh[c] = (sperry1[3]*LAIphe[c])*V1[c];
        Eplant2Coh[c] = (sperry2[3]*LAIphe[c])*V2[c];
        Eplant3Coh[c] = (sperry3[3]*LAIphe[c])*V3[c];
        
        //Whole plant relative conductance can be calculated from the supply function
        //but we calculate it from the comparison of TmaxCoh and actual transpiration
        //for comparison with the simple hydraulic model
        if((V1[c]>0.0)) K1[c] = sperry1[3]/TmaxperLAI;
        else K1[c] = 0.0;
        if((V2[c]>0.0)) K2[c] = sperry2[3]/TmaxperLAI;
        else K2[c] = 0.0;
        if((V3[c]>0.0)) K3[c] = sperry3[3]/TmaxperLAI;
        else K3[c] = 0.0;
        
        //Average of water potential (kPa) according to proportion of roots (need to improve!)
        PlantPsi[c] = 1000.0*((sperry1[1]*V1[c])+(sperry2[1]*V2[c])+(sperry3[1]*V3[c]));
      } else {
        Eplant1Coh[c] = 0.0;
        Eplant2Coh[c] = 0.0;
        Eplant3Coh[c] = 0.0;
        PlantPsi[c] = 0.0;
        K1[c] = 0.0;
        K2[c] = 0.0;
        K3[c] = 0.0;
      }
    }
  }
  NumericVector Eplant = Eplant1Coh+Eplant2Coh+Eplant3Coh;

  NumericVector EplantVec = NumericVector::create(std::accumulate(Eplant1Coh.begin(),Eplant1Coh.end(),0.0),
                            std::accumulate(Eplant2Coh.begin(),Eplant2Coh.end(),0.0),
                            std::accumulate(Eplant3Coh.begin(),Eplant3Coh.end(),0.0));

  //Evaporation from bare soil
  double Gsoil = soil["Gsoil"];
  double Ksoil = soil["Ksoil"];
  double Esoil = soilevaporation((Water_FC[0]*(1.0 - W[0])), PETsoil, Gsoil);

  //Exponential decay to divide bare soil evaporation among layers
  NumericVector EsoilVec = NumericVector::create(Esoil*(1.0-exp(-Ksoil*d1)),
                                                 Esoil*(exp(-Ksoil*d1)-exp(-Ksoil*(d1+d2))),
                                                 Esoil*(exp(-Ksoil*(d1+d2))));

  //Apply decrease in soil layers
  W[0] = std::max(W[0]-(EplantVec[0]+EsoilVec[0])/Water_FC[0],0.0);
  if(d2>0) W[1] = std::max(std::min(W[1]-(EplantVec[1]+EsoilVec[1])/Water_FC[1],1.0),0.0);
  if(d3>0) W[2] = std::max(std::min(W[2]-(EplantVec[2]+EsoilVec[2])/Water_FC[2],1.0),0.0);

  NumericVector DDS = (Phe*(V1*(1.0 - K1)+V2*(1.0 - K2)+V3*(1.0 - K3)));

  //For comunication with photosynthesis
  for(int c=0;c<numCohorts;c++) transpiration[c] = Eplant[c];

  List l = List::create(_["NetPrec"] = NetPrec, _["Runon"] = runon, _["Infiltration"] = Infiltration, _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage,
                        _["LAIcell"] = LAIcell, _["Cm"] = Cm, _["Lground"] = LgroundPAR, _["Tmax"] = Tmax,
                        _["EsoilVec"] = EsoilVec, _["TmaxCoh"] = TmaxCoh, _["EplantVec"] = EplantVec, _["psiVec"] = psi,
                        _["EplantCoh"] = Eplant, _["psiCoh"] = PlantPsi, _["DDS"] = DDS);
  return(l);
}


// Soil water balance with PET calculated using Penman-Monteith
// [[Rcpp::export(".swbDay2")]]
List swbDay2(DataFrame x, List soil, double gdd, double tmin, double tmax, double rhmin, double rhmax, double rad, double wind, 
             double latitude, double elevation,  int J, 
             double rain, double er, double runon=0.0, 
             String hydraulicMode = "Simple") {
  
  
  NumericVector LAI = Rcpp::as<Rcpp::NumericVector>(x["LAI"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(x["H"]);
  NumericVector V1 = Rcpp::as<Rcpp::NumericVector>(x["V.1"]);
  NumericVector V2 = Rcpp::as<Rcpp::NumericVector>(x["V.2"]);
  NumericVector V3 = Rcpp::as<Rcpp::NumericVector>(x["V.3"]);
  NumericVector Sgdd = Rcpp::as<Rcpp::NumericVector>(x["Sgdd"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(x["k"]);
  NumericVector gRainIntercept = Rcpp::as<Rcpp::NumericVector>(x["g"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(x["CR"]);
  NumericVector transpiration = Rcpp::as<Rcpp::NumericVector>(x["Transpiration"]);
  
  NumericVector W = soil["W"];
  NumericVector psi = soil["psi"];
  NumericVector dVec = soil["dVec"];
  NumericVector Theta_FC = soil["Theta_FC"];
  NumericVector Water_FC = soil["Water_FC"];
  NumericVector macro = soil["macro"];
  NumericVector clay = soil["clay"];
  NumericVector sand = soil["sand"];
  
  double d1 = dVec[0], d2 = dVec[1], d3 = dVec[2];
  
  int numCohorts = LAI.size();
  //Determine whether leaves are out (phenology) and the adjusted Leaf area
  NumericVector LAIphe(numCohorts);
  NumericVector Phe = pmin(pmax(gdd/Sgdd,0.0),1.0);
  double s = 0.0, LAIcell = 0.0, Cm = 0.0, canopyHeight = 0.0;
  for(int c=0;c<numCohorts;c++) {
    if(Sgdd[c]==0.0) Phe[c]=1.0;
    LAIphe[c] = LAI[c]*Phe[c]; //LAI modified by phenology
    s += (kPAR[c]*LAIphe[c]);
    LAIcell += LAIphe[c];
    Cm += LAIphe[c]*gRainIntercept[c];
    if(canopyHeight<H[c]) canopyHeight = H[c];
  }
  //Light extinction
  NumericVector CohASWRF = cohortAbsorbedSWRFraction(LAIphe,  H, CR, kPAR);
  NumericVector CohPAR = parcohortC(H, LAIphe, kPAR, CR)/100.0;
  double LgroundPAR = exp((-1)*s);
  double LgroundSWR = 1.0 - std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  
  //Hydrologic input
  double NetPrec = 0.0, Infiltration= 0.0, Runoff= 0.0, DeepDrainage= 0.0;
  if(rain>0.0) {
    //Interception
    NetPrec = rain - interceptionGashDay(rain,Cm,LgroundPAR,er);
    //Net Runoff and infiltration
    Infiltration = infiltrationDay(NetPrec+runon, soil["Ssoil"]);
    Runoff = (NetPrec+runon) - Infiltration;
    //Input of the first soil layer is infiltration
    double VI = Infiltration;
    double Wn;
    //Update topsoil layer
    if(VI>0) {
      Wn = W[0]*Water_FC[0] + VI*(1.0-macro[0]); //Update water volume
      VI = VI*macro[0] + std::max(Wn - Water_FC[0],0.0); //Add the excess to the infiltrating water (saturated flow)
      W[0] = std::min(std::max(0.0,Wn/Water_FC[0]),1.0); //Update theta
    }
    //Update subsoil layer
    if((d2>0) & (VI>0)){
      Wn = W[1]*Water_FC[1] + VI*(1.0-macro[1]); //Update water volume
      VI = VI*macro[1] + std::max(Wn - Water_FC[1], 0.0); //Add the excess to the percolating water (saturated flow)
      W[1] = std::min(std::max(0.0,Wn/Water_FC[1]),1.0); //Update theta
    }
    // Update rocky layer
    if((d3>0) & (VI>0)){
      Wn = W[2]*Water_FC[2] + VI*(1.0-macro[2]);
      DeepDrainage = VI*macro[2]+std::max(Wn - Water_FC[2],0.0); //Set the excess as deep drainage
      W[2] = std::min(std::max(0.0,Wn/Water_FC[2]),1.0); //Update theta
    } else if(VI>0) {
      DeepDrainage = VI;
    }
  }
  
  psi[0] = theta2psi(clay[0], sand[0], W[0]*Theta_FC[0]);
  psi[1] = theta2psi(clay[1], sand[1], W[1]*Theta_FC[1]);
  psi[2] = theta2psi(clay[2], sand[2], W[2]*Theta_FC[2]);
  // Rcout<<psi[0]<<" "<<psi[1]<<" "<< psi[2]<<"\n";
  
  
  //Net radiation
  double Rn = netRadiation(latitude,elevation, J, tmin, tmax, rhmin,rhmax,rad);
 
  double condsum = 0.0;
  NumericVector RC_min = x["RC_min"];
  for(int c=0;c<numCohorts;c++) {
    condsum = condsum + (1.0/RC_min[c])*LAIphe[c];
  }
  //Maximum canopy transpiration using Penman-Monteith
  double Rncanopy = Rn*(1.0-LgroundSWR);
  double Tmax = PenmanMonteithPET((1.0/condsum), elevation, tmin, tmax, rhmin, rhmax, Rncanopy, wind);
  
  //Fraction of Tmax attributed to each plant cohort
  double pabs = std::accumulate(CohASWRF.begin(),CohASWRF.end(),0.0);
  NumericVector TmaxCoh(numCohorts);
  if(pabs>0.0) TmaxCoh = Tmax*(CohASWRF/pabs);
  
  //Potential soil evaporation
  double PETsoil = PenmanMonteithPET(200.0, elevation, tmin, tmax, rhmin, rhmax, Rn, wind);
  
  //Total PET (for comparison purposes only)
  double pet = PETsoil + Tmax;
  
  //Plant transpiration for each cohort
  NumericVector Eplant1Coh(numCohorts);
  NumericVector Eplant2Coh(numCohorts);
  NumericVector Eplant3Coh(numCohorts);
  NumericVector PlantPsi(numCohorts);
  NumericVector K1(numCohorts);
  NumericVector K2(numCohorts);
  NumericVector K3(numCohorts);
  if(hydraulicMode =="Simple") {
    int WeibullShape=3;
    NumericVector Psi_Extract = Rcpp::as<Rcpp::NumericVector>(x["Psi_Extract"]);
    K1 = Psi2K(psi[0], Psi_Extract, WeibullShape);
    K2 = Psi2K(psi[1], Psi_Extract, WeibullShape);
    K3 = Psi2K(psi[2], Psi_Extract, WeibullShape);
    Eplant1Coh = pmax(TmaxCoh*K1*V1,0.0);
    Eplant2Coh = pmax(TmaxCoh*K2*V2,0.0);
    Eplant3Coh = pmax(TmaxCoh*K3*V3,0.0);
  } else if(hydraulicMode =="Sperry") {
    NumericVector VG_ksmax = soil["VG_ksmax"];
    NumericVector VG_n = soil["VG_n"];
    NumericVector VG_alpha = soil["VG_alpha"];
    NumericVector VC_kxmax = x["VC_kxmax"];
    NumericVector VC_c = x["VC_c"];
    NumericVector VC_d = x["VC_d"];
    for(int c=0;c<numCohorts;c++){
      double TmaxperLAI = TmaxCoh[c]/LAIphe[c];
      if(LAIphe[c]<=0.0) TmaxperLAI = 0.0;
      if(TmaxperLAI>0.0) { 
        NumericVector sperry1 = regulatedPsiTwoElements(TmaxperLAI, psi[0]/1000.0, VG_ksmax[0], VC_kxmax[c], 
                                                        VG_n[0], VG_alpha[0], VC_c[c], VC_d[c], 0.1, -10.0);
        NumericVector sperry2 = regulatedPsiTwoElements(TmaxperLAI, psi[1]/1000.0, VG_ksmax[1], VC_kxmax[c],
                                                        VG_n[1], VG_alpha[1], VC_c[c], VC_d[c], 0.1, -10.0);
        NumericVector sperry3 = regulatedPsiTwoElements(TmaxperLAI, psi[2]/1000.0, VG_ksmax[2], VC_kxmax[c],
                                                        VG_n[2], VG_alpha[2], VC_c[c], VC_d[c], 0.1, -10.0);
        
        //Actual transpirated water from each soil layer depending on the amount of roots
        Eplant1Coh[c] = (sperry1[3]*LAIphe[c])*V1[c];
        Eplant2Coh[c] = (sperry2[3]*LAIphe[c])*V2[c];
        Eplant3Coh[c] = (sperry3[3]*LAIphe[c])*V3[c];
        
        //Whole plant relative conductance can be calculated from the supply function
        //but we calculate it from the comparison of TmaxCoh and actual transpiration
        //for comparison with the simple hydraulic model
        if((V1[c]>0.0)) K1[c] = sperry1[3]/TmaxperLAI;
        else K1[c] = 0.0;
        if((V2[c]>0.0)) K2[c] = sperry2[3]/TmaxperLAI;
        else K2[c] = 0.0;
        if((V3[c]>0.0)) K3[c] = sperry3[3]/TmaxperLAI;
        else K3[c] = 0.0;
        
        //Average of water potential (kPa) according to proportion of roots (need to improve!)
        PlantPsi[c] = 1000.0*((sperry1[1]*V1[c])+(sperry2[1]*V2[c])+(sperry3[1]*V3[c]));
      } else {
        Eplant1Coh[c] = 0.0;
        Eplant2Coh[c] = 0.0;
        Eplant3Coh[c] = 0.0;
        PlantPsi[c] = 0.0;
        K1[c] = 0.0;
        K2[c] = 0.0;
        K3[c] = 0.0;
      }
    }
  }
  NumericVector Eplant = Eplant1Coh+Eplant2Coh+Eplant3Coh;
    
  NumericVector EplantVec = NumericVector::create(std::accumulate(Eplant1Coh.begin(),Eplant1Coh.end(),0.0),
                                                  std::accumulate(Eplant2Coh.begin(),Eplant2Coh.end(),0.0),
                                                  std::accumulate(Eplant3Coh.begin(),Eplant3Coh.end(),0.0));
  
  //Evaporation from bare soil
  double Gsoil = soil["Gsoil"];
  double Ksoil = soil["Ksoil"];
  
  double Esoil = soilevaporation((Water_FC[0]*(1.0 - W[0])), PETsoil, Gsoil);
  
  //Exponential decay to divide bare soil evaporation among layers
  NumericVector EsoilVec = NumericVector::create(Esoil*(1.0-exp(-Ksoil*d1)),
                                                 Esoil*(exp(-Ksoil*d1)-exp(-Ksoil*(d1+d2))),
                                                 Esoil*(exp(-Ksoil*(d1+d2))));
  
  W[0] = std::max(W[0]-(EplantVec[0]+EsoilVec[0])/Water_FC[0],0.0);
  if(d2>0) W[1] = std::max(std::min(W[1]-(EplantVec[1]+EsoilVec[1])/Water_FC[1],1.0),0.0);
  if(d3>0) W[2] = std::max(std::min(W[2]-(EplantVec[2]+EsoilVec[2])/Water_FC[2],1.0),0.0);
  
  
  NumericVector DDS = (Phe*((V1*(1.0 - K1)+V2*(1.0 - K2)+V3*(1.0 - K3))));

  //For comunication with photosynthesis
  for(int c=0;c<numCohorts;c++) {
    transpiration[c] = Eplant[c];
  }
  
  List l = List::create(_["PET"] = pet, _["NetPrec"] = NetPrec, _["Runon"] = runon, _["Infiltration"] = Infiltration, _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage,
                        _["LAIcell"] = LAIcell, _["Cm"] = Cm, _["Lground"] = LgroundPAR, _["Tmax"] = Tmax,
                        _["EsoilVec"] = EsoilVec, _["EplantVec"] = EplantVec, _["psiVec"] = psi,
                        _["TmaxCoh"] = TmaxCoh, _["EplantCoh"] = Eplant, _["psiCoh"] = PlantPsi, _["DDS"] = DDS);
  return(l);
}

NumericVector getTrackSpeciesTranspiration( NumericVector trackSpecies, NumericVector Eplant, DataFrame x) {
  int nTrackSpecies = trackSpecies.size();
  NumericVector Eplantsp(nTrackSpecies, 0.0);
  NumericVector SP = x["SP"];
  int nCoh = SP.size();
  int ts;
  for(int its =0;its<nTrackSpecies;its++) {
    ts = trackSpecies[its];
    for(int i=0;i<nCoh;i++) {
      if(SP[i]==ts) {
        Eplantsp[its] += Eplant[i];
      }
    }
  }
  return(Eplantsp);
}

NumericVector getTrackSpeciesDDS(NumericVector trackSpecies, NumericVector DDS, DataFrame x) {
  int nTrackSpecies = trackSpecies.size();
  NumericVector DDSsp(nTrackSpecies, 0.0);
  NumericVector LAI = x["LAI"];
  NumericVector SP = x["SP"];
  int nCoh = LAI.size();
  int ts;
  double laiSum;
  for(int its =0;its<nTrackSpecies;its++) {
    ts = trackSpecies[its];
    laiSum = 0.0;
    for(int i=0;i<nCoh;i++) {
      if(SP[i]==ts) {
        DDSsp[its] += DDS[i]*LAI[i];
        laiSum +=LAI[i];
      }
    }
    DDSsp = DDSsp/laiSum;
  }
  return(DDSsp);
}


// [[Rcpp::export(".swbgridDay")]]
List swbgridDay(CharacterVector lct, List xList, List soilList, 
                IntegerVector waterO, List queenNeigh, List waterQ,
                NumericVector gddVec, NumericVector petVec, NumericVector rainVec, 
                NumericVector erVec, NumericVector trackSpecies) {
  int nX = xList.size();
  int nTrackSpecies = trackSpecies.size();
  NumericVector NetPrec(nX,NA_REAL), Runon(nX,0.0), Infiltration(nX,NA_REAL);
  NumericVector Runoff(nX,NA_REAL), DeepDrainage(nX,NA_REAL), W1(nX,NA_REAL), W2(nX,NA_REAL);
  NumericVector W3(nX,NA_REAL), Esoil(nX,NA_REAL), Eplant(nX,NA_REAL);
  NumericMatrix Transpiration(nX, nTrackSpecies), DDS(nX, nTrackSpecies);
  double runoffExport = 0.0;
  for(int i=0;i<nX;i++) {
    //get next cell in order
    int iCell = waterO[i]-1; //Decrease index!!!!
    if((lct[iCell]=="Wildland") || (lct[iCell]=="Agriculture") ) {
      DataFrame x = Rcpp::as<Rcpp::DataFrame>(xList[iCell]);
      List soil = Rcpp::as<Rcpp::List>(soilList[iCell]);
      //Run daily soil water balance for the current cell
      List res = swbDay1(x, soil, gddVec[iCell], petVec[iCell], rainVec[iCell], erVec[iCell], Runon[iCell]);
      NetPrec[iCell] = res["NetPrec"];
      Runon[iCell] = res["Runon"];
      Infiltration[iCell] = res["Infiltration"];
      Runoff[iCell] = res["Runoff"];
      DeepDrainage[iCell] = res["DeepDrainage"];
      Esoil[iCell] = sum(Rcpp::as<Rcpp::NumericVector>(res["EsoilVec"]));
      NumericVector EplantCoh = res["EplantCoh"];
      NumericVector DDScell = res["DDS"];
      Eplant[iCell] = sum(EplantCoh);
      if(nTrackSpecies>0) {
        Transpiration(iCell,_) = getTrackSpeciesTranspiration(trackSpecies, EplantCoh, x);
        DDS(iCell,_) = getTrackSpeciesDDS(trackSpecies, DDScell, x);
      }
      NumericVector W = soil["W"];
      W1[iCell] = W[0];
      W2[iCell] = W[1];
      W3[iCell] = W[2];

      //Assign runoff to runon of neighbours
      double ri =  Runoff[iCell];
      if(ri>0.0) {
        IntegerVector ni = Rcpp::as<Rcpp::IntegerVector>(queenNeigh[iCell]);
        NumericVector qi = Rcpp::as<Rcpp::NumericVector>(waterQ[iCell]);
        if(ni.size()>0) {
          for(int j=0;j<ni.size();j++) Runon[ni[j]-1] += qi[j]*ri; //decrease index
        } else {
          runoffExport += ri; //If no suitable neighbours add ri to landscape export via runoff
        }
      }
    } else if(lct[iCell]=="Rock") {//all Precipitation becomes runoff if cell is rock outcrop
      Runoff[iCell] =  Runon[iCell]+rainVec[iCell];
      double ri = Runoff[iCell];
      if(ri>0.0) {
        IntegerVector ni = Rcpp::as<Rcpp::IntegerVector>(queenNeigh[iCell]);
        NumericVector qi = Rcpp::as<Rcpp::NumericVector>(waterQ[iCell]);
        if(ni.size()>0) {
          for(int j=0;j<ni.size();j++) Runon[ni[j]-1] += qi[j]*ri;//decrease index
        } else {
          runoffExport += ri; //If no suitable neighbours add ri to landscape export via runoff
        }
      }
    } else if(lct[iCell]=="Static") {
      // static cells receive water from other cells or Precipitation
      // but do not export to the atmosphere contribute nor to other cells.
      // Hence, water balance over the landscape is achieved by
      // adding this water to the landscape export via runoff.
      runoffExport += Runon[iCell] + rainVec[iCell];
    }
  }
  DataFrame waterBalance = DataFrame::create(_["NetPrec"] = NetPrec, _["Runon"] = Runon, _["Infiltration"] = Infiltration,
                                   _["Runoff"] = Runoff, _["DeepDrainage"] = DeepDrainage,
                                   _["W1"] = W1, _["W2"] = W2, _["W3"] = W3,
                                   _["Esoil"] = Esoil, _["Eplant"] = Eplant);
  return(List::create(_["WaterBalance"] = waterBalance,
                      _["RunoffExport"] = runoffExport,
                      _["Transpiration"] = Transpiration,
                      _["DDS"] = DDS));
}



void checkswbInput(DataFrame x, List soil, String petMode, String hydraulicMode) {
  if(!x.containsElementNamed("LAI")) stop("LAI missing in swbInput");
  if(!x.containsElementNamed("H")) stop("H missing in swbInput");
  if(!x.containsElementNamed("V.1")) stop("V.1 missing in swbInput");
  if(!x.containsElementNamed("V.2")) stop("V.2 missing in swbInput");
  if(!x.containsElementNamed("V.3")) stop("V.3 missing in swbInput");
  if(!x.containsElementNamed("Sgdd")) stop("Sgdd missing in swbInput");
  if(petMode=="PenmanMonteith") {
    if(!x.containsElementNamed("RC_min")) stop("RC_min missing in swbInput");
  }
  if(hydraulicMode=="Simple") {
    if(!x.containsElementNamed("Psi_Extract")) stop("Psi_Extract missing in swbInput");
  } else if(hydraulicMode=="Sperry") {
    if(!x.containsElementNamed("VC_kxmax")) stop("VC_kxmax missing in swbInput");
    if(!x.containsElementNamed("VC_c")) stop("VC_c missing in swbInput");
    if(!x.containsElementNamed("VC_d")) stop("VC_d missing in swbInput");
    if(!soil.containsElementNamed("VG_ksmax")) stop("VG_ksmax missing in soil");
    if(!soil.containsElementNamed("VG_n")) stop("VG_n missing in soil");
    if(!soil.containsElementNamed("VG_alpha")) stop("VG_alpha missing in soil");
  }
  if(!x.containsElementNamed("k")) stop("k missing in swbInput");
  if(!x.containsElementNamed("g")) stop("g missing in swbInput");
  if(!x.containsElementNamed("CR")) stop("CR missing in swbInput");
  if(!soil.containsElementNamed("W")) stop("W missing in soil");
  if(!soil.containsElementNamed("psi")) stop("psi missing in soil");
  if(!soil.containsElementNamed("dVec")) stop("dVec missing in soil");
  if(!soil.containsElementNamed("Theta_FC")) stop("Theta_FC missing in soil");
  if(!soil.containsElementNamed("Water_FC")) stop("Water_FC missing in soil");
  if(!soil.containsElementNamed("macro")) stop("macro missing in soil");
  if(!soil.containsElementNamed("clay")) stop("clay missing in soil");
  if(!soil.containsElementNamed("sand")) stop("sand missing in soil");
}

// [[Rcpp::export("swb")]]
List swb(DataFrame x, List soil, DataFrame meteo, String petMode = "Input", String hydraulicMode = "Simple", double latitude = NA_REAL, double elevation = NA_REAL, bool verbose = true) {
  if((petMode!="Input") & (petMode!="PenmanMonteith")) stop("Wrong PET mode ('petMode' should be either 'Input' or 'PenmanMonteith')");
  if((hydraulicMode!="Simple") & (hydraulicMode!="Sperry")) stop("Wrong Transpiration mode ('hydraulicMode' should be either 'Simple' or 'Sperry')");
  checkswbInput(x, soil, petMode, hydraulicMode);
    
  NumericVector Precipitation = meteo["Precipitation"];
  NumericVector Temperature = meteo["MeanTemperature"];
  IntegerVector DOY = meteo["DOY"];
  NumericVector GDD = gdd(DOY, Temperature, 5.0);
  NumericVector ER = er(DOY);
  NumericVector SP = x["SP"];
  int numCohorts = x.nrows();
  int numDays = Precipitation.size();
  NumericVector PET(numDays);
  NumericVector Tmax(numDays);
  
  //Water balance variables
  NumericVector Esoil(numDays);
  NumericVector Eplant1(numDays);
  NumericVector Eplant2(numDays);
  NumericVector Eplant3(numDays);
  NumericVector LAIcell(numDays);
  NumericVector Cm(numDays);
  NumericVector Lground(numDays);
  NumericVector Runoff(numDays);
  NumericVector NetPrec(numDays);
  NumericVector Interception(numDays);
  NumericVector Infiltration(numDays);
  NumericVector DeepDrainage(numDays);
  NumericVector W1(numDays); //Soil moisture content in relation to field capacity
  NumericVector W2(numDays);
  NumericVector W3(numDays);
  NumericVector psi1(numDays);
  NumericVector psi2(numDays);
  NumericVector psi3(numDays);
  
  NumericMatrix PlantTmax(numDays, numCohorts);
  NumericMatrix PlantPsi(numDays, numCohorts);
  NumericMatrix PlantStress(numDays, numCohorts);
  NumericMatrix PlantTranspiration(numDays, numCohorts);
  NumericVector EplantCohTot(numCohorts, 0.0);
  
  
  NumericVector Water_FC = soil["Water_FC"];
  NumericVector W = soil["W"];
  W1[0] = W[0];
  W2[0] = W[1];
  W3[0] = W[2];
  
  if(verbose) Rcout << "W1i:"<< round(100*W1[0])/100<<
    " W2i:"<< round(100*W2[0])/100<<
      " W3i:"<< round(100*W3[0])/100<<"\n";
  
  if(verbose) Rcout << "Starting daily balance...\n";
  if(petMode=="Input") {
    PET = meteo["PET"];
    for(int i=0;i<numDays;i++) {
      List s = swbDay1(x, soil, GDD[i], PET[i], Precipitation[i], ER[i], 0.0, hydraulicMode, false); //No Runon in simulations for a single cell

      Tmax[i] = s["Tmax"];
      Lground[i] = s["Lground"];
      Esoil[i] = sum(Rcpp::as<Rcpp::NumericVector>(s["EsoilVec"]));
      LAIcell[i] = s["LAIcell"];
      Cm[i] = s["Cm"];
      DeepDrainage[i] = s["DeepDrainage"];
      Infiltration[i] = s["Infiltration"];
      Runoff[i] = s["Runoff"];
      NetPrec[i] = s["NetPrec"];
      Interception[i] = Precipitation[i]-NetPrec[i];
      
      NumericVector psi = s["psiVec"];
      psi1[i] = psi[0];
      psi2[i] = psi[1];
      psi3[i] = psi[2];
      
      NumericVector EplantCoh = s["EplantCoh"];
      NumericVector EplantVec = s["EplantVec"];
      Eplant1[i] = EplantVec[0];
      Eplant2[i] = EplantVec[1];
      Eplant3[i] = EplantVec[2];
      PlantTranspiration(i,_) = EplantCoh;
      PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(s["DDS"]);
      PlantTmax(i,_) = Rcpp::as<Rcpp::NumericVector>(s["TmaxCoh"]);
      PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(s["psiCoh"]);
      EplantCohTot = EplantCohTot + EplantCoh;
      
      //Stress of each plant on each layer
      if(i<(numDays-1)){
        W1[i+1] = W[0];
        W2[i+1] = W[1];
        W3[i+1] = W[2];
      }
    }
  } else if(petMode=="PenmanMonteith") {
    if(NumericVector::is_na(latitude)) stop("Value for 'latitude' should not be missing.");
    if(NumericVector::is_na(elevation)) stop("Value for 'elevation' should not be missing.");
    NumericVector MinTemperature = meteo["MinTemperature"];
    NumericVector MaxTemperature = meteo["MaxTemperature"];
    NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
    NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
    NumericVector Radiation = meteo["Radiation"];
    NumericVector WindSpeed = meteo["WindSpeed"];
    for(int i=0;i<numDays;i++) {
      List s = swbDay2(x, soil, GDD[i], MinTemperature[i], MaxTemperature[i], 
                       MinRelativeHumidity[i], MaxRelativeHumidity[i], Radiation[i], WindSpeed[i], 
                       latitude, elevation, DOY[i], Precipitation[i], ER[i], 0.0, hydraulicMode);
      
      PET[i] = s["PET"];
      Tmax[i] = s["Tmax"];
      Lground[i] = s["Lground"];
      Esoil[i] = sum(Rcpp::as<Rcpp::NumericVector>(s["EsoilVec"]));
      LAIcell[i] = s["LAIcell"];
      Cm[i] = s["Cm"];
      DeepDrainage[i] = s["DeepDrainage"];
      Infiltration[i] = s["Infiltration"];
      Runoff[i] = s["Runoff"];
      NetPrec[i] = s["NetPrec"];
      Interception[i] = Precipitation[i]-NetPrec[i];
      
      NumericVector psi = s["psiVec"];
      psi1[i] = psi[0];
      psi2[i] = psi[1];
      psi3[i] = psi[2];
      
      NumericVector EplantCoh = s["EplantCoh"];
      NumericVector EplantVec = s["EplantVec"];
      Eplant1[i] = EplantVec[0];
      Eplant2[i] = EplantVec[1];
      Eplant3[i] = EplantVec[2];
      PlantTmax(i,_) = Rcpp::as<Rcpp::NumericVector>(s["TmaxCoh"]);
      PlantTranspiration(i,_) = EplantCoh;
      PlantPsi(i,_) = Rcpp::as<Rcpp::NumericVector>(s["psiCoh"]);
      PlantStress(i,_) = Rcpp::as<Rcpp::NumericVector>(s["DDS"]);
      EplantCohTot = EplantCohTot + EplantCoh;
      
      if(i<(numDays-1)){
        W1[i+1] = W[0];
        W2[i+1] = W[1];
        W3[i+1] = W[2];
      }
    }
  }
  if(verbose) Rcout << "done\n";
  
  NumericVector ML1 = W1*Water_FC[0];
  NumericVector ML2 = W2*Water_FC[1];
  NumericVector ML3 = W3*Water_FC[2];
  NumericVector MLTot = ML1+ML2+ML3;
  
  NumericVector Eplanttot = Eplant1+Eplant2+Eplant3;
  NumericVector Etot = Eplanttot+Esoil;
  
  if(verbose) {
    double PETsum = sum(PET);
    double Precipitationsum = sum(Precipitation);
    double NetPrecsum = sum(NetPrec);
    double Interceptionsum = sum(Interception);
    double Esoilsum = sum(Esoil);
    double Runoffsum  = sum(Runoff);
    double Infiltrationsum  = sum(Infiltration);
    double DeepDrainagesum = sum(DeepDrainage);
    double Eplant1sum = sum(Eplant1);
    double Eplant2sum = sum(Eplant2);
    double Eplant3sum = sum(Eplant3);
    double Eplantsum = Eplant1sum+Eplant2sum+Eplant3sum;
    
    Rcout<<"Total PET (mm) "  <<round(PETsum) <<"\n";
    Rcout<<"Total Precipitation (mm) "  <<round(Precipitationsum) <<"\n";
    Rcout<<"Interception (mm) " << round(Interceptionsum)  <<" Net Prec (mm) " << round(NetPrecsum) <<"\n";
    Rcout<<"Infiltration (mm) " << round(Infiltrationsum)  <<
      " Runoff (mm) " << round(Runoffsum) <<
        " Deep drainage (mm) "  << round(DeepDrainagesum)  <<"\n";
    Rcout<<"Esoil (mm) " << round(Esoilsum) <<
      " Eplant (mm) "  <<round(Eplantsum) <<"\n";
    Rcout<<"W1f: "<< round(100*W1[numDays-1])/100  <<
      " W2f: "<< round(100*W2[numDays-1])/100<<
        " W3f: "<< round(100*W3[numDays-1])/100<<"\n";
    Rcout<<"Final volume: "<< round(MLTot[numDays-1])<<"\n\n";
    
  }
  if(verbose) Rcout<<"Building SWB and DWB output ...";
  
  Rcpp::DataFrame SWB = DataFrame::create(_["W1"]=W1, _["W2"]=W2, _["W3"]=W3, _["ML1"]=ML1, _["ML2"]=ML2, _["ML3"]=ML3,_["MLTot"]=MLTot,_["psi1"]=psi1,_["psi2"]=psi2, _["psi3"]=psi3);
  Rcpp::DataFrame DWB = DataFrame::create(_["LAIcell"]=LAIcell, _["Cm"]=Cm, _["Lground"] = Lground, _["PET"]=PET, _["Tmax"] = Tmax,
                                          _["Precipitation"] = Precipitation, _["NetPrec"]=NetPrec,_["Infiltration"]=Infiltration, _["Runoff"]=Runoff, _["DeepDrainage"]=DeepDrainage, 
                                          _["Etot"]=Etot,_["Esoil"]=Esoil, _["Eplanttot"]=Eplanttot,_["Eplant1"]=Eplant1, _["Eplant2"]=Eplant2,_["Eplant3"]=Eplant3);
  
  SWB.attr("row.names") = meteo.attr("row.names") ;
  DWB.attr("row.names") = meteo.attr("row.names") ;
  if(verbose) Rcout<<"plant output ...";
  
  PlantTranspiration.attr("dimnames") = List::create(meteo.attr("row.names"), SP.attr("names"));
  PlantStress.attr("dimnames") = List::create(meteo.attr("row.names"), SP.attr("names")) ;
  
  if(verbose) Rcout<<"list ...";
  List l = List::create(Named("petMode") = petMode, Named("hydraulicMode") = hydraulicMode,
                        Named("DailyBalance")=DWB, Named("SoilWaterBalance")=SWB,
                        Named("PlantTmax") = PlantTmax, Named("PlantTranspiration") = PlantTranspiration,
                        Named("PlantPsi") = PlantPsi, 
                        Named("PlantStress") = PlantStress);
  l.attr("class") = CharacterVector::create("swb","list");
  if(verbose) Rcout<<"done.\n";
  return(l);
}
