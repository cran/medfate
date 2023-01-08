#define STRICT_R_HEADERS
#include <Rcpp.h>
#include <numeric>
#include "lightextinction.h"
#include "windextinction.h"
#include "windKatul.h"
#include "hydraulics.h"
#include "biophysicsutils.h"
#include "phenology.h"
#include "forestutils.h"
#include "tissuemoisture.h"
#include "carbon.h"
#include "photosynthesis.h"
#include "root.h"
#include "soil.h"
#include <meteoland.h>
using namespace Rcpp;

const double SIGMA_Wm2 = 5.67*1e-8;
const double Cp_JKG = 1013.86; // J * kg^-1 * ºC^-1
const double Cp_Jmol = 29.37152; // J * mol^-1 * ºC^-1
const double eps_xylem = 1e3; // xylem elastic modulus (1 GPa = 1000 MPa)


List profitMaximization2(List supplyFunction, int initialPos,
                         double Catm, double Patm, double Tair, double vpa, double u, 
                         double SWRabs, double LWRnet, double Q, double Vmax298, double Jmax298, 
                         double leafWidth, double refLeafArea,
                         double Gswmin, double Gswmax) {
  
  NumericVector E = supplyFunction["E"];
  NumericVector dEdP = supplyFunction["dEdP"];
  NumericVector leafPsi = supplyFunction["psiLeaf"];
  int nsteps = E.size();
  
  double maxdEdp = 0.0, mindEdp = 99999999.0;

  int valid = 1;
  for(int i=0;i<nsteps-1;i++) {
    if((i>0) && (E[i] > E[i-1])) valid++; 
    mindEdp = std::min(mindEdp, dEdP[i]);
    maxdEdp = std::max(maxdEdp, dEdP[i]);
  }
  nsteps = valid;
  // Rcout<< valid<< " "<< maxdEdp << " "<<mindEdp<< " "<< mindEdp/maxdEdp<<"\n";
  // initial pos cannot be over the valid steps
  initialPos = std::min(initialPos, nsteps-1);
  
  //Evaluate photosynthesis and profit at Agmax
  NumericVector photoAgMax = leafPhotosynthesisOneFunction2(E[nsteps-1], leafPsi[nsteps - 1], Catm, Patm, Tair, vpa, u, 
                                                          SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                                          leafWidth, refLeafArea, false);
  double Agmax = photoAgMax["GrossPhotosynthesis"];
  double profitAgMax = (1.0 - ((maxdEdp - dEdP[nsteps-1])/(maxdEdp - mindEdp))); 
  

  //Photosynthesis and profit maximization at current value of initialPos
  NumericVector photoInitial = leafPhotosynthesisOneFunction2(E[initialPos], leafPsi[initialPos], Catm, Patm, Tair, vpa, u, 
                                                           SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                                           leafWidth, refLeafArea, false);
  double AgInitial = photoInitial["GrossPhotosynthesis"];
  double profitInitial = (AgInitial/Agmax) - ((maxdEdp - dEdP[initialPos])/(maxdEdp - mindEdp)); 
  if(Agmax ==0.0) profitInitial = - ((maxdEdp - dEdP[initialPos])/(maxdEdp - mindEdp)); //To avoid 0/0 division
    
  NumericVector photoPrev, photoNext, photoCurrent;
  int iPos = initialPos;
  double profitPrev, profitNext, profitCurrent;
  
  photoCurrent = photoInitial;
  profitCurrent = profitInitial;
  
  int prevStep = 0;
  bool cont = true;
  // Rcout<< " initialPos " << initialPos;
  while(cont) {
    // Rcout<<".";
    if((iPos>0) && (prevStep <= 0)) {
      photoPrev = leafPhotosynthesisOneFunction2(E[iPos-1], leafPsi[iPos-1], Catm, Patm, Tair, vpa, u, 
                                                 SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                                 leafWidth, refLeafArea, false);
      double AgPrev = photoPrev["GrossPhotosynthesis"];
      profitPrev = (AgPrev/Agmax) - ((maxdEdp - dEdP[iPos-1])/(maxdEdp - mindEdp)); 
      if(Agmax ==0.0) profitPrev = - ((maxdEdp - dEdP[iPos-1])/(maxdEdp - mindEdp)); 
    } else {
      photoPrev = photoCurrent;
      profitPrev = profitCurrent;
    }
    if((iPos < (nsteps-1)) && (prevStep >= 0)) {
      photoNext = leafPhotosynthesisOneFunction2(E[iPos+1], leafPsi[iPos+1], Catm, Patm, Tair, vpa, u, 
                                                 SWRabs, LWRnet, Q, Vmax298, Jmax298, 
                                                 leafWidth, refLeafArea, false);
      double AgNext = photoNext["GrossPhotosynthesis"];
      profitNext = (AgNext/Agmax) - ((maxdEdp - dEdP[iPos+1])/(maxdEdp - mindEdp)); 
      if(Agmax ==0.0) profitNext = - ((maxdEdp - dEdP[iPos+1])/(maxdEdp - mindEdp)); 
    } else {
      photoNext = photoAgMax;
      profitNext = profitAgMax;
    }
    bool selDecr = ((profitPrev >= profitCurrent) && (photoPrev["Gsw"] >= Gswmin)) || (photoCurrent["Gsw"]>Gswmax);
    selDecr = selDecr && (prevStep<=0) && (iPos>0);
    bool selIncr = ((profitNext > profitCurrent) && (photoNext["Gsw"] <= Gswmax)) || (photoCurrent["Gsw"]<Gswmin);
    selIncr = selIncr && (prevStep>=0) && (iPos<(nsteps-1));
    if(selDecr) {
      profitCurrent = profitPrev;
      photoCurrent = photoPrev;
      iPos = iPos-1;
      prevStep = -1;
      // Rcout <<"-";
    } else if(selIncr) {
      profitCurrent = profitNext;
      photoCurrent = photoNext;
      iPos = iPos+1;
      // Rcout <<"+";
      prevStep = 1;
    } else {
      cont = false;
    }
  }
  // Rcout <<"\n";
  // Rcout<< " finalPos " << iPos<< " final profit "<< profitCurrent <<"\n"; 
  
  return(List::create(Named("photosynthesisFunction") = photoCurrent,
                      Named("Profit") = profitCurrent,
                      Named("iMaxProfit")=iPos));
}

//' Stomatal regulation
//' 
//' Set of high-level functions used in the calculation of stomatal conductance and transpiration. 
//' Function \code{transp_profitMaximization} calculates gain and cost functions, 
//' as well as profit maximization from supply and photosynthesis input functions. 
//' Function \code{transp_stomatalRegulationPlot} produces a plot with the cohort supply functions against water potential 
//' and a plot with the cohort photosynthesis functions against water potential, both with the maximum profit values indicated.
//' 
//' @param supplyFunction Water supply function (see \code{\link{hydraulics_supplyFunctionNetwork}}).
//' @param photosynthesisFunction Function returned by \code{photo_photosynthesisFunction()}.
//' @param Gswmin,Gswmax Minimum and maximum stomatal conductance to water vapour (mol·m-2·s-1).
//' 
//' @return
//' Function \code{transp_profitMaximization} returns a list with the following elements:
//' \itemize{
//'   \item{\code{Cost}: Cost function [0-1].}
//'   \item{\code{Gain}: Gain function [0-1].}
//'   \item{\code{Profit}: Profit function [0-1].}
//'   \item{\code{iMaxProfit}: Index corresponding to maximum profit (starting from 0).}
//' }
//' 
//' @references
//' Sperry, J. S., M. D. Venturas, W. R. L. Anderegg, M. Mencuccini, D. S. Mackay, Y. Wang, and D. M. Love. 2017. Predicting stomatal responses to the environment from the optimization of photosynthetic gain and hydraulic cost. Plant Cell and Environment 40, 816-830 (doi: 10.1111/pce.12852).
//' 
//' @author Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
//' 
//' @seealso
//' \code{\link{transp_transpirationSperry}}, \code{\link{hydraulics_supplyFunctionNetwork}}, \code{\link{biophysics_leafTemperature}}, \code{\link{photo_photosynthesis}}, \code{\link{spwb_day}}, \code{\link{plot.spwb_day}}
//' 
//' @examples
//' #Load example daily meteorological data
//' data(examplemeteo)
//' 
//' #Load example plot plant data
//' data(exampleforestMED)
//' 
//' #Default species parameterization
//' data(SpParamsMED)
//' 
//' #Initialize soil with default soil params (4 layers)
//' examplesoil = soil(defaultSoilParams(4))
//' 
//' #Initialize control parameters
//' control = defaultControl(transpirationMode="Sperry")
//' 
//' #Initialize input
//' x2 = forest2spwbInput(exampleforestMED,examplesoil, SpParamsMED, control)
//' 
//' # Stomatal VPD curve and chosen value for the 12th time step at day 100
//' transp_stomatalRegulationPlot(x2, examplemeteo, day=100, timestep = 12,
//'                               latitude = 41.82592, elevation = 100, type="VPD")
//'  
//' @name transp_stomatalregulation
// [[Rcpp::export("transp_profitMaximization")]]
List profitMaximization(List supplyFunction, DataFrame photosynthesisFunction, double Gswmin, double Gswmax) {
  NumericVector supplyE = supplyFunction["E"];
  NumericVector supplydEdp = supplyFunction["dEdP"];
  NumericVector Ag = photosynthesisFunction["GrossPhotosynthesis"];
  NumericVector leafTemp = photosynthesisFunction["LeafTemperature"];
  NumericVector leafVPD = photosynthesisFunction["LeafVPD"];
  NumericVector Gsw = photosynthesisFunction["Gsw"];
  NumericVector supplyKterm = supplyFunction["kterm"];
  int nsteps = supplydEdp.size();
  double maxdEdp = 0.0, mindEdp = 99999999.0;
  double Agmax = 0.0;
  //Find valid limits according to stomatal conductance
  int ini = 0, fin = nsteps-1;
 
  for(int i=ini;i<fin;i++) {
    mindEdp = std::min(mindEdp, supplydEdp[i]);
    maxdEdp = std::max(maxdEdp, supplydEdp[i]);
    Agmax = std::max(Agmax, Ag[i]);
  }
  
  //Evaluate profit for valid steps
  NumericVector profit(nsteps, NA_REAL);
  NumericVector cost(nsteps, NA_REAL);
  NumericVector gain(nsteps, NA_REAL);
  for(int i=ini;i<fin;i++) {
    gain[i] = Ag[i]/Agmax;
    cost[i] = (maxdEdp-supplydEdp[i])/(maxdEdp-mindEdp); 
    profit[i] = gain[i]-cost[i];
  }
  
  while((Gsw[ini]<=Gswmin) && (ini<fin)) ini++;
  while((Gsw[fin]>=Gswmax) && (fin>ini)) fin--; 
  
  //Ensure that ini <=fin
  ini = std::min(ini, fin);
  fin = std::max(ini,fin);
  
  int imaxprofit=ini;
  double maxprofit=profit[ini];
  if(fin>ini) {
    for(int i=ini+1;i<=fin;i++){
      if((profit[i]>maxprofit)) {
        maxprofit = profit[i];
        imaxprofit = i;
      }
    }
  }
  // Rcout<<ini<< " "<< fin<< " Gsw= " << Gsw[imaxprofit] <<" Gswmax= "<<Gswmax<<" Gswmin "<<Gswmin<<" iPM="<< imaxprofit<<" Eini=" <<supplyE[ini]<<" Efin=" <<supplyE[fin]<<" E[iPM]=" <<supplyE[imaxprofit]<<"\n";
  if((Gsw[imaxprofit] > Gswmax) && (imaxprofit>ini)) {
    Rcout<<ini<< " "<< fin<< " Gsw= " << Gsw[imaxprofit] <<" Gswmax= "<<Gswmax<<" Gswmin "<<Gswmin<<" iPM="<< imaxprofit<<" Eini=" <<supplyE[ini]<<" Efin=" <<supplyE[fin]<<" E[iPM]=" <<supplyE[imaxprofit]<<"\n";
    for(int i=0;i<Gsw.size();i++) {
      Rcout<< i << " Gsw "<< Gsw[i] << " supplyE "<< supplyE[i] << " leafT "<< leafTemp[i]<< " leafVPD "<< leafVPD[i]  << "\n";
    }
    stop("Gsw > Gswmax");
  }
  return(List::create(Named("Cost") = cost,
                      Named("Gain") = gain,
                      Named("Profit") = profit,
                      Named("iMaxProfit")=imaxprofit));
}


void copyRhizoPsi(int c, int iPM, 
                  NumericMatrix RhizoPsi, NumericMatrix RhizoPsiMAT,
                  LogicalMatrix layerConnected, 
                  List RHOP, List layerConnectedPools,
                  NumericVector VCroot_c, NumericVector VCroot_d,  
                  bool plantWaterPools) {
  int nlayers = layerConnected.ncol();
  int numCohorts = layerConnected.nrow();
  
  if(!plantWaterPools) {
    int cl = 0;
    for(int l=0;l<nlayers;l++) {
      if(layerConnected(c,l)) {
        RhizoPsiMAT(c,l) = RhizoPsi(iPM,cl);
        cl++;
      } 
    }
  } else {
    NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
    LogicalMatrix layerConnectedCoh = Rcpp::as<Rcpp::LogicalMatrix>(layerConnectedPools[c]);
    NumericVector rplv(numCohorts,NA_REAL);
    NumericVector vplv(numCohorts,NA_REAL);
    int cl = 0;
    for(int l=0;l<nlayers;l++) {
      int clj = 0;
      for(int j=0;j<numCohorts;j++) {
        if(layerConnectedCoh(j,l)) {
          rplv[clj] = RhizoPsi(iPM,cl);
          vplv[clj] = RHOPcoh(c,l);
          cl++;
          clj++;
        }
      }
      NumericVector pv(clj,NA_REAL);
      NumericVector vv(clj,NA_REAL);
      for(int j=0;j<clj;j++) {
        pv[j] = rplv[j];
        vv[j] = vplv[j];
      }
      RhizoPsiMAT(c,l) = averagePsi(pv,vv,VCroot_c[c], VCroot_d[c]);
    }
  }
}

List transpirationSperry(List x, NumericVector meteovec, 
                  double latitude, double elevation, double slope, double aspect, 
                  double solarConstant, double delta,
                  double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                  bool verbose = false, int stepFunctions = NA_INTEGER, 
                  bool modifyInput = true) {
  //Control parameters
  List control = x["control"];
  String soilFunctions = control["soilFunctions"];
  List numericParams = control["numericParams"];
  int ntrial = numericParams["ntrial"];
  int maxNsteps  = numericParams["maxNsteps"];
  double psiTol = numericParams["psiTol"];
  double ETol = numericParams["ETol"];
  bool capacitance = control["capacitance"];
  String cavitationRefill = control["cavitationRefill"];
  double refillMaximumRate = control["refillMaximumRate"];
  double klatleaf = control["klatleaf"];
  double klatstem = control["klatstem"];
  int ntimesteps = control["ndailysteps"];
  int nsubsteps = control["nsubsteps"];
  String rhizosphereOverlap = control["rhizosphereOverlap"];
  bool plantWaterPools = (rhizosphereOverlap!="total");
  double verticalLayerSize = control["verticalLayerSize"];
  double windMeasurementHeight  = control["windMeasurementHeight"];
  double thermalCapacityLAI = control["thermalCapacityLAI"];
  bool multiLayerBalance = control["multiLayerBalance"];
  double defaultWindSpeed = control["defaultWindSpeed"];
  double nonSugarConcentration = control["nonSugarConcentration"];
  
  //Meteo input
  double tmin = meteovec["tmin"];
  double tmax = meteovec["tmax"];
  double tminPrev = meteovec["tminPrev"];
  double tmaxPrev = meteovec["tmaxPrev"];
  double tminNext = meteovec["tminNext"];
  double prec = meteovec["prec"];
  double rhmin = meteovec["rhmin"];
  double rhmax = meteovec["rhmax"];
  double rad = meteovec["rad"];
  double wind = meteovec["wind"];
  double Catm = meteovec["Catm"];

  //Vegetation input
  DataFrame cohorts = Rcpp::as<Rcpp::DataFrame>(x["cohorts"]);
  DataFrame above = Rcpp::as<Rcpp::DataFrame>(x["above"]);
  NumericVector LAIlive = Rcpp::as<Rcpp::NumericVector>(above["LAI_live"]);
  NumericVector LAIphe = Rcpp::as<Rcpp::NumericVector>(above["LAI_expanded"]);
  NumericVector LAIdead = Rcpp::as<Rcpp::NumericVector>(above["LAI_dead"]);
  NumericVector H = Rcpp::as<Rcpp::NumericVector>(above["H"]);
  NumericVector CR = Rcpp::as<Rcpp::NumericVector>(above["CR"]);
  NumericVector N = Rcpp::as<Rcpp::NumericVector>(above["N"]);
  
  int numCohorts = LAIlive.size();
  
  //Soil input
  List soil = x["soil"];
  NumericVector dVec = soil["dVec"];
  int nlayers = dVec.length();
  NumericVector Water_FC = waterFC(soil, soilFunctions);
  NumericVector Theta_FC = thetaFC(soil, soilFunctions);
  NumericVector VG_n = Rcpp::as<Rcpp::NumericVector>(soil["VG_n"]);
  NumericVector VG_alpha = Rcpp::as<Rcpp::NumericVector>(soil["VG_alpha"]);
  NumericVector Tsoil = soil["Temp"]; 
  NumericVector sand = soil["sand"];
  NumericVector clay = soil["clay"];
  NumericVector Ws = soil["W"]; //Access to soil state variable
  double SWE = soil["SWE"];
  NumericVector psiSoil = psi(soil, soilFunctions); //Get soil water potential
  
  //Canopy params
  DataFrame canopyParams = Rcpp::as<Rcpp::DataFrame>(x["canopy"]);
  NumericVector zlow = canopyParams["zlow"];
  NumericVector zmid = canopyParams["zmid"];
  NumericVector zup = canopyParams["zup"];
  NumericVector Tair = canopyParams["Tair"];
  NumericVector VPair = canopyParams["VPair"];
  NumericVector Cair = canopyParams["Cair"];
  int ncanlayers = Tair.size(); //Number of canopy layers
  for(int l=0;l<ncanlayers;l++) { //If canopy layers have missing values, then initialize with Catm
    if(!multiLayerBalance) Cair[l] = Catm;
    else {
      if(NumericVector::is_na(Cair[l])) Cair[l] = Catm; 
    }
  }
  if(multiLayerBalance) Cair[ncanlayers-1] = Catm;
  
  //Root distribution input
  DataFrame belowdf = Rcpp::as<Rcpp::DataFrame>(x["below"]);
  List belowLayers = Rcpp::as<Rcpp::List>(x["belowLayers"]);
  NumericMatrix V =Rcpp::as<Rcpp::NumericMatrix>(belowLayers["V"]);
  NumericMatrix VCroot_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VCroot_kmax"]);
  NumericMatrix VGrhizo_kmax= Rcpp::as<Rcpp::NumericMatrix>(belowLayers["VGrhizo_kmax"]);
  
  //Water pools
  NumericMatrix Wpool = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["Wpool"]);
  List RHOP;
  NumericVector poolProportions(numCohorts);
  if(plantWaterPools) {
    RHOP = belowLayers["RHOP"];
    poolProportions = belowdf["poolProportions"];
  }
  
  //Base parameters
  DataFrame paramsInterception = Rcpp::as<Rcpp::DataFrame>(x["paramsInterception"]);
  NumericVector alphaSWR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["alphaSWR"]);
  NumericVector gammaSWR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["gammaSWR"]);
  NumericVector kPAR = Rcpp::as<Rcpp::NumericVector>(paramsInterception["kPAR"]);
  
  //Anatomy parameters
  DataFrame paramsAnatomy = Rcpp::as<Rcpp::DataFrame>(x["paramsAnatomy"]);
  NumericVector leafWidth = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["LeafWidth"]);
  NumericVector Al2As = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["Al2As"]);
  NumericVector r635 = Rcpp::as<Rcpp::NumericVector>(paramsAnatomy["r635"]);
  
  //Transpiration parameters
  DataFrame paramsTransp = Rcpp::as<Rcpp::DataFrame>(x["paramsTranspiration"]);
  NumericVector Gswmin = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmin"]);
  NumericVector Gswmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Gswmax"]);
  NumericVector Plant_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["Plant_kmax"]);
  NumericVector VCstem_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_kmax"]);
  NumericVector VCstem_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_c"]);
  NumericVector VCstem_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCstem_d"]);
  NumericVector VCleaf_kmax = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_kmax"]);
  NumericVector VCleaf_c = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_c"]);
  NumericVector VCleaf_d = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCleaf_d"]);
  NumericVector VCroot_kmax_sum = Rcpp::as<Rcpp::NumericVector>(paramsTransp["VCroot_kmax"]);
  NumericVector VCroot_c = paramsTransp["VCroot_c"];
  NumericVector VCroot_d = paramsTransp["VCroot_d"];
  NumericVector Vmax298 = paramsTransp["Vmax298"];
  NumericVector Jmax298 = paramsTransp["Jmax298"];

  //Water storage parameters
  DataFrame paramsWaterStorage = Rcpp::as<Rcpp::DataFrame>(x["paramsWaterStorage"]);
  NumericVector maxFMC = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["maxFMC"]);
  NumericVector StemPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemPI0"]);
  NumericVector StemEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemEPS"]);
  NumericVector StemAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["StemAF"]);
  NumericVector Vsapwood = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vsapwood"]); //l·m-2 = mm
  NumericVector LeafPI0 = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafPI0"]);
  NumericVector LeafEPS = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafEPS"]);
  NumericVector LeafAF = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["LeafAF"]);
  NumericVector Vleaf = Rcpp::as<Rcpp::NumericVector>(paramsWaterStorage["Vleaf"]); //l·m-2 = mm
  

  //Comunication with outside
  DataFrame internalWater = Rcpp::as<Rcpp::DataFrame>(x["internalWater"]);
  NumericVector StemPLCVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemPLC"]);
  NumericVector StemSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["StemSympPsi"]);
  NumericVector LeafSympPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafSympPsi"]);
  NumericVector LeafPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["LeafPsi"]);
  NumericVector Stem1PsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem1Psi"]);
  NumericVector Stem2PsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Stem2Psi"]);
  NumericVector RootCrownPsiVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["RootCrownPsi"]);
  NumericMatrix RhizoPsiMAT = Rcpp::as<Rcpp::NumericMatrix>(belowLayers["RhizoPsi"]);
  NumericVector EinstVEC = Rcpp::as<Rcpp::NumericVector>(internalWater["Einst"]);

  if(NumericVector::is_na(aspect)) aspect = 0.0;
  if(NumericVector::is_na(slope)) slope = 0.0;
  double latrad = latitude * (M_PI/180.0);
  double asprad = aspect * (M_PI/180.0);
  double slorad = slope * (M_PI/180.0);
  
  //Step in seconds
  double tstep = 86400.0/((double) ntimesteps);

  //Atmospheric pressure
  double Patm = meteoland::utils_atmosphericPressure(elevation);
  
 
  //Daily average water vapor pressure at the atmosphere (kPa)
  double vpatm = meteoland::utils_averageDailyVP(tmin, tmax, rhmin,rhmax);
  //If canopy VP is missing or not multilayer initiate it to vpatm
  if(NumericVector::is_na(VPair[0]) || (!multiLayerBalance)){
    for(int i=0;i<ncanlayers;i++) VPair[i] = vpatm;
  }
  //Daily cloud cover
  double cloudcover = 0.0;
  if(prec >0.0) cloudcover = 1.0;
  bool clearday = (prec==0);
  
  
  //1. Leaf Phenology: Adjusted leaf area index
  NumericVector Phe(numCohorts);
  double LAIcell = 0.0, LAIcelldead = 0.0, LAIcelllive = 0.0, LAIcellexpanded = 0.0;
  double canopyHeight = 100.0; //Minimum canopy height of 1 m
  for(int c=0;c<numCohorts;c++) {
    Phe[c]=LAIphe[c]/LAIlive[c]; //Phenological status
    if(LAIlive[c]==0.0) Phe[c] = 0.0;
    LAIcell += (LAIphe[c]+LAIdead[c]);
    LAIcelldead += LAIdead[c];
    LAIcelllive += LAIlive[c];
    LAIcellexpanded +=LAIphe[c];
    if((canopyHeight<H[c]) && ((LAIphe[c]+LAIdead[c])>0.0)) canopyHeight = H[c];
  }
  //Create z vector with all layer height limits
  NumericVector z(ncanlayers+1,0.0);
  for(int i=1;i<=ncanlayers;i++) z[i] = z[i-1] + verticalLayerSize;
  //LAI distribution per layer and cohort
  NumericMatrix LAIme = LAIdistributionVectors(z, LAIphe, H, CR); //Expanded leaves
  NumericMatrix LAImd = LAIdistributionVectors(z, LAIdead, H, CR); //Dead (standing) leaves
  NumericMatrix LAImx = LAIdistributionVectors(z, LAIlive, H, CR); //Maximum leaf expansion
  //LAI profile per layer
  NumericVector LAIpx = LAIprofileVectors(z, LAIlive, H, CR);
  NumericVector LAIpe = LAIprofileVectors(z, LAIphe, H, CR);
  NumericVector LAIpd = LAIprofileVectors(z, LAIdead, H, CR);
  NumericVector lad = 100.0*(LAIpe + LAIpd)/verticalLayerSize;
  //3. Wind extinction profile
  if(NumericVector::is_na(wind)) wind = defaultWindSpeed; //set to default if missing
  wind = std::min(10.0, std::max(wind, 0.1)); //Bound between 0.1 m/s (0.36 km/h)  and 10 m/s (36 km/h)
  DataFrame canopyTurbulence = NA_REAL;
  NumericVector zWind(ncanlayers,wind), dU(ncanlayers, 0.0), uw(ncanlayers, 0.0);
  if(canopyHeight>0.0) {
    canopyTurbulence = windCanopyTurbulence(zmid, lad,  canopyHeight, 
                                                      wind, windMeasurementHeight);
    zWind = canopyTurbulence["u"]; 
    dU = Rcpp::as<Rcpp::NumericVector>(canopyTurbulence["du"]);
    uw = canopyTurbulence["uw"];
  } 
  //4a. Instantaneous direct and diffuse shorwave radiation
  DataFrame ddd = meteoland::radiation_directDiffuseDay(solarConstant, latrad, slorad, asprad, delta,
                                                        rad, clearday, ntimesteps);
  NumericVector solarHour = ddd["SolarHour"]; //in radians
  
  //4b. Instantaneous air temperature (above canopy) and longwave radiation
  NumericVector Tatm(ntimesteps), lwdr(ntimesteps), Tcan(ntimesteps, NA_REAL), Tsunrise(ntimesteps);
  NumericVector net_LWR_can(ntimesteps),LEcan_heat(ntimesteps), Hcan_heat(ntimesteps), Ebal(ntimesteps);
  NumericVector net_LWR_soil(ntimesteps), Ebalsoil(ntimesteps), Hcansoil(ntimesteps), LEsoil_heat(ntimesteps);
  NumericMatrix Tsoil_mat(ntimesteps, nlayers);
  NumericMatrix Tcan_mat(ntimesteps, ncanlayers);
  NumericMatrix VPcan_mat(ntimesteps, ncanlayers);
  //Daylength in seconds (assuming flat area because we want to model air temperature variation)
  double tauday = meteoland::radiation_daylengthseconds(latrad,0.0,0.0, delta); 
  for(int n=0;n<ntimesteps;n++) {
    //From solar hour (radians) to seconds from sunrise
    Tsunrise[n] = (solarHour[n]*43200.0/M_PI)+ (tauday/2.0) +(tstep/2.0); 
    //Calculate instantaneous temperature and light conditions
    Tatm[n] = temperatureDiurnalPattern(Tsunrise[n], tmin, tmax, tminPrev, tmaxPrev, tminNext, tauday);
    //Longwave sky diffuse radiation (W/m2)
    lwdr[n] = meteoland::radiation_skyLongwaveRadiation(Tatm[n], vpatm, cloudcover);
  }
  if(NumericVector::is_na(Tair[0])) {//If missing initialize canopy profile with atmospheric air temperature 
    for(int i=0;i<ncanlayers;i++) Tair[i] = Tatm[0];
  }
  if(NumericVector::is_na(Tsoil[0])) {//If missing initialize soil temperature with atmospheric air temperature 
    for(int l=0;l<nlayers; l++) Tsoil[l] = Tatm[0];
  }
  //Take initial canopy air temperature from previous day
  Tcan[0] = sum(Tair*LAIpx)/sum(LAIpx);
  for(int j=0;j<ncanlayers; j++) {
    Tcan_mat(0,j) = Tair[j];
    VPcan_mat(0,j) = VPair[j];
  }
  //Take temperature soil vector 
  Tsoil_mat(0,_) = Tsoil; 
  
  
  
  //4c. Light extinction and absortion by time steps
  List lightExtinctionAbsortion = instantaneousLightExtinctionAbsortion(LAIme, LAImd, LAImx,
                                                                        kPAR, alphaSWR, gammaSWR,
                                                                        ddd, 
                                                                        ntimesteps, 0.1);
  List sunshade = lightExtinctionAbsortion["sunshade"];
  List abs_PAR_SL_COH_list = sunshade["PAR_SL"];
  List abs_PAR_SH_COH_list = sunshade["PAR_SH"];
  List abs_SWR_SL_COH_list = sunshade["SWR_SL"];
  List abs_SWR_SH_COH_list = sunshade["SWR_SH"];
  List multilayer = lightExtinctionAbsortion["multilayer"];
  List abs_SWR_SL_ML_list = multilayer["SWR_SL"];
  List abs_SWR_SH_ML_list = multilayer["SWR_SH"];
  NumericVector fsunlit = lightExtinctionAbsortion["fsunlit"];
  NumericVector abs_SWR_can = lightExtinctionAbsortion["SWR_can"];
  NumericVector abs_SWR_soil = lightExtinctionAbsortion["SWR_soil"];

  NumericVector LAI_SL(numCohorts,0.0);
  NumericVector LAI_SH(numCohorts,0.0);
  NumericVector Vmax298SL(numCohorts,0.0);
  NumericVector Vmax298SH(numCohorts,0.0);
  NumericVector Jmax298SL(numCohorts,0.0);
  NumericVector Jmax298SH(numCohorts,0.0);

  for(int c=0;c<numCohorts;c++) {
    // Rcout<<"cohort "<<c<<":\n";
    //Constant properties through time steps
    NumericVector Vmax298layer(ncanlayers), Jmax298layer(ncanlayers);
    NumericVector SLarealayer(ncanlayers), SHarealayer(ncanlayers);
    double sn =0.0;
    for(int i=(ncanlayers-1);i>=0.0;i--) {
      //Effect of nitrogen concentration decay through the canopy
      double fn = exp(-0.713*(sn+LAIme(i,c)/2.0)/sum(LAIme(_,c)));
      // Rcout<<" l"<<i<<" fsunlit: "<< fsunlit[i]<<" lai: "<< LAIme(i,c)<<" fn: "<< fn <<"\n";
      sn+=LAIme(i,c);
      SLarealayer[i] = LAIme(i,c)*fsunlit[i];
      SHarealayer[i] = LAIme(i,c)*(1.0-fsunlit[i]);
      Vmax298layer[i] = Vmax298[c]*fn;
      Jmax298layer[i] = Jmax298[c]*fn;
    }
    for(int i=0;i<ncanlayers;i++) {
      LAI_SL[c] +=SLarealayer[i];
      LAI_SH[c] +=SHarealayer[i];
      Vmax298SL[c] +=Vmax298layer[i]*LAIme(i,c)*fsunlit[i];
      Jmax298SL[c] +=Jmax298layer[i]*LAIme(i,c)*fsunlit[i];
      Vmax298SH[c] +=Vmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
      Jmax298SH[c] +=Jmax298layer[i]*LAIme(i,c)*(1.0-fsunlit[i]);
    }
    
  }
  // double kb = lightExtinctionAbsortion["kb"];  //Proportion of sunlit extinction coefficient
  // double gbf = lightExtinctionAbsortion["gbf"]; //Ground fractions
  // double gdf = lightExtinctionAbsortion["gdf"];
  
  //Determine canopy vertical layer corresponding to cohort canopy, sunlit and shade leaves for each cohort
  IntegerVector iLayerCohort(numCohorts), iLayerSunlit(numCohorts), iLayerShade(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    double num = 0.0, den = 0.0, numsl=0.0, densl =0.0, numsh = 0.0, densh=0.0;
    for(int i=0;i<ncanlayers;i++) {
      num += LAIme(i,c)*zmid[i];
      den += LAIme(i,c);
      numsl += LAIme(i,c)*zmid[i]*fsunlit[i];
      densl += LAIme(i,c)*fsunlit[i];
      numsh += LAIme(i,c)*zmid[i]*(1.0 - fsunlit[i]);
      densh += LAIme(i,c)*(1.0-fsunlit[i]);
    }
    double hc_sl = numsl/densl;
    double hc_sh = numsh/densh;
    double hc  = num/den;
    for(int i=0;i<ncanlayers;i++) {
      if((hc > zlow[i]) && (hc <=zup[i])) iLayerCohort[c] = i;
      if((hc_sl > zlow[i]) && (hc_sl <=zup[i])) iLayerSunlit[c] = i;
      if((hc_sh > zlow[i]) && (hc_sh <=zup[i])) iLayerShade[c] = i;
    }
    // Rcout << c << " "<< hc_sl<<" "<< iLayerSunlit[c]<< " "<< hc_sh<<" "<< iLayerShade[c]<<"\n";
  }
  
  //Hydraulics: determine layers where the plant is connected
  IntegerVector nlayerscon(numCohorts,0);
  LogicalMatrix layerConnected(numCohorts, nlayers);
  List layerConnectedPools(numCohorts);
  
  NumericMatrix SoilWaterExtract(numCohorts, nlayers);
  std::fill(SoilWaterExtract.begin(), SoilWaterExtract.end(), 0.0);
  NumericMatrix soilLayerExtractInst(nlayers, ntimesteps);
  std::fill(soilLayerExtractInst.begin(), soilLayerExtractInst.end(), 0.0);


  //Average sap fluidity
  double sapFluidityDay = 1.0/waterDynamicViscosity((tmin+tmax)/2.0);
  
  //Hydraulics: supply functions
  List hydraulicNetwork(numCohorts);
  List supply(numCohorts);
  List supplyAboveground(numCohorts);
  supply.attr("names") = above.attr("row.names");
  for(int c=0;c<numCohorts;c++) {
    
    if(!plantWaterPools) {
      //Determine connected layers (non-zero fine root abundance)
      nlayerscon[c] = 0;
      for(int l=0;l<nlayers;l++) {
        if(V(c,l)>0.0) {
          layerConnected(c,l)= true;
          nlayerscon[c]=nlayerscon[c]+1;
        } else {
          layerConnected(c,l) = false;
        }
      }
      // Rcout<<c<<" "<< nlayerscon[c]<<"\n";
      if(nlayerscon[c]==0) stop("Plant cohort not connected to any soil layer!");
      
      // Copy values from connected layers
      NumericVector Vc = NumericVector(nlayerscon[c]);
      NumericVector VCroot_kmaxc = NumericVector(nlayerscon[c]);
      NumericVector VGrhizo_kmaxc = NumericVector(nlayerscon[c]);
      NumericVector psic = NumericVector(nlayerscon[c]);
      NumericVector VG_nc = NumericVector(nlayerscon[c]);
      NumericVector VG_alphac= NumericVector(nlayerscon[c]);
      int cnt=0;
      for(int l=0;l<nlayers;l++) {
        if(layerConnected(c,l)) {
          Vc[cnt] = V(c,l);
          VCroot_kmaxc[cnt] = VCroot_kmax(c,l);
          VGrhizo_kmaxc[cnt] = VGrhizo_kmax(c,l);
          psic[cnt] = psiSoil[l];
          VG_nc[cnt] = VG_n[l];
          VG_alphac[cnt] = VG_alpha[l];
          cnt++;
        }
      }
      
      
      //Build supply function networks 
      if(!capacitance) {
        hydraulicNetwork[c] = List::create(_["psisoil"] = psic,
                                 _["krhizomax"] = VGrhizo_kmaxc,_["nsoil"] = VG_nc,_["alphasoil"] = VG_alphac,
                                 _["krootmax"] = sapFluidityDay*VCroot_kmaxc, _["rootc"] = VCroot_c[c], _["rootd"] = VCroot_d[c],
                                 _["kstemmax"] = sapFluidityDay*VCstem_kmax[c], _["stemc"] = VCstem_c[c], _["stemd"] = VCstem_d[c],
                                 _["kleafmax"] = sapFluidityDay*VCleaf_kmax[c], _["leafc"] = VCleaf_c[c], _["leafd"] = VCleaf_d[c],
                                 _["PLCstem"] = NumericVector::create(StemPLCVEC[c],StemPLCVEC[c]));
        supply[c] = supplyFunctionNetwork(hydraulicNetwork[c],
                                          0.0, maxNsteps,
                                          ntrial, psiTol, ETol, 0.001); 
      } else {
        hydraulicNetwork[c] = List::create(_["psisoil"] = psic,
                                           _["krhizomax"] = VGrhizo_kmaxc,_["nsoil"] = VG_nc,_["alphasoil"] = VG_alphac,
                                           _["krootmax"] = sapFluidityDay*VCroot_kmaxc, _["rootc"] = VCroot_c[c], _["rootd"] = VCroot_d[c],
                                           _["kstemmax"] = sapFluidityDay*VCstem_kmax[c], _["stemc"] = VCstem_c[c], _["stemd"] = VCstem_d[c],
                                           _["PLCstem"] = NumericVector::create(0.0));
        supply[c] = supplyFunctionNetworkStem1(hydraulicNetwork[c],
                                               0.0, maxNsteps,
                                               ntrial, psiTol, ETol, 0.001); 
      }
    } else {
      // Plant water pools
      
      //Copy soil water potentials from pools
      NumericMatrix psiSoilM(numCohorts, nlayers);
      List soil_pool = clone(soil);
      NumericVector Ws_pool = soil_pool["W"];
      for(int j = 0; j<numCohorts;j++) {
        //Copy values of soil moisture from pool of cohort j
        for(int l = 0; l<nlayers;l++) Ws_pool[l] = Wpool(j,l);
        //Calculate soil water potential
        psiSoilM(j,_) = psi(soil_pool, soilFunctions);
      }
      
      //Determine connected layers (non-zero fine root abundance)
      NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
      LogicalMatrix layerConnectedCoh(numCohorts, nlayers);
      
      nlayerscon[c] = 0;
      for(int j=0; j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          if((V(c,l)>0.0) && (RHOPcoh(j,l)>0.0)) {
            layerConnectedCoh(j,l)= true;
            nlayerscon[c]=nlayerscon[c] + 1;
          } else {
            layerConnectedCoh(j,l) = false;
          }
        }
      }
      if(nlayerscon[c]==0) stop("Plant cohort not connected to any soil layer!");
      //Store in list
      layerConnectedPools[c] = layerConnectedCoh; 
      
      // Copy values from connected layers
      NumericVector Vc = NumericVector(nlayerscon[c]);
      NumericVector VCroot_kmaxc = NumericVector(nlayerscon[c]);
      NumericVector VGrhizo_kmaxc = NumericVector(nlayerscon[c]);
      NumericVector psic = NumericVector(nlayerscon[c]);
      NumericVector VG_nc = NumericVector(nlayerscon[c]);
      NumericVector VG_alphac= NumericVector(nlayerscon[c]);
      int cnt=0;
      for(int j=0; j<numCohorts;j++) {
        for(int l=0;l<nlayers;l++) {
          if(layerConnectedCoh(j,l)) {
            Vc[cnt] = V(c,l)*RHOPcoh(j,l);
            VCroot_kmaxc[cnt] = VCroot_kmax(c,l)*RHOPcoh(j,l);
            VGrhizo_kmaxc[cnt] = VGrhizo_kmax(c,l)*RHOPcoh(j,l);
            psic[cnt] = psiSoilM(j,l);
            VG_nc[cnt] = VG_n[l];
            VG_alphac[cnt] = VG_alpha[l];
            cnt++;
          }
        }
      }
      
      //Build supply function networks 
      if(!capacitance) {
        hydraulicNetwork[c] = List::create(_["psisoil"] = psic,
                                           _["krhizomax"] = VGrhizo_kmaxc,_["nsoil"] = VG_nc,_["alphasoil"] = VG_alphac,
                                           _["krootmax"] = sapFluidityDay*VCroot_kmaxc, _["rootc"] = VCroot_c[c], _["rootd"] = VCroot_d[c],
                                           _["kstemmax"] = sapFluidityDay*VCstem_kmax[c], _["stemc"] = VCstem_c[c], _["stemd"] = VCstem_d[c],
                                           _["kleafmax"] = sapFluidityDay*VCleaf_kmax[c], _["leafc"] = VCleaf_c[c], _["leafd"] = VCleaf_d[c],
                                           _["PLCstem"] = NumericVector::create(StemPLCVEC[c],StemPLCVEC[c]));
        supply[c] = supplyFunctionNetwork(hydraulicNetwork[c],
                                          0.0, maxNsteps,
                                          ntrial, psiTol, ETol, 0.001); 
      } else {
        hydraulicNetwork[c] = List::create(_["psisoil"] = psic,
                                           _["krhizomax"] = VGrhizo_kmaxc,_["nsoil"] = VG_nc,_["alphasoil"] = VG_alphac,
                                           _["krootmax"] = sapFluidityDay*VCroot_kmaxc, _["rootc"] = VCroot_c[c], _["rootd"] = VCroot_d[c],
                                           _["kstemmax"] = sapFluidityDay*VCstem_kmax[c], _["stemc"] = VCstem_c[c], _["stemd"] = VCstem_d[c],
                                           _["PLCstem"] = NumericVector::create(0.0));
        supply[c] = supplyFunctionNetworkStem1(hydraulicNetwork[c],
                                               0.0, maxNsteps,
                                               ntrial, psiTol, ETol, 0.001); 
      }
    }
  }
  

  
  //Transpiration and photosynthesis
  NumericMatrix K(numCohorts, nlayers);
  NumericVector Eplant(numCohorts, 0.0), Anplant(numCohorts, 0.0), Agplant(numCohorts, 0.0);
  NumericMatrix Rninst(numCohorts,ntimesteps);
  NumericMatrix dEdPInst(numCohorts, ntimesteps);
  NumericMatrix Qinst(numCohorts,ntimesteps);
  NumericMatrix Einst(numCohorts, ntimesteps);
  NumericMatrix Aninst(numCohorts, ntimesteps), Aginst(numCohorts, ntimesteps);
  NumericMatrix LeafPsiInst(numCohorts, ntimesteps), StemPsiInst(numCohorts, ntimesteps);
  NumericMatrix LeafSympPsiInst(numCohorts, ntimesteps), StemSympPsiInst(numCohorts, ntimesteps);
  NumericMatrix LeafRWCInst(numCohorts, ntimesteps), StemRWCInst(numCohorts, ntimesteps);
  NumericMatrix LeafSympRWCInst(numCohorts, ntimesteps), StemSympRWCInst(numCohorts, ntimesteps);
  NumericMatrix RootPsiInst(numCohorts, ntimesteps);
  NumericMatrix PWBinst(numCohorts, ntimesteps);
  NumericMatrix E_SL(numCohorts, ntimesteps);
  NumericMatrix E_SH(numCohorts, ntimesteps);
  NumericMatrix An_SL(numCohorts, ntimesteps), Ag_SL(numCohorts, ntimesteps);
  NumericMatrix An_SH(numCohorts, ntimesteps), Ag_SH(numCohorts, ntimesteps);
  NumericMatrix Psi_SL(numCohorts, ntimesteps);
  NumericMatrix Psi_SH(numCohorts, ntimesteps);
  NumericMatrix Ci_SL(numCohorts, ntimesteps);
  NumericMatrix Ci_SH(numCohorts, ntimesteps);
  NumericMatrix SWR_SL(numCohorts, ntimesteps);
  NumericMatrix SWR_SH(numCohorts, ntimesteps);
  NumericMatrix PAR_SL(numCohorts, ntimesteps);
  NumericMatrix PAR_SH(numCohorts, ntimesteps);
  NumericMatrix LWR_SL(numCohorts, ntimesteps);
  NumericMatrix LWR_SH(numCohorts, ntimesteps);
  NumericMatrix GSW_SH(numCohorts, ntimesteps);
  NumericMatrix GSW_SL(numCohorts, ntimesteps);
  NumericMatrix VPD_SH(numCohorts, ntimesteps);
  NumericMatrix VPD_SL(numCohorts, ntimesteps);
  NumericMatrix Temp_SH(numCohorts, ntimesteps);
  NumericMatrix Temp_SL(numCohorts, ntimesteps);
  NumericVector minLeafPsi(numCohorts,0.0), maxLeafPsi(numCohorts,-99999.0); 
  NumericVector maxGSW_SL(numCohorts,-99999.0), maxGSW_SH(numCohorts,-99999.0); 
  NumericVector minGSW_SL(numCohorts,99999.0), minGSW_SH(numCohorts,99999.0); 
  NumericVector maxTemp_SL(numCohorts,-99999.0), maxTemp_SH(numCohorts,-99999.0); 
  NumericVector minTemp_SL(numCohorts,99999.0), minTemp_SH(numCohorts,99999.0); 
  NumericVector minLeafPsi_SL(numCohorts,0.0), maxLeafPsi_SL(numCohorts,-99999.0); 
  NumericVector minLeafPsi_SH(numCohorts,0.0), maxLeafPsi_SH(numCohorts,-99999.0);
  NumericVector minStemPsi(numCohorts, 0.0), minRootPsi(numCohorts,0.0); //Minimum potentials experienced
  NumericMatrix minPsiRhizo(numCohorts, nlayers);
  if(numCohorts>0) std::fill(minPsiRhizo.begin(), minPsiRhizo.end(), 0.0);
  NumericMatrix PLC(numCohorts, ntimesteps);
  NumericVector PLCm(numCohorts), RWCsm(numCohorts), RWClm(numCohorts),RWCssm(numCohorts), RWClsm(numCohorts);
  NumericVector dEdPm(numCohorts);
  NumericVector PWB(numCohorts,0.0);
  
  IntegerVector iPMSunlit(numCohorts,0), iPMShade(numCohorts,0); //Initial values set to closed stomata
    
  List outPhotoSunlit(numCohorts);
  List outPhotoShade(numCohorts);
  List outPMSunlit(numCohorts);
  List outPMShade(numCohorts);
  outPhotoSunlit.attr("names") = above.attr("row.names");
  outPhotoShade.attr("names") = above.attr("row.names");
  outPMSunlit.attr("names") = above.attr("row.names");
  outPMShade.attr("names") = above.attr("row.names");
  
  List lwrExtinctionList(ntimesteps);
  
  for(int n=0;n<ntimesteps;n++) { //Time loop
    //Longwave radiation
    List lwrExtinction = longwaveRadiationSHAW(LAIme, LAImd, LAImx, 
                                               lwdr[n], Tsoil[0], Tair);
    lwrExtinctionList[n] = lwrExtinction;
    net_LWR_soil[n] = lwrExtinction["Lnet_ground"];
    net_LWR_can[n]= lwrExtinction["Lnet_canopy"];
    NumericMatrix Lnet_cohort_layer = lwrExtinction["Lnet_cohort_layer"];
    
    //Retrieve radiation absorbed for the current time step
    NumericVector absPAR_SL_COH = abs_PAR_SL_COH_list[n];
    NumericVector absPAR_SH_COH = abs_PAR_SH_COH_list[n];
    NumericVector absSWR_SL_COH = abs_SWR_SL_COH_list[n];
    NumericVector absSWR_SH_COH = abs_SWR_SH_COH_list[n];
    NumericMatrix absSWR_SL_ML = abs_SWR_SL_ML_list[n];
    NumericMatrix absSWR_SH_ML = abs_SWR_SH_ML_list[n];

    for(int c=0;c<numCohorts;c++) { //Plant cohort loop
      
      //default values
      dEdPInst(c,n) = 0.0;
      Einst(c,n) = 0.0;
      Aginst(c,n) = 0.0;
      Aninst(c,n) = 0.0;
      
      if(LAIphe[c]>0.0) { //Process transpiration and photosynthesis only if there are some leaves
        PAR_SL(c,n) = absPAR_SL_COH[c];
        PAR_SH(c,n) = absPAR_SH_COH[c];
        SWR_SL(c,n) = absSWR_SL_COH[c];
        SWR_SH(c,n) = absSWR_SH_COH[c];
        // for(int j=0;j<ncanlayers;j++) Rcout<< n << " "<< c<< " "<<j<<" " << Lnet_cohort_layer(j,c)<<"\n";
        // Rcout<< n << " "<< c<< " LAIsl: " << LAI_SL[c]<< " LAIsh: " << LAI_SH[c]<< " LWRnet: "<< sum(Lnet_cohort_layer(_,c))<<" "<< sum(Lnet_cohort_layer(_,c)*fsunlit)<< " "<<sum(Lnet_cohort_layer(_,c)*(1.0 - fsunlit))<<"\n";
        LWR_SL(c,n) = sum(Lnet_cohort_layer(_,c)*fsunlit);
        LWR_SH(c,n) = sum(Lnet_cohort_layer(_,c)*(1.0 - fsunlit));
        
        //NumericVector PLCStemPrev = NumericVector::create(StemPLCVEC[c],StemPLCVEC[c]);
        //NumericVector psiStemPrev = NumericVector::create(Stem1PsiVEC[c],Stem2PsiVEC[c]);
        //double LeafPsiPrev = LeafPsiVEC[c];

        // Rcout<<c<<" E "<<EinstPrev<<" PR "<< RootPsiPrev<<" PL "<<LeafPsiPrev<< " PS "<<psiStemPrev[0]<< " "<<rwcsleafPrev<< " "<<RWCStemPrev[0]<<"\n";
        NumericVector fittedE, dEdP;
        NumericVector LeafPsi, psiRootCrown;
        
        //Retrieve supply functions
        List sFunctionBelow, sFunctionAbove;
        if(!capacitance) {
          sFunctionAbove = supply[c];
          sFunctionBelow = supply[c];
        } else {
          double psiPLCStem = apoplasticWaterPotential(1.0-StemPLCVEC[c], VCstem_c[c], VCstem_d[c]);
          double psiRootCrownFake = std::min(0.0,E2psiXylemUp(EinstVEC[c], Stem1PsiVEC[c],VCstem_kmax[c]*2.0, VCstem_c[c], VCstem_d[c], psiPLCStem));
          if(NumericVector::is_na(psiRootCrownFake)) psiRootCrownFake = 0.0;
          double psiFineRootFake= std::min(0.0,E2psiXylemUp(EinstVEC[c], psiRootCrownFake,VCroot_kmax_sum[c], VCroot_c[c], VCroot_d[c]));
          if(NumericVector::is_na(psiFineRootFake)) psiFineRootFake = 0.0;
          // Rcout<< c << " EinstVEC[c] "<< EinstVEC[c] << " Stem1PsiVEC[c] "<< Stem1PsiVEC[c]<<" psiFineRootFake "<< psiFineRootFake << " psiRootCrownFake "<< psiRootCrownFake<<"\n";
          // sFunctionAbove = supplyAboveground[c];
          double sapFluidityBelow = 1.0/waterDynamicViscosity(Tsoil[0]);
          double sapFluidityAbove = 1.0/waterDynamicViscosity(Tair[iLayerCohort[c]]);
          List hn = List::create(_["krootmax"] = NumericVector::create(sapFluidityBelow*VCroot_kmax_sum[c]), _["rootc"] = VCroot_c[c], _["rootd"] = VCroot_d[c],
                                 _["kstemmax"] = sapFluidityAbove*VCstem_kmax[c], _["stemc"] = VCstem_c[c], _["stemd"] = VCstem_d[c],
                                 _["kleafmax"] = sapFluidityAbove*VCleaf_kmax[c], _["leafc"] = VCleaf_c[c], _["leafd"] = VCleaf_d[c],
                                 _["PLCstem"] = NumericVector::create(StemPLCVEC[c]));
          
          sFunctionAbove = supplyFunctionFineRootLeaf(psiFineRootFake, hn, 
                                                      0.0, maxNsteps, 
                                                      ETol, 0.001);
          
          sFunctionBelow = supply[c];
        }

        //Determine turgor loss point (as proxy of stomatal closure)
        double psiTlp = turgorLossPoint(LeafPI0[c], LeafEPS[c]);
        
        //Retrieve transpiration, LeafPsi and dEdP vectors
        fittedE = sFunctionAbove["E"];
        dEdP = sFunctionAbove["dEdP"];
        LeafPsi = sFunctionAbove["psiLeaf"];
        
        //Get info from sFunctionAbove
        psiRootCrown = sFunctionAbove["psiRootCrown"];
        
        if(fittedE.size()>0) {
          //Photosynthesis function for sunlit and shade leaves
          List PMSunlit = profitMaximization2(sFunctionAbove, iPMSunlit[c], 
                                          Cair[iLayerSunlit[c]], Patm,
                                          Tair[iLayerSunlit[c]], VPair[iLayerSunlit[c]], zWind[iLayerSunlit[c]], 
                                          SWR_SL(c,n), LWR_SL(c,n), irradianceToPhotonFlux(PAR_SL(c,n)), 
                                          Vmax298SL[c], Jmax298SL[c], leafWidth[c], LAI_SL[c],
                                          Gswmin[c], Gswmax[c]);
          NumericVector photoSunlit = PMSunlit["photosynthesisFunction"];
          iPMSunlit[c] = PMSunlit["iMaxProfit"];
          List PMShade = profitMaximization2(sFunctionAbove, iPMShade[c], 
                                             Cair[iLayerShade[c]], Patm,
                                             Tair[iLayerShade[c]], VPair[iLayerShade[c]], zWind[iLayerShade[c]], 
                                             SWR_SH(c,n), LWR_SH(c,n), irradianceToPhotonFlux(PAR_SH(c,n)), 
                                             Vmax298SH[c], Jmax298SH[c], leafWidth[c], LAI_SH[c],
                                             Gswmin[c], Gswmax[c]);          
          NumericVector photoShade = PMShade["photosynthesisFunction"];
          iPMShade[c] = PMShade["iMaxProfit"];
          
          
          
          //Store?
          if(!IntegerVector::is_na(stepFunctions)) {
            if(n==stepFunctions) {
              outPhotoSunlit[c] = leafPhotosynthesisFunction2(fittedE, LeafPsi, Cair[iLayerSunlit[c]], Patm,
                                                              Tair[iLayerSunlit[c]], VPair[iLayerSunlit[c]], 
                                                              zWind[iLayerSunlit[c]], 
                                                              SWR_SL(c,n), LWR_SL(c,n), 
                                                              irradianceToPhotonFlux(PAR_SL(c,n)), 
                                                              Vmax298SL[c], 
                                                              Jmax298SL[c], 
                                                              leafWidth[c], LAI_SL[c]);
              outPhotoShade[c] = leafPhotosynthesisFunction2(fittedE, LeafPsi, Cair[iLayerShade[c]], Patm,
                                                             Tair[iLayerShade[c]], VPair[iLayerShade[c]], 
                                                             zWind[iLayerShade[c]], 
                                                             SWR_SH(c,n), LWR_SH(c,n), 
                                                             irradianceToPhotonFlux(PAR_SH(c,n)),
                                                             Vmax298SH[c], 
                                                             Jmax298SH[c], 
                                                             leafWidth[c], LAI_SH[c]);
              outPMSunlit[c] = PMSunlit;
              outPMShade[c] = PMShade;
            }
          }
          // Rcout<<iPMSunlit[c]<<" "<<iPMShade[c] <<" "<<GwSunlit[iPMSunlit[c]]<<" "<<GwShade[iPMShade[c]]<<" "<<fittedE[iPMSunlit[c]]<<" "<<fittedE[iPMShade[c]]<<"\n";
          //Get leaf status
          E_SH(c,n) = fittedE[iPMShade[c]];
          E_SL(c,n) = fittedE[iPMSunlit[c]];
          Psi_SH(c,n) = LeafPsi[iPMShade[c]];
          Psi_SL(c,n) = LeafPsi[iPMSunlit[c]];
          An_SH(c,n) = photoShade["NetPhotosynthesis"];
          An_SL(c,n) = photoSunlit["NetPhotosynthesis"];
          Ag_SH(c,n) = photoShade["GrossPhotosynthesis"];
          Ag_SL(c,n) = photoSunlit["GrossPhotosynthesis"];
          Ci_SH(c,n) = photoShade["Ci"];
          Ci_SL(c,n) = photoSunlit["Ci"];
          GSW_SH(c,n)= photoShade["Gsw"];
          GSW_SL(c,n)= photoSunlit["Gsw"];
          VPD_SH(c,n)= photoShade["LeafVPD"];
          VPD_SL(c,n)= photoSunlit["LeafVPD"];
          Temp_SH(c,n)= photoShade["LeafTemperature"];
          Temp_SL(c,n)= photoSunlit["LeafTemperature"];

          //Scale photosynthesis
          double Agsum = Ag_SL(c,n)*LAI_SL[c] + Ag_SH(c,n)*LAI_SH[c];
          double Ansum = An_SL(c,n)*LAI_SL[c] + An_SH(c,n)*LAI_SH[c];
          Aginst(c,n) = (1e-6)*12.01017*Agsum*tstep;
          Aninst(c,n) = (1e-6)*12.01017*Ansum*tstep;
          
          //Average flow from sunlit and shade leaves
          double Eaverage = (fittedE[iPMSunlit[c]]*LAI_SL[c] + fittedE[iPMShade[c]]*LAI_SH[c])/(LAI_SL[c] + LAI_SH[c]);
          
          
          //Find iPM for  flow corresponding to the  average flow
          double absDiff = 99999999.9;
          int iPM = -1;
          for(int k=0;k<fittedE.size();k++){ //Only check up to the size of fittedE
            double adk = std::abs(fittedE[k]-Eaverage);
            if(adk<absDiff) {
              absDiff = adk;
              iPM = k;
            }
          }
          if(iPM==-1) {
            Rcout<<"\n iPM -1! Eaverage="<< Eaverage << " fittedE.size= "<< fittedE.size()<<" iPMSunlit[c]="<< iPMSunlit[c]<< " fittedE[iPMSunlit[c]]="<<fittedE[iPMSunlit[c]]<<" iPMShade[c]="<<iPMShade[c]<<" fittedE[iPMShade[c]]="<<fittedE[iPMShade[c]]<<"\n";
            stop("");
          }

          //Store instantaneous total conductance
          dEdPInst(c,n) = dEdP[iPM];
          
          //Store instantaneous flow and leaf water potential
          EinstVEC[c] = fittedE[iPM];
          LeafPsiVEC[c] = LeafPsi[iPM];
          RootCrownPsiVEC[c] = psiRootCrown[iPM]; 
          
          //Scale from instantaneous flow to water volume in the time step
          Einst(c,n) = fittedE[iPM]*0.001*0.01802*LAIphe[c]*tstep; 
          
          NumericVector Esoilcn(nlayerscon[c],0.0);
          NumericVector ElayersVEC(nlayerscon[c],0.0);
          
          
          //Get info from sFunctionBelow (this will be different depending on wether capacitance is considered)
          NumericMatrix ERhizo = Rcpp::as<Rcpp::NumericMatrix>(sFunctionBelow["ERhizo"]);
          NumericMatrix RhizoPsi = Rcpp::as<Rcpp::NumericMatrix>(sFunctionBelow["psiRhizo"]);
          
          if(!capacitance) {
            //Store steady state stem and rootcrown and root surface water potential values
            NumericMatrix newStemPsi = Rcpp::as<Rcpp::NumericMatrix>(sFunctionAbove["psiStem"]);
            Stem1PsiVEC[c] = newStemPsi(iPM,0); 
            Stem2PsiVEC[c] = newStemPsi(iPM,1);
            for(int lc=0;lc<nlayerscon[c];lc++) {
              ElayersVEC[lc] = ERhizo(iPM,lc)*tstep; //Scale according to the time step
            }
            
            //Copy RhizoPsi and from connected layers to RhizoPsi from soil layers
            copyRhizoPsi(c,iPM, 
                         RhizoPsi, RhizoPsiMAT,
                         layerConnected, 
                         RHOP, layerConnectedPools,
                         VCroot_c, VCroot_d,  
                         plantWaterPools);
            
            StemSympPsiVEC[c] = Stem1PsiVEC[c]; //Stem symplastic compartment coupled with apoplastic compartment
            LeafSympPsiVEC[c] = LeafPsiVEC[c]; //Leaf symplastic compartment coupled with apoplastic compartment
            
            // Store the PLC corresponding to stem1 water potential
            if(cavitationRefill!="total") {
              StemPLCVEC[c] = std::max(StemPLCVEC[c], 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c])); 
            } else { //Immediate refilling
              StemPLCVEC[c] = 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c]); 
            }
            
          } else {
            //Store steady state stem2 water potential
            NumericVector newStemPsi2 = Rcpp::as<Rcpp::NumericVector>(sFunctionAbove["psiStem2"]);
            Stem2PsiVEC[c] = newStemPsi2[iPM];
            
            NumericVector newStemPsi1 = Rcpp::as<Rcpp::NumericVector>(sFunctionBelow["psiStem1"]);

            int iPMB = -1;
            
            //TO DO: stem segment water balance
            //Estimate current apoplastic and symplastic volumes
            // NOTE: Vsapwood and Vleaf are in l·m-2
            double VLeafSymp_mmolmax = 1000.0*((Vleaf[c]*(1.0-LeafAF[c]))/0.018); //mmol·m-2
            double VStemSymp_mmolmax = 1000.0*((Vsapwood[c]*(1.0-StemAF[c]))/0.018); //mmol·m-2
            //Substract from maximum apoplastic compartment embolized conduits
            double VStemApo_mmolmax = 1000.0*((Vsapwood[c]*StemAF[c])/0.018); //mmol·m-2
            double RWCLeafSymp = symplasticRelativeWaterContent(LeafSympPsiVEC[c], LeafPI0[c], LeafEPS[c]); //mmol·m-2
            double RWCStemSymp = symplasticRelativeWaterContent(StemSympPsiVEC[c], StemPI0[c], StemEPS[c]); //mmol·m-2
            double VLeafSymp_mmol = VLeafSymp_mmolmax * RWCLeafSymp;
            double VStemSymp_mmol = VStemSymp_mmolmax * RWCStemSymp;
            double Vcav = 0.0;
            //Perform water balance
            // Rcout<<"\n"<<c<<" Before - iPM " << iPM<< " EinstVEC[c]: "<< EinstVEC[c]<<" Vol: "<<VStemApo_mmol<<" RWC:"<< RWCStemApo <<" Psi: "<< Stem1PsiVEC[c]<< " LeafPsiVEC[c]: "<<LeafPsiVEC[c]<<"\n";
            for(double scnt=0.0; scnt<tstep;scnt += 1.0) {
              //Find flow corresponding to Stem1PsiVEC[c]
              //Find iPM for water potential corresponding to the current water potential
              double absDiff = 99999999.9;
              iPMB = -1;
              for(int k=0;k<newStemPsi1.size();k++){ //Only check up to the size of fittedE
                double adk = std::abs(newStemPsi1[k]-Stem1PsiVEC[c]);
                if(adk<absDiff) {
                  absDiff = adk;
                  iPMB = k;
                }
              }
              if(iPMB==-1) {
                Rcout<<"\n Stem1PsiVEC[c]="<< Stem1PsiVEC[c] << " newStemPsi1.size= "<< newStemPsi1.size()<<"\n";
                stop("iPMB = -1");
              }
              // Stem1PsiVEC[c] = newStemPsi1[iPMB];
                
              //Add flow from soil to ElayersVEC
              for(int lc=0;lc<nlayerscon[c];lc++) ElayersVEC[lc] += ERhizo(iPMB,lc); 
              
              //Calculate stem and leaf lateral flows
              double Flatstem = (StemSympPsiVEC[c] - Stem1PsiVEC[c])*klatstem;
              double Flatleaf = (LeafSympPsiVEC[c] - LeafPsiVEC[c])*klatleaf;


              //Leaf symplastic water balance
              VLeafSymp_mmol += (-Flatleaf);
              RWCLeafSymp = VLeafSymp_mmol/VLeafSymp_mmolmax;
              LeafSympPsiVEC[c] = symplasticWaterPotential(std::min(1.0,RWCLeafSymp), LeafPI0[c], LeafEPS[c]);
              if(NumericVector::is_na(LeafSympPsiVEC[c]))  LeafSympPsiVEC[c] = -40.0;
              
              //Stem symplastic water balance
              VStemSymp_mmol += (-Flatstem);
              RWCStemSymp = VStemSymp_mmol/VStemSymp_mmolmax;
              StemSympPsiVEC[c] = symplasticWaterPotential(std::min(1.0,RWCStemSymp), StemPI0[c], StemEPS[c]);
              if(NumericVector::is_na(StemSympPsiVEC[c]))  StemSympPsiVEC[c] = -40.0;
              
              //Stem apoplastic water balance
              double Vchange = (Flatstem + sum(ERhizo(iPMB,_)) - (EinstVEC[c] - Flatleaf)) + Vcav;
              
              Stem1PsiVEC[c] = Stem1PsiVEC[c] + eps_xylem*(Vchange/VStemApo_mmolmax);
              
              // VStemApo_mmol += (Flatstem + sum(ERhizo(iPMB,_)) - (EinstVEC[c] - Flatleaf));
              // RWCStemApo = VStemApo_mmol/VStemApo_mmolmax;
              // Stem1PsiVEC[c] = apoplasticWaterPotential(std::min(1.0,RWCStemApo), VCstem_c[c], VCstem_d[c]);
              // if(NumericVector::is_na(Stem1PsiVEC[c]))  Stem1PsiVEC[c] = -40.0;
              

              //Recalculate PLC and calculate volume corresponding to new cavitation
              double plc_old = StemPLCVEC[c];
              if(cavitationRefill!="total") {
                StemPLCVEC[c] = std::max(StemPLCVEC[c], 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c])); 
                Vcav = VStemApo_mmolmax*(StemPLCVEC[c]-plc_old);
              } else { //Immediate refilling
                StemPLCVEC[c] = 1.0 - xylemConductance(Stem1PsiVEC[c], 1.0, VCstem_c[c], VCstem_d[c]); 
                Vcav = 0.0;
              }

            }
            
            // Rcout<<c<<" after - EinstVEC: "<<EinstVEC[c] << " RWCStemApo: " << RWCStemApo << "  Stem1PsiVEC:"<< Stem1PsiVEC[c]<<" StemSympPsiVEC: "<< StemSympPsiVEC[c]<<" LeafSympPsiVEC: "<< LeafSympPsiVEC[c] <<"\n";

            //Copy RhizoPsi and from connected layers to RhizoPsi from soil layers
            copyRhizoPsi(c,iPMB, 
                         RhizoPsi, RhizoPsiMAT,
                         layerConnected, 
                         RHOP, layerConnectedPools,
                         VCroot_c, VCroot_d,  
                         plantWaterPools);
          }
          
          //Scale soil water extracted from leaf to cohort level
          for(int lc=0;lc<nlayerscon[c];lc++) {
            Esoilcn[lc] = ElayersVEC[lc]*0.001*0.01802*LAIphe[c]; //Scale from flow to water volume in the time step
          }
          
          //Balance between extraction and transpiration
          PWBinst(c,n) = sum(Esoilcn) - Einst(c,n);
          
          //Add step transpiration to daily plant cohort transpiration
          Eplant[c] += Einst(c,n);
          Anplant[c] += Aninst(c,n);
          Agplant[c] += Aginst(c,n);
          //Add PWB
          PWB[c] += PWBinst(c,n); 
          

          
          //Copy transpiration and from connected layers to transpiration from soil layers
          if(!plantWaterPools) {
            int cl = 0;
            for(int l=0;l<nlayers;l++) {
              if(layerConnected(c,l)) {
                SoilWaterExtract(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
                soilLayerExtractInst(l,n) += Esoilcn[cl];
                cl++;
              } 
            }
          } else {
            NumericMatrix RHOPcoh = Rcpp::as<Rcpp::NumericMatrix>(RHOP[c]);
            LogicalMatrix layerConnectedCoh = Rcpp::as<Rcpp::LogicalMatrix>(layerConnectedPools[c]);
            int cl = 0;
            for(int j = 0;j<numCohorts;j++) {
              for(int l=0;l<nlayers;l++) {
                if(layerConnectedCoh(j,l)) {
                  SoilWaterExtract(c,l) += Esoilcn[cl]; //Add to cummulative transpiration from layers
                  soilLayerExtractInst(l,n) += Esoilcn[cl];
                  //Apply extraction to soil layer
                  if(modifyInput) Wpool(j,l) = Wpool(j,l) - (Esoilcn[cl]/(Water_FC[l]*poolProportions[j])); //Apply extraction from pools
                  cl++;
                }
              }
            }
          }
        } else {
          if(verbose) Rcout<<"NS!";
          Psi_SH(c,n) = NA_REAL;
          Psi_SL(c,n) = NA_REAL;
          GSW_SH(c,n)= NA_REAL;
          GSW_SL(c,n)= NA_REAL;
          VPD_SH(c,n)= NA_REAL;
          VPD_SL(c,n)= NA_REAL;
          Temp_SH(c,n)= NA_REAL;
          Temp_SL(c,n)= NA_REAL;
        }        
      } else if(N[c]>0.0) { //Cohorts with living individuals but no LAI should be in equilibrium with soil (i.e. no transpiration)
        List sFunctionBelow = supply[c];
        NumericVector  psiRootCrown = sFunctionBelow["psiRootCrown"];
        RootCrownPsiVEC[c] = psiRootCrown[0];
        if(!capacitance) {
          NumericVector  psiStem1 = sFunctionBelow["psiStem"];
          Stem1PsiVEC[c] = psiStem1[0];
          StemSympPsiVEC[c] = psiStem1[0];
          NumericVector LeafPsi = sFunctionBelow["psiLeaf"];
          LeafPsiVEC[c] = LeafPsi[0];
          LeafSympPsiVEC[c] = LeafPsi[0];
        } else {
          NumericVector  psiStem1 = sFunctionBelow["psiStem1"];
          Stem1PsiVEC[c] = psiStem1[0];
          StemSympPsiVEC[c] = psiStem1[0];
          LeafPsiVEC[c] = Stem1PsiVEC[c];
          LeafSympPsiVEC[c] = StemSympPsiVEC[c];
        }
      }
      
      if(N[c]>0.0) {
        //Store (for output) instantaneous leaf, stem and root potential, plc and rwc values
        PLC(c,n) = StemPLCVEC[c];
        StemSympRWCInst(c,n) = symplasticRelativeWaterContent(StemSympPsiVEC[c], StemPI0[c], StemEPS[c]);
        LeafSympRWCInst(c,n) = symplasticRelativeWaterContent(LeafSympPsiVEC[c], LeafPI0[c], LeafEPS[c]);
        StemRWCInst(c,n) = StemSympRWCInst(c,n)*(1.0 - StemAF[c]) + apoplasticRelativeWaterContent(Stem1PsiVEC[c], VCstem_c[c], VCstem_d[c])*StemAF[c];
        LeafRWCInst(c,n) = LeafSympRWCInst(c,n)*(1.0 - LeafAF[c]) + apoplasticRelativeWaterContent(LeafPsiVEC[c], VCleaf_c[c], VCleaf_d[c])*LeafAF[c];
        StemPsiInst(c,n) = Stem1PsiVEC[c]; 
        LeafPsiInst(c,n) = LeafPsiVEC[c]; //Store instantaneous (average) leaf potential
        RootPsiInst(c,n) = RootCrownPsiVEC[c]; //Store instantaneous root crown potential
        LeafSympPsiInst(c,n) = LeafSympPsiVEC[c];
        StemSympPsiInst(c,n) = StemSympPsiVEC[c];
        
        //Store the minimum water potential of the day (i.e. mid-day)
        minGSW_SL[c] = std::min(minGSW_SL[c], GSW_SL(c,n));
        minGSW_SH[c] = std::min(minGSW_SH[c], GSW_SH(c,n));
        maxGSW_SL[c] = std::max(maxGSW_SL[c], GSW_SL(c,n));
        maxGSW_SH[c] = std::max(maxGSW_SH[c], GSW_SH(c,n));
        minTemp_SL[c] = std::min(minTemp_SL[c], Temp_SL(c,n));
        minTemp_SH[c] = std::min(minTemp_SH[c], Temp_SH(c,n));
        maxTemp_SL[c] = std::max(maxTemp_SL[c], Temp_SL(c,n));
        maxTemp_SH[c] = std::max(maxTemp_SH[c], Temp_SH(c,n));
        minLeafPsi_SL[c] = std::min(minLeafPsi_SL[c],Psi_SL(c,n));
        minLeafPsi_SH[c] = std::min(minLeafPsi_SH[c],Psi_SH(c,n));
        maxLeafPsi_SL[c] = std::max(maxLeafPsi_SL[c],Psi_SL(c,n));
        maxLeafPsi_SH[c] = std::max(maxLeafPsi_SH[c],Psi_SH(c,n));
        minLeafPsi[c] = std::min(minLeafPsi[c],LeafPsiInst(c,n));
        maxLeafPsi[c] = std::max(maxLeafPsi[c],LeafPsiInst(c,n));
        minStemPsi[c] = std::min(minStemPsi[c],StemPsiInst(c,n));
        minRootPsi[c] = std::min(minRootPsi[c],RootPsiInst(c,n));
        for(int l=0;l<nlayers;l++) {
          minPsiRhizo(c,l) = std::min(minPsiRhizo(c,l),RhizoPsiMAT(c,l));
        }
      }
    } //End of cohort loop
    
    //CANOPY AND SOIL ENERGY BALANCE
    
    //Soil latent heat (soil evaporation)
    //Latent heat (snow fusion) as J/m2/s
    double soilEvapStep = abs_SWR_soil[n]*(soilEvaporation/sum(abs_SWR_soil));
    double snowMeltStep = abs_SWR_soil[n]*(snowMelt/sum(abs_SWR_soil));
    if(sum(abs_SWR_soil)==0.0) { // avoid zero sums
      soilEvapStep = 0.0; 
      snowMeltStep = 0.0;
    }
    if(SWE>0.0) {
      abs_SWR_soil[n] = 0.0; //Set SWR absorbed by soil to zero (for energy balance) if snow pack is present 
    }
    double LEsoilevap = (1e6)*meteoland::utils_latentHeatVaporisation(Tsoil[0])*soilEvapStep/tstep;
    double LEsnow = (1e6)*(snowMeltStep*0.33355)/tstep; // 0.33355 = latent heat of fusion
    LEsoil_heat[n] = LEsoilevap + LEsnow;
    
    
    //Canopy evaporation (mm) in the current step
    double canEvapStep = canopyEvaporation*(abs_SWR_can[n]/sum(abs_SWR_can));

    //Canopy convective heat exchange
    double RAcan = aerodynamicResistance(canopyHeight,std::max(wind,1.0)); //Aerodynamic resistance to convective heat transfer
    Hcan_heat[n] = (meteoland::utils_airDensity(Tatm[n],Patm)*Cp_JKG*(Tcan[n]-Tatm[n]))/RAcan;
    
    if(!multiLayerBalance) {//Canopy balance assuming a single layer
      //Soil-canopy turbulent heat exchange
      double wind2m = windSpeedMassmanExtinction(200.0, wind, LAIcell, canopyHeight);
      double RAsoil = aerodynamicResistance(200.0, std::max(wind2m,1.0)); //Aerodynamic resistance to convective heat transfer from soil
      Hcansoil[n] = (meteoland::utils_airDensity(Tcan[n],Patm)*Cp_JKG*(Tcan[n]-Tsoil[0]))/RAsoil;
      //Latent heat (evaporation + transpiration)
      double LEwat = (1e6)*meteoland::utils_latentHeatVaporisation(Tcan[n])*(sum(Einst(_,n)) + canEvapStep)/tstep;
      LEcan_heat[n] = LEwat; 
      //Canopy temperature changes
      Ebal[n] = abs_SWR_can[n]+ net_LWR_can[n] - LEcan_heat[n] - Hcan_heat[n] - Hcansoil[n];
      double canopyAirThermalCapacity = meteoland::utils_airDensity(Tcan[n],Patm)*Cp_JKG;
      double canopyThermalCapacity =  canopyAirThermalCapacity + (0.5*(0.8*LAIcelllive + 1.2*LAIcell) + LAIcelldead)*thermalCapacityLAI; //Avoids zero capacity for winter deciduous
      double Tcannext = Tcan[n]+ std::max(-3.0, std::min(3.0, tstep*Ebal[n]/canopyThermalCapacity)); //Avoids changes in temperature that are too fast
      if(n<(ntimesteps-1)) Tcan[n+1] = Tcannext;
      for(int i=0;i<ncanlayers;i++) Tair[i] = Tcannext;
      
      
      //Soil energy balance including exchange with canopy
      Ebalsoil[n] = abs_SWR_soil[n] + net_LWR_soil[n] + Hcansoil[n] - LEsoil_heat[n]; //Here we use all energy escaping to atmosphere
      
      
      //Soil temperature changes
      NumericVector soilTchange = temperatureChange(dVec, Tsoil, sand, clay, Ws, Theta_FC, Ebalsoil[n]);
      for(int l=0;l<nlayers;l++) Tsoil[l] = Tsoil[l] + (soilTchange[l]*tstep);
      if(n<(ntimesteps-1)) Tsoil_mat(n+1,_)= Tsoil;
      
    } else { //Multilayer canopy balance
      double moistureAtm = 0.622*(vpatm/Patm)*meteoland::utils_airDensity(Tatm[n],Patm);
      double CO2Atm = 0.409*Catm*44.01; //mg/m3
        
      double tsubstep = tstep/((double) nsubsteps); 
      double maxTchange = 3.0/((double) nsubsteps);
      double maxMoistureChange = 0.001/((double)nsubsteps); //=0.16 kPa per step
      double maxCO2Change = 180.0/((double)nsubsteps); //= 10 ppm per step
      double deltaZ = (verticalLayerSize/100.0); //Vertical layer size in m
      DataFrame LWR_layer = Rcpp::as<Rcpp::DataFrame>(lwrExtinction["LWR_layer"]);
      NumericVector LWRnet_layer = LWR_layer["Lnet"];
      Ebal[n] = 0.0;
      LEcan_heat[n] = 0.0;
      NumericVector Tairnext(ncanlayers), LElayer(ncanlayers), absSWRlayer(ncanlayers), Rnlayer(ncanlayers), Hleaflayer(ncanlayers);
      NumericVector layerThermalCapacity(ncanlayers);
      NumericVector moistureET(ncanlayers), rho(ncanlayers), moistureLayer(ncanlayers), moistureLayernext(ncanlayers);
      NumericVector CO2An(ncanlayers), CO2Layer(ncanlayers), CO2Layernext(ncanlayers);
      for(int i=0;i<ncanlayers;i++) {
        rho[i] = meteoland::utils_airDensity(Tair[i],Patm);
        absSWRlayer[i] = sum(absSWR_SL_ML(i,_)) + sum(absSWR_SH_ML(i,_));
        //Radiation balance
        Rnlayer[i] = absSWRlayer[i] + LWRnet_layer[i];
        NumericVector pLayer = LAIme(i,_)/LAIphe; //Proportion of each cohort LAI in layer i
        //Instantaneous layer transpiration
        //from mmolH2O/m2/s to kgH2O/m2/s
        double ElayerInst = 0.001*0.01802*sum(LAIme(i,_)*(E_SL(_,n)*fsunlit[i] + E_SH(_,n)*(1.0-fsunlit[i])));
        //Assumes Layers contribute to evaporation proportionally to their LAI fraction
        double layerEvapInst = (canEvapStep/tstep)*(LAIpe[i]/LAIcellexpanded);
        //Estimate instantaneous mgCO2/m2 absorption for the layer, taking into account the proportion of sunlit and shade leaves of each cohort
        //from micro.molCO2/m2/s to mgCO2/m2/s
        double Anlayer =(1e-3)*44.01*sum(LAIme(i,_)*(An_SL(_,n)*fsunlit[i] + An_SH(_,n)*(1.0-fsunlit[i])));
        // 1000.0*(44.01/12.0)*sum(Aninst(_,n)*pLayer); 
        double LEwat = (1e6)*meteoland::utils_latentHeatVaporisation(Tair[i])*(ElayerInst + layerEvapInst);
        LElayer[i] = LEwat; //Energy spent in vaporisation
        LEcan_heat[n] = LElayer[i];
        // layerThermalCapacity[i] = (0.5*(0.8*LAIcelllive + 1.2*LAIcell) + LAIcelldead)*thermalCapacityLAI/((double) ncanlayers);
        layerThermalCapacity[i] =  (0.5*(0.8*LAIpx[i] + 1.2*LAIpe[i]) + LAIpd[i])*thermalCapacityLAI; //Avoids zero capacity for winter deciduous
        
        moistureLayer[i] = 0.622*(VPair[i]/Patm)*rho[i]; //kg water vapour/m3
        moistureET[i] = (ElayerInst + layerEvapInst)/(deltaZ); //kg water vapour /m3/s
        
        CO2Layer[i] = 0.409*Cair[i]*44.01; //mg/m3
        CO2An[i] = -1.0*Anlayer/(deltaZ); //mg/m3/s
        // Rcout<<n<< " "<< i<< " - Rn: "<<Rnlayer[i]<<" LE: "<<LElayer[i]<<" Hleaf: "<<Hleaflayer[i]<< " Tini: "<< Tair[i]<<"\n";
        // Rcout<<n<< " "<< i<< " - moistureET: "<<moistureET[i]<<" moistureLayer: "<<moistureLayer[i]<<" CO2An: "<<CO2An[i]<< " CO2Layer: "<< CO2Layer[i]<<"\n";
      }
      //Add soil moisture evaporation
      moistureET[0] += soilEvapStep/(deltaZ*tstep); //kg/m3/s
      
      
      for(int s=0;s<nsubsteps;s++) {
        double RAsoil = aerodynamicResistance(200.0, std::max(zWind[0],1.0)); //Aerodynamic resistance to convective heat transfer from soil
        double Hcansoils = Cp_JKG*meteoland::utils_airDensity(Tair[0],Patm)*(Tair[0]-Tsoil[0])/RAsoil;
        // double Hcan_heats = (meteoland::utils_airDensity(Tatm[n],Patm)*Cp_JKG*(Tair[ncanlayers-1]-Tatm[n]))/RAcan;
        for(int i=0;i<ncanlayers;i++) {
          double deltaH = 0.0;
          double deltaMoisture = 0.0;
          double deltaCO2 = 0.0;
          // double Hlayers = 0.0;
          //Add turbulent heat flow (positive gradient when temperature is larger above)
          if(i==0) { //Lower layer
            deltaH -= (Cp_JKG*rho[i]*(Tair[i+1] - Tair[i])*uw[i])/(deltaZ*dU[i]);
            deltaH -= Hcansoils;
            deltaMoisture -= ((moistureLayer[i+1] - moistureLayer[i])*uw[i])/(deltaZ*dU[i]);
            deltaCO2 -= ((CO2Layer[i+1] - CO2Layer[i])*uw[i])/(deltaZ*dU[i]);
          } else if((i > 0) && (i<(ncanlayers-1))) { //Intermediate layers
            deltaH -= (Cp_JKG*rho[i]*(Tair[i+1] - Tair[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaMoisture -= ((moistureLayer[i+1] - moistureLayer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaCO2 -= ((CO2Layer[i+1] - CO2Layer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
          } else if(i==(ncanlayers-1)){ //Upper layer
            deltaH -= (Cp_JKG*rho[i]*(Tatm[n] - Tair[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaMoisture -= ((moistureAtm - moistureLayer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
            deltaCO2 -= ((CO2Atm - CO2Layer[i-1])*uw[i])/(2.0*deltaZ*dU[i]);
          }
          Hleaflayer[i] = 0.0;
          for(int c=0;c<numCohorts;c++) {
            double gHa = 0.189*pow(std::max(zWind[i],0.1)/(leafWidth[c]*0.0072), 0.5);
            double Hsunlit = 2.0*Cp_Jmol*rho[i]*(Temp_SL(c, n)-Tair[i])*gHa;
            double Hshade = 2.0*Cp_Jmol*rho[i]*(Temp_SH(c, n)-Tair[i])*gHa;
            // Rcout<<c<<" " << Hsunlit<< " "<<Hshade<<" \n";
            Hleaflayer[i] +=(Hsunlit*fsunlit[i] + Hshade*(1.0-fsunlit[i]))*LAIme(i,c);
          }
          double EbalLayer = Rnlayer[i] - LElayer[i] + Hleaflayer[i] + deltaH; 
          //Instantaneous changes in temperature due to internal energy balance
          double deltaT = EbalLayer/(rho[i]*Cp_JKG + layerThermalCapacity[i]); 
          
          // if(s==0) Rcout<<n<< " "<< i<< " "<< s <<" - Rn: "<<Rnlayer[i]<<" LE: "<<LElayer[i]<<" Hleaf: "<<Hleaflayer[i]<<" H: "<<Hlayers<< " Ebal: "<<EbalLayer<< " LTC: " << rholayer*Cp_JKG + layerThermalCapacity[i]<< " Tini: "<< Tair[i]<< " deltaT: "<<deltaT<<"\n";
          Tairnext[i] = Tair[i] +  std::max(-1.0*maxTchange, std::min(maxTchange, tsubstep*deltaT)); //Avoids changes in temperature that are too fast
          //Changes in water vapour
          moistureLayernext[i] = moistureLayer[i] + std::max(-1.0*maxMoistureChange, std::min(maxMoistureChange, tsubstep*(moistureET[i]+ deltaMoisture)));
          // if(i==0) Rcout<<n<< " "<< i<< " "<< s <<" - moisture: "<<moistureLayer[i]<<" delta: "<<deltaMoisture<<" next: "<< moistureLayernext[i]<<"\n";
          //Changes in CO2
          CO2Layernext[i] = CO2Layer[i] + std::max(-1.0*maxCO2Change, std::min(maxCO2Change, tsubstep*(CO2An[i] + deltaCO2)));
        }
        for(int i=0;i<ncanlayers;i++) {
          Tair[i] = Tairnext[i]; 
          moistureLayer[i] = moistureLayernext[i];
          CO2Layer[i] = CO2Layernext[i];
        }
        
        //Soil energy balance including exchange with canopy
        double Ebalsoils =  abs_SWR_soil[n] + net_LWR_soil[n]  - LEsoil_heat[n] + Hcansoils; //Here we use all energy escaping to atmosphere
        Ebalsoil[n] +=Ebalsoils;
        Hcansoil[n] +=Hcansoils;
        //Soil temperature changes
        NumericVector soilTchange = temperatureChange(dVec, Tsoil, sand, clay, Ws, Theta_FC, Ebalsoils);
        for(int l=0;l<nlayers;l++) Tsoil[l] = Tsoil[l] + (soilTchange[l]*tsubstep);
      }
      Hcansoil[n] = Hcansoil[n]/((double) nsubsteps);
      Ebalsoil[n] = Ebalsoil[n]/((double) nsubsteps);
      for(int i=0;i<ncanlayers;i++) {
        VPair[i] = ((moistureLayer[i]/rho[i])*Patm)/0.622;
        Cair[i] = CO2Layer[i]/(0.409*44.01);
        // Rcout<< n << " "<<i << " - " << moistureLayer[i]<< " "<< VPair[i]<<"\n";
      }
      // stop("kk");
      // Rcout<< n << " "<<abs_SWR_can[n]<< " = "<< sum(absSWRlayer)<<" = "<< (sum(absSWR_SL_ML) + sum(absSWR_SH_ML))<<" "<<net_LWR_can[n]<< " = "<< sum(LWRnet_layer)<<"\n";
      //Canopy energy balance
      Ebal[n] = abs_SWR_can[n]+ net_LWR_can[n] - LEcan_heat[n] - Hcan_heat[n] - Hcansoil[n];
      if(n<(ntimesteps-1)) {
        Tcan[n+1] = sum(Tair*LAIpx)/sum(LAIpx); 
        Tsoil_mat(n+1,_)= Tsoil;
      }
    }
    if(n<(ntimesteps-1)) for(int i=0;i<ncanlayers;i++) {
      Tcan_mat(n+1,i) = Tair[i];
      VPcan_mat(n+1,i) = VPair[i];
    }
  } //End of timestep loop

  //4z. Plant daily drought stress (from root collar mid-day water potential)
  NumericVector SoilExtractCoh(numCohorts,0.0);
  NumericVector DDS(numCohorts, 0.0), LFMC(numCohorts, 0.0);
  for(int c=0;c<numCohorts;c++) {
    SoilExtractCoh[c] =  sum(SoilWaterExtract(c,_));
    PLCm[c] = sum(PLC(c,_))/((double)PLC.ncol());
    RWCsm[c] = sum(StemRWCInst(c,_))/((double)StemRWCInst.ncol());
    RWClm[c] = sum(LeafRWCInst(c,_))/((double)LeafRWCInst.ncol());
    LFMC[c] = maxFMC[c]*((1.0/r635[c])*RWClm[c]+(1.0 - (1.0/r635[c]))*RWCsm[c]);
    dEdPm[c] = sum(dEdPInst(c,_))/((double)dEdPInst.ncol());  
    double maxConductance = maximumSoilPlantConductance(VGrhizo_kmax(c,_), VCroot_kmax(c,_), VCstem_kmax[c], VCleaf_kmax[c]);
    DDS[c] = Phe[c]*(1.0 - (dEdPm[c]/(sapFluidityDay*maxConductance)));
    
    if(cavitationRefill=="rate") {
      double SAmax = 10e4/Al2As[c]; //cm2·m-2 of leaf area
      double r = refillMaximumRate*std::max(0.0, (StemSympPsiVEC[c] + 1.5)/1.5);
      StemPLCVEC[c] = std::max(0.0, StemPLCVEC[c] - (r/SAmax));
    }
  }
  
  
  //B.3 - Substract  extracted water from soil moisture 
  if(modifyInput){
    for(int l=0;l<nlayers;l++) {
      Ws[l] = std::max(Ws[l] - (sum(soilLayerExtractInst(l,_))/Water_FC[l]),0.0);
    } 
    if(!plantWaterPools) { //copy soil to the pools of all cohorts
      for(int c=0;c<numCohorts;c++) {
        for(int l=0;l<nlayers;l++) {
          Wpool(c,l) = Ws[l];
        }
      }
    }
  }
  //6. Copy LAIexpanded for output
  NumericVector LAIcohort(numCohorts);
  for(int c=0;c<numCohorts;c++) {
    LAIcohort[c]= LAIphe[c];
  }
  LAIcohort.attr("names") = above.attr("row.names");
  
  DataFrame Tinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                      _["Tatm"] = Tatm, _["Tcan"] = Tcan, _["Tsoil"] = Tsoil_mat);
  Tcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  VPcan_mat.attr("dimnames") = List::create(seq(1,ntimesteps), seq(1,ncanlayers));
  DataFrame CEBinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                        _["SWRcan"] = abs_SWR_can, _["LWRcan"] = net_LWR_can,
                                        _["LEcan"] = LEcan_heat, _["Hcan"] = Hcan_heat, _["Ebalcan"] = Ebal);
  DataFrame SEBinst = DataFrame::create(_["SolarHour"] = solarHour, 
                                        _["Hcansoil"] = Hcansoil, _["LEsoil"] = LEsoil_heat, _["SWRsoil"] = abs_SWR_soil, _["LWRsoil"] = net_LWR_soil,
                                        _["Ebalsoil"] = Ebalsoil);
  List EB = List::create(_["Temperature"]=Tinst, _["CanopyEnergyBalance"] = CEBinst, _["SoilEnergyBalance"] = SEBinst,
                         _["TemperatureLayers"] = NA_REAL, _["VaporPressureLayers"] = NA_REAL);
  if(multiLayerBalance) {
    EB["TemperatureLayers"] = Tcan_mat;
    EB["VaporPressureLayers"] = VPcan_mat;
  }
  E_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  E_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Psi_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Psi_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  An_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  An_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ag_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ag_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  SWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  SWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PAR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PAR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LWR_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ci_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Ci_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GSW_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  GSW_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Temp_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SH.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  VPD_SL.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));

  Einst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  dEdPInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemSympPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemSympRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  RootPsiInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Aginst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  Aninst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PLC.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  LeafRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  StemRWCInst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  PWBinst.attr("dimnames") = List::create(above.attr("row.names"), seq(1,ntimesteps));
  minPsiRhizo.attr("dimnames") = List::create(above.attr("row.names"), seq(1,nlayers));
  soilLayerExtractInst.attr("dimnames") = List::create(seq(1,nlayers), seq(1,ntimesteps));
  for(int c=0;c<numCohorts;c++) {
    Vmax298SL[c] = Vmax298SL[c]/LAI_SL[c];
    Jmax298SL[c] = Jmax298SL[c]/LAI_SL[c];
    Vmax298SH[c] = Vmax298SH[c]/LAI_SH[c];
    Jmax298SH[c] = Jmax298SH[c]/LAI_SH[c];
  }
  DataFrame Sunlit = DataFrame::create(
    _["LAI"] = LAI_SL, 
    _["Vmax298"] = Vmax298SL,
    _["Jmax298"] = Jmax298SL,
    _["LeafPsiMin"] = minLeafPsi_SL, 
    _["LeafPsiMax"] = maxLeafPsi_SL, 
    _["GSWMin"] = minGSW_SL,
    _["GSWMax"] = maxGSW_SL,
    _["TempMin"] = minTemp_SL,
    _["TempMax"] = maxTemp_SL  
  );
  DataFrame Shade = DataFrame::create(
    _["LAI"] = LAI_SH, 
    _["Vmax298"] = Vmax298SH,
    _["Jmax298"] = Jmax298SH,
    _["LeafPsiMin"] = minLeafPsi_SH, 
    _["LeafPsiMax"] = maxLeafPsi_SH, 
    _["GSWMin"] = minGSW_SH,
    _["GSWMax"] = maxGSW_SH,
    _["TempMin"] = minTemp_SH,
    _["TempMax"] = maxTemp_SH  
  );
  Sunlit.attr("row.names") = above.attr("row.names");
  Shade.attr("row.names") = above.attr("row.names");
  
  List ShadeInst = List::create(
    _["Abs_SWR"] = SWR_SH,
    _["Abs_PAR"]=PAR_SH,
    _["Net_LWR"] = LWR_SH,
    _["Ag"] = Ag_SH,
    _["An"] = An_SH,
    _["Ci"] = Ci_SH,
    _["E"] = E_SH,
    _["Gsw"] = GSW_SH,
    _["VPD"] = VPD_SH,
    _["Temp"] = Temp_SH,
    _["Psi"] = Psi_SH);
  List SunlitInst = List::create(
    _["Abs_SWR"]=SWR_SL,
    _["Abs_PAR"]=PAR_SL,
    _["Net_LWR"] = LWR_SL,
    _["Ag"] = Ag_SL,
    _["An"] = An_SL,
    _["Ci"] = Ci_SL,
    _["E"] = E_SL,
    _["Gsw"] = GSW_SL,
    _["VPD"] = VPD_SL,
    _["Temp"] = Temp_SL,
    _["Psi"] = Psi_SL);
  
  List PlantsInst = List::create(
    _["E"]=Einst, _["Ag"]=Aginst, _["An"]=Aninst,
    _["dEdP"] = dEdPInst,
    _["RootPsi"] = RootPsiInst, 
    _["StemPsi"] = StemPsiInst,
    _["LeafPsi"] = LeafPsiInst,
    _["StemSympPsi"] = StemSympPsiInst,
    _["LeafSympPsi"] = LeafSympPsiInst,
    _["StemPLC"] = PLC, 
    _["StemRWC"] = StemRWCInst,
    _["LeafRWC"] = LeafRWCInst,
    _["StemSympRWC"] = StemSympRWCInst,
    _["LeafSympRWC"] = LeafSympRWCInst,
    _["PWB"] = PWBinst);
  DataFrame Plants = DataFrame::create(_["LAI"] = LAIcohort,
                                       _["LAIlive"] = clone(LAIlive),
                                       _["Extraction"] = SoilExtractCoh,
                                       _["Transpiration"] = Eplant,
                                       _["GrossPhotosynthesis"] = Agplant,
                                       _["NetPhotosynthesis"] = Anplant,
                                       _["RootPsi"] = minRootPsi, 
                                       _["StemPsi"] = minStemPsi, 
                                       _["StemPLC"] = PLCm, //Average daily stem PLC
                                       _["LeafPsiMin"] = minLeafPsi, 
                                       _["LeafPsiMax"] = maxLeafPsi, 
                                       _["dEdP"] = dEdPm,//Average daily soilplant conductance
                                       _["DDS"] = DDS, //Daily drought stress is the ratio of average soil plant conductance over its maximum value
                                       _["StemRWC"] = RWCsm,
                                       _["LeafRWC"] = RWClm,
                                       _["LFMC"] = LFMC,
                                       _["WaterBalance"] = PWB);
  Plants.attr("row.names") = above.attr("row.names");
  NumericVector Stand = NumericVector::create(_["LAI"] = LAIcell, 
                                              _["LAIlive"] = LAIcelllive, 
                                              _["LAIexpanded"] = LAIcellexpanded, 
                                              _["LAIdead"] = LAIcelldead);
  
  List l;
  if(!IntegerVector::is_na(stepFunctions)){
    l = List::create(_["cohorts"] = clone(cohorts),
                     _["EnergyBalance"] = EB,
                     _["Stand"] = Stand,
                     _["Extraction"] = SoilWaterExtract,
                     _["ExtractionInst"] = soilLayerExtractInst,
                     _["RhizoPsi"] = minPsiRhizo,
                     _["Plants"] = Plants,
                     _["SunlitLeaves"] = Sunlit,
                     _["ShadeLeaves"] = Shade,
                     _["PlantsInst"] = PlantsInst,
                     _["SunlitLeavesInst"] = SunlitInst,
                     _["ShadeLeavesInst"] = ShadeInst,
                     _["LightExtinction"] = lightExtinctionAbsortion,
                     _["LWRExtinction"] = lwrExtinctionList,
                     _["CanopyTurbulence"] = canopyTurbulence,
                     _["SupplyFunctions"] = supply,
                     _["PhotoSunlitFunctions"] = outPhotoSunlit,
                     _["PhotoShadeFunctions"] = outPhotoShade,
                     _["PMSunlitFunctions"] = outPMSunlit,
                     _["PMShadeFunctions"] = outPMShade);
  } else {
    l = List::create(_["cohorts"] = clone(cohorts),
                     _["EnergyBalance"] = EB,
                     _["Extraction"] = SoilWaterExtract,
                     _["RhizoPsi"] = minPsiRhizo,
                     _["Stand"] = Stand,
                     _["Plants"] = Plants,
                     _["SunlitLeaves"] = Sunlit,
                     _["ShadeLeaves"] = Shade,
                     _["ExtractionInst"] = soilLayerExtractInst,
                     _["PlantsInst"] = PlantsInst,
                     _["SunlitLeavesInst"] = SunlitInst,
                     _["ShadeLeavesInst"] = ShadeInst,
                     _["LightExtinction"] = lightExtinctionAbsortion,
                     _["LWRExtinction"] = lwrExtinctionList,
                     _["CanopyTurbulence"] = canopyTurbulence,
                     _["SupplyFunctions"] = supply);
    
  } 
  l.attr("class") = CharacterVector::create("pwb_day","list");
  return(l);
}

//' @rdname transp_modes
//' 
//' @param canopyEvaporation Canopy evaporation (from interception) for \code{day} (mm).
//' @param soilEvaporation Bare soil evaporation for \code{day} (mm).
//' @param snowMelt Snow melt values  for \code{day} (mm).
//' @param stepFunctions An integer to indicate a simulation step for which photosynthesis and profit maximization functions are desired.
//' 
// [[Rcpp::export("transp_transpirationSperry")]]
List transpirationSperry(List x, DataFrame meteo, int day,
                        double latitude, double elevation, double slope, double aspect,
                        double canopyEvaporation = 0.0, double snowMelt = 0.0, double soilEvaporation = 0.0,
                        int stepFunctions = NA_INTEGER, 
                        bool modifyInput = true) {
  List control = x["control"];
  if(!meteo.containsElementNamed("MinTemperature")) stop("Please include variable 'MinTemperature' in weather input.");
  NumericVector MinTemperature = meteo["MinTemperature"];
  if(!meteo.containsElementNamed("MaxTemperature")) stop("Please include variable 'MaxTemperature' in weather input.");
  NumericVector MaxTemperature = meteo["MaxTemperature"];
  if(!meteo.containsElementNamed("MinRelativeHumidity")) stop("Please include variable 'MinRelativeHumidity' in weather input.");
  NumericVector MinRelativeHumidity = meteo["MinRelativeHumidity"];
  if(!meteo.containsElementNamed("MaxRelativeHumidity")) stop("Please include variable 'MaxRelativeHumidity' in weather input.");
  NumericVector MaxRelativeHumidity = meteo["MaxRelativeHumidity"];
  if(!meteo.containsElementNamed("Radiation")) stop("Please include variable 'Radiation' in weather input.");
  NumericVector Radiation = meteo["Radiation"];
  if(!meteo.containsElementNamed("Precipitation")) stop("Please include variable 'Precipitation' in weather input.");
  NumericVector Precipitation = meteo["Precipitation"];
  NumericVector WindSpeed(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("WindSpeed")) WindSpeed = meteo["WindSpeed"];
  NumericVector CO2(Precipitation.length(), NA_REAL);
  if(meteo.containsElementNamed("CO2")) CO2 = meteo["CO2"];
  
  CharacterVector dateStrings = meteo.attr("row.names");
  std::string c = as<std::string>(dateStrings[day-1]);
  int J = meteoland::radiation_julianDay(std::atoi(c.substr(0, 4).c_str()),std::atoi(c.substr(5,2).c_str()),std::atoi(c.substr(8,2).c_str()));
  double delta = meteoland::radiation_solarDeclination(J);
  double solarConstant = meteoland::radiation_solarConstant(J);
  
  double prec = Precipitation[day-1];
  double rad = Radiation[day-1];
  double tmax = MaxTemperature[day-1];
  double tmin = MinTemperature[day-1];
  double tmaxPrev = tmax;
  double tminPrev = tmin;
  double tminNext = tmin;
  if(day>1) {
    tmaxPrev = MaxTemperature[day-2];
    tminPrev = MinTemperature[day-2];
  }
  if(day<(MaxTemperature.length()-1)) tminNext = MinTemperature[day];
  double rhmax = MaxRelativeHumidity[day-1];
  double rhmin = MinRelativeHumidity[day-1];
  double wind = WindSpeed[day-1];
  double Catm = CO2[day-1];
  if(NumericVector::is_na(Catm)) Catm = control["defaultCO2"];

  NumericVector meteovec = NumericVector::create(
    Named("tmin") = tmin, 
    Named("tmax") = tmax,
    Named("tminPrev") = tminPrev, 
    Named("tmaxPrev") = tmaxPrev, 
    Named("tminNext") = tminNext, 
    Named("prec") = prec,
    Named("rhmin") = rhmin, 
    Named("rhmax") = rhmax, 
    Named("rad") = rad, 
    Named("wind") = wind, 
    Named("Catm") = Catm);
  return(transpirationSperry(x, meteovec,
                     latitude, elevation, slope, aspect,
                     solarConstant, delta,
                     canopyEvaporation, snowMelt, soilEvaporation,
                     false, stepFunctions, 
                     modifyInput));
} 




