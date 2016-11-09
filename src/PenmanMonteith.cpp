#include <Rcpp.h>
using namespace Rcpp;

const double SIGMA = 4.903*pow(10,-9.0); //Stefan-Boltzmann constant MJ/K^4/m2/day
const double GSC = 0.0820; //solar constant in MJ/m2/min

/*
 * latitude - Latitude (radians)
 * elevation - Elevation (m)
 * Tmax - Maximum temperature (Celsius)
 * Tmin - Minimum temperature (Celsius)
 * RHmax - Maximum relative humidity (%)
 * RHmin - Minimum relative humidity (%)
 * J - Julian day
 * R_s - Incident solar radiation (MJ/m2)
 * alpha - Albedo (from 0 to 1)
 */
double netRadiation(double latitude, double elevation,  int J, 
                    double Tmin, double Tmax, double RHmin, double RHmax, double R_s, 
                    double alpha = 0.08) {
  double d_r2 = 1.0 + 0.033 * cos(2.0 * PI/365.0 * ((double)J));
  double delta2 = 0.409 * sin(2.0 * PI/365.0 * ((double)J) - 1.39);
  double w_s = acos(-tan(latitude) * tan(delta2));
  // double N = 24.0/PI*w_s; (N not used)
  double R_a = (1440.0/PI) * d_r2 * GSC * (w_s * sin(latitude) * sin(delta2) + cos(latitude) * cos(delta2) * sin(w_s));
  double R_so = (0.75 + (2.0 * 0.00001) * elevation) * R_a;
  //Saturation vapour pressure at temperature kPa
  double vs_Tmax = 0.6108 * exp(17.27 * Tmax/(Tmax + 237.3));
  double vs_Tmin = 0.6108 * exp(17.27 * Tmax/(Tmax + 237.3));
  //Mean daily actual vapour pressure, kPa
  double vabar = (vs_Tmin * (RHmax/100.0) + vs_Tmax * (RHmin/100.0))/2.0;
  //Net longwave radiation (MJ.m^-2.day^-1)
  double R_nl = SIGMA * (0.34 - 0.14 * sqrt(vabar)) * (pow(Tmax + 273.2,4.0) + pow(Tmin + 273.2,4.0))/2.0 * (1.35 * R_s/R_so - 0.35);
  //Net shortwave radiation (MJ.m^-2.day^-1, after acounting for surface albedo)
  double R_ns = (1.0 - alpha) * R_s; 
  //Net radiation
  return(std::max(0.0,R_ns - R_nl)); 
}
/*
 *   Calculates daily potential evapotranspiration (mm.m^-2.day^-1) for a surface given microenvironmental conditions
 *   (absorbed radiation, wind, temperature, relative humidity) and given a canopy stomatal resistance.
 *   Basically is the Penman-Monteith equation for daily steps, where aerodynamic resistance to heat flux
 *   is calculated internally (from wind and temperature) and the canopy (stomatal) resistance is an input.
 *   
 *  
 *   REF:  McMahon et al. (2013) 
 *   Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis
 *   Hydrology & Earth System Sciences
 *   See also:  http://www.fao.org/docrep/x0490e/x0490e06.htm
 *   
 *   rc[s/m] - Canopy vapour flux (stomatal) resistance
 *   elevation [m] - Elevation over sea level
 *   Tmin, Tmax [Celsius] - Minimum and maximum temperature
 *   RHmin, RHmax [%] - Minimum and maximum relative humidity
 *   Rn [MJ.m^-2.day^-1] - daily net radiation
 *   u [m.s^-1] - wind speed
 */
// [[Rcpp::export(".PenmanMonteithPET")]]
double PenmanMonteithPET(double rc, double elevation, 
                         double Tmin, double Tmax,
                         double RHmin, double RHmax,
                         double Rn, double u = NA_REAL) {
  double Tday = (Tmax + Tmin)/2.0;
  //Saturation vapour pressure at temperature kPa
  double vs_Tmax = 0.6108 * exp(17.27 * Tmax/(Tmax + 237.3));
  double vs_Tmin = 0.6108 * exp(17.27 * Tmin/(Tmin + 237.3));
  //Daily saturation vapour pressure kPa
  double vas = (vs_Tmax + vs_Tmin)/2.0;
  //Mean daily actual vapour pressure, kPa
  double vabar = (vs_Tmin * RHmax/100.0 + vs_Tmax * RHmin/100.0)/2.0;
  //Atmospheric pressure kPa
  double Pa = 101.32500*pow((293.0-0.065*elevation)/293.0,5.2553);
  //Slope of the saturation vapor pressure curve kPa.Celsius^-1
  double delta = 4098.0 * (0.6108 * exp((17.27 * Tday)/(Tday + 237.3)))/pow(Tday + 237.3,2.0);
  //Vapour pressure deficit, kPa
  double vpd = (vas - vabar);
  //Latent heat of vaporisation MJ.kg^-1
  double lambda = 2.5023 - 0.00243054*Tday;
  // Rcout<<"Lambda: "<< lambda<<"\n";
  //Psychrometric constant kPa.Celsius^-1
  double gamma = (0.00163*Pa)/lambda;
  // Rcout<<"Gamma: "<< gamma<<"\n";
  //Density of air kg.m^-3
  double rho = (Pa/(1.01*(Tday+273.16)*0.287));
  //Unit conversion from s^-1 to day^-1
  double Kt = 86400.0;
  //Specific heat of air MJ.kg^1.Celsius^-1
  double Cp = (gamma*0.622*lambda)/Pa; 

  //Radiative resistance to heat flux s.m^-1
  // double rr = Kt*(rho*Cp)/(4.0*SIGMA*pow(Tday+273.0,3.0));
  //Convective resistance to heat flux s.m^-1
  if(NumericVector::is_na(u)) u = 2.0;
  u = std::max(u, 0.000001);
  double rh = 208.0/u; //Aerodynamic resistance to convective heat transfer
  //Resistance to sensible heat flux (from Biome BGC) s.m^-1
  // double ra = (rh*rr)/(rh+rr);
  rc = std::max(rc, 70.0);
  double ra = rh;
//   double N1fao = (0.408*delta*(Rn));
//   double N2fao = gamma*(900.0/(Tday+273.0))*u*vpd;
//   double Dfao = delta + gamma*(1.0+0.34*u);
//   double Efao = (N1fao+N2fao)/(Dfao);
//   Rcout<<"N1(FAO) : "<<N1fao<< " N2(FAO) : "<<N2fao<<" D(FAO) "<<Dfao<<" E(FAO) "<< Efao <<"\n";
  double N1 = ((delta/lambda)*Rn);
  double N2 = Kt*(vpd*(rho*Cp)/(lambda*ra));
  double D = delta + gamma*(1.0+(rc/ra));
  double E = (N1+ N2)/(D);
  // Rcout<<"N1 : "<<N1<< " N2 : "<<N2<<" D "<<D<<" E "<< E <<"\n";
  return(E);
}
