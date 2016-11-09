#include <Rcpp.h>

#ifndef SWB_H
#define SWB_H
#endif
using namespace Rcpp;

NumericVector K2Psi(NumericVector K, NumericVector psi_extr, int ws);
NumericVector er(IntegerVector DOY, double ERconv=0.05, double ERsyn = 0.2);
NumericVector gdd(IntegerVector DOY, NumericVector Temp, double Tbase = 5.0);
double soilevaporation(double DEF,double PETs, double Gsoil);
double infiltrationDay(double NetPrec, double Ssoil);

List swbDay1(DataFrame x, List soil, double gdd, double pet, double rain, double er, double runon=0.0, String hydraulicMode = "Simple", bool verbose=false);
