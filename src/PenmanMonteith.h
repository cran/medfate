#include <Rcpp.h>

#ifndef PENMANMONTEITH_H
#define PENMANMONTEITH_H
#endif
using namespace Rcpp;

double netRadiation(double latitude, double elevation,  int J, double Tmin, double Tmax, double RHmin, double RHmax, double R_s, double alpha = 0.08);
double PenmanMonteithPET(double rc, double elevation, double Tmin, double Tmax, double RHmin, double RHmax, double Rn, double u = 5);
