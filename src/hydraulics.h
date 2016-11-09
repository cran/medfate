#include <Rcpp.h>

#ifndef HYDRAULICS_H
#define HYDRAULICS_H
#endif
using namespace Rcpp;

NumericVector regulatedPsiTwoElements(double Emax, double psiSoil, double ksmax, double kmax, double n, double alpha, double c, double d, double dE = 0.1, double psiMax = -10.0);
