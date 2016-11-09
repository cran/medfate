#include <Rcpp.h>

#ifndef FUELMOISTURE_H
#define FUELMOISTURE_H
#endif
using namespace Rcpp;

NumericVector sunRiseSet(double L0, double A, double I, double delta);
NumericVector fuelConditions(double airTemp, double airHumidity, double fuelRadiation, double fuelWindSpeed);
double fine1hday(double m0, double fuelTemp, double fuelHumidity, double fuelWindSpeed, double effRain);
double coarse10hday(double m0, 
                    double prevFuelTempMax, double prevFuelHumidityMin,
                    double currFuelTempMin, double currFuelHumidityMax,
                    double rainDuration);
double coarse100hday(double m0, double fuelTempMin, double fuelHumidityMax, 
                     double fuelTempMax, double fuelHumidityMin, 
                     double numSunHours, double rainDuration);
DataFrame deadFuelMoisture(NumericVector m0, NumericVector airTemp, NumericVector airHumidity, NumericVector fuelRadiation, NumericVector fuelWindSpeed, NumericVector effRain, NumericVector rainDuration);

NumericVector cohortFuelMoistureContent(List swbDay, DataFrame swbInput, DataFrame SpParams, int WeibullShape=3);
double canopyLiveFuelMoisture(double canopyBaseHeight, double canopyTopHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR);
double fuelbedLiveFuelMoisture(double fuelbedHeight, NumericVector cohortFMC, NumericVector cohortLoading, NumericVector H, NumericVector CR);
