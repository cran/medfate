#include <Rcpp.h>

#ifndef FORESTUTILS_H
#define FORESTUTILS_H
#endif
using namespace Rcpp;

void checkSpeciesParameters(DataFrame SpParams, CharacterVector params);

NumericVector leafDevelopmentStatus(NumericVector Sgdd, double gdd);
double leafDevelopmentStatus(double Sgdd, double gdd);

NumericVector ldrRS_one(double Z50, double Z95, NumericVector d);
NumericVector conicRS_one(double Z, double d1, double d2, double d3);
NumericMatrix conicRS(NumericVector Z, NumericVector d);
NumericMatrix ldrRS(NumericVector Z50, NumericVector Z95, NumericVector d);
NumericMatrix ldrProfile(NumericVector Z50, NumericVector Z95, NumericVector d);

NumericVector cohortParameter(List x, DataFrame SpParams, String parName);
CharacterVector cohortCharacterParameter(List x, DataFrame SpParams, String parName);

NumericVector cohortHeight(List x);

NumericVector treeBasalArea(NumericVector N, NumericVector dbh);
NumericVector treeCohortBasalArea(List x);
NumericVector cohortBasalArea(List x);
NumericVector dbhClassBasalArea(List x, NumericVector DBHbreaks);
double forestBasalArea(List x);
double forestBasalAreaForMinDBH(List x, double minDBH);

double treeDensity(List x);
double minDBHDensity(List x, double minDBH);
NumericVector dbhClassDensity(List x, NumericVector DBHbreaks);

NumericVector treeFuel(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true);
NumericVector shrubFuel(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true);
NumericVector cohortFuel(List x, DataFrame SpParams, double gdd = NA_REAL, bool includeDead = true);

NumericVector cohortCrownRatio(List x, DataFrame SpParams);
NumericVector cohortCrownBaseHeight(List x, DataFrame SpParams);
NumericVector cohortCrownLength(List x, DataFrame SpParams);

NumericVector cohortFoliarBiomass(List x, DataFrame SpParams, double gdd = NA_REAL);
NumericVector cohortEquilibriumLeafLitter(List x, DataFrame SpParams, double AET = 800);
NumericVector cohortEquilibriumSmallBranchLitter(List x, DataFrame SpParams, double smallBranchDecompositionRate = 0.81);

NumericVector treeLAI(IntegerVector SP, NumericVector N, NumericVector dbh, DataFrame SpParams, NumericVector pEmb=NumericVector(0), double gdd = NA_REAL);
NumericVector shrubLAI(IntegerVector SP, NumericVector Cover, NumericVector H, DataFrame SpParams, double gdd = NA_REAL);
NumericVector cohortLAI(List x, DataFrame SpParams, double gdd = NA_REAL);
NumericMatrix LAIdistribution(NumericVector z, NumericVector LAI, NumericVector H, NumericVector CR);
NumericMatrix LAIdistribution(NumericVector z, List x, DataFrame SpParams, double gdd = NA_REAL);

double shrubCover(List x, double excludeMinHeight = 0.0);

DataFrame forest2swbInput(List x, DataFrame SpParams, NumericVector d, double gdd = NA_REAL, String petMode = "Input", String hydraulicMode = "Simple");

void deleteTreeCohort(List x, int treeCohort);
void deleteShrubCohort(List x, int shrubCohort);

int minDBHTreeCohort(List x, double excludeMin = 0.0);
