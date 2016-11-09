#include <Rcpp.h>

using namespace Rcpp;
using namespace std;



// [[Rcpp::export(".gammds")]]
double gammds ( double x, double p)
  
  //****************************************************************************80
  //
  //  Purpose:
  //
  //    GAMMDS computes the incomplete Gamma integral.
  //
  //  Discussion:
  //
  //    The parameters must be positive.  An infinite series is used.
  //
  //  Licensing:
  //
  //    This code is distributed under the GNU LGPL license. 
  //
  //  Modified:
  //
  //    22 January 2008
  //
  //  Author:
  //
  //    Original FORTRAN77 version by Chi Leung Lau.
  //    C++ version by John Burkardt.
  //
  //  Reference:
  //
  //    Chi Leung Lau,
  //    Algorithm AS 147:
  //    A Simple Series for the Incomplete Gamma Integral,
  //    Applied Statistics,
  //    Volume 29, Number 1, 1980, pages 113-114.
  //
  //  Parameters:
  //
  //    Input, double X, P, the arguments of the incomplete
  //    Gamma integral.  X and P must be greater than 0.
  //
  //    Output, double GAMMDS, the value of the incomplete
  //    Gamma integral.
  //
{
  double a;
  double arg;
  double c;
  double e = 1.0E-09;
  double f;
  double uflo = 1.0E-37;
  double value;
  //
  //  Check the input.
  //
  if ( x <= 0.0 )
  {
    stop("x <= 0.0 in gammds");
    value = 0.0;
    return value;
  }
  
  if ( p <= 0.0 ) 
  {
    stop("p <= 0.0 in gammds");
    value = 0.0;
    return value;
  }
  //
  //  LGAMMA is the natural logarithm of the gamma function.
  //
  arg = p * log ( x ) - lgamma ( p + 1.0 ) - x;
  
  if ( arg < log ( uflo ) )
  {
    // stop("underflow during the computation in gammds");
    value = NA_REAL;
    return value;
  }
  
  f = exp ( arg );
  
  if ( f == 0.0 )
  {
    // stop("underflow during the computation in gammds");
    value = NA_REAL;
    return value;
  }
  
  //
  //  Series begins.
  //
  c = 1.0;
  value = 1.0;
  a = p;
  
  for ( ; ; )
  {
    a = a + 1.0;
    c = c * x / a;
    value = value + c;
    
    if ( c <= e * value )
    {
      break;
    }
  }
  
  value = value * f;
  
  return value;
}

// [[Rcpp::export(".Egamma")]]
double Egamma(double psi, double kxmax, double c, double d) {
  if(psi>0.0) stop("psi has to be negative");
  double h = 1.0/c;
  double z = pow(-psi/d,c);
  double g = exp(lgamma(h))*gammds(z,h); //Upper incomplete gamma, without the normalizing factor
  return(kxmax*(d/c)*g);
}

// [[Rcpp::export("hydraulics.EXylem")]]
double EXylem(double psiPlant, double psiUpstream, double kxmax, double c, double d) {
  if(psiPlant > psiUpstream) return(NA_REAL);
  return(Egamma(psiPlant, kxmax, c, d)-Egamma(psiUpstream, kxmax, c,d));
}

// [[Rcpp::export("hydraulics.xylemConductance")]]
double xylemConductance(double psi, double kxmax, double c, double d) {
  if(psi>0.0) {
    Rcout<< psi<<"\n";
    stop("psi has to be negative"); 
  }
  return(kxmax*exp(-pow(-psi/d,c)));
}

// [[Rcpp::export("hydraulics.psiCrit")]]
double psiCrit(double c, double d) {
  return(-d * pow(-log(0.01), 1.0/c));
}

// [[Rcpp::export("hydraulics.E2psiXylem")]]
double E2psiXylem(double E, double psiUpstream, double kxmax, double c, double d, double psiStep = -0.01, double psiMax = -10.0) {
  if(E<0.0) stop("E has to be positive");
  if(E==0) return(psiUpstream);
  double psi = psiUpstream;
  double Eg = 0.0;
  double psiPrev = psi;
  while(Eg<E) {
    psiPrev = psi;
    psi = psi + psiStep;
    Eg = EXylem(psi,psiUpstream, kxmax, c, d);
    if(psi<psiMax) return(NA_REAL);
    if(NumericVector::is_na(Eg)) return(NA_REAL);
  }
  return(psiPrev);
}

// [[Rcpp::export("hydraulics.Ecrit")]]
double Ecrit(double psiUpstream, double kxmax, double c, double d) {
  return(EXylem(psiCrit(c,d), psiUpstream, kxmax, c, d));
}

// [[Rcpp::export("hydraulics.regulatedPsiXylem")]]
NumericVector regulatedPsiXylem(double E, double psiUpstream, double kxmax, double c, double d, double psiStep = -0.01) {
  //If Ein > Ecrit then set Ein to Ecrit
  double psiUnregulated = E2psiXylem(E, psiUpstream, kxmax, c, d, psiStep);
  double Ec = Ecrit(psiUpstream, kxmax,c,d);
  double Ein = E;
  if(Ein > Ec) {
    Ein = Ec;
    psiUnregulated = psiCrit(c,d);
  }
  double deltaPsiUnregulated = psiUnregulated - psiUpstream;
  double kp = xylemConductance(psiUpstream, kxmax, c, d);
  double deltaPsiRegulated = deltaPsiUnregulated*(xylemConductance(psiUnregulated, kxmax, c, d)/kp);
  //replace by maximum if found for lower psi values
  // Rcout <<"Initial "<<psiUnregulated << " "<< deltaPsiRegulated <<"\n";
  for(double psi = psiUpstream; psi > psiUnregulated; psi +=psiStep) {
    double deltaPsi = (psi-psiUpstream)*(xylemConductance(psi, kxmax, c,d)/kp);
    // Rcout <<psi << " "<< deltaPsi<< " "<< deltaPsiRegulated <<"\n";
    if(NumericVector::is_na(deltaPsiRegulated)) deltaPsiRegulated = deltaPsi;
    else if(deltaPsi < deltaPsiRegulated) deltaPsiRegulated = deltaPsi;
  }
  double psiRegulated = psiUpstream + deltaPsiRegulated;
  double Efin = EXylem(psiRegulated, psiUpstream, kxmax, c, d);
  double relativeConductance1 = Efin/Ein;
  double relativeConductance2 = Efin/E;
  return(NumericVector::create(psiUnregulated, psiRegulated, Ein, Efin, relativeConductance1, relativeConductance2));
}
  
// [[Rcpp::export("hydraulics.vanGenuchtenConductance")]]
double vanGenuchtenConductance(double psi, double ksmax, double n, double alpha) {
  double v = 1.0/(pow(alpha*abs(psi),n)+1.0);
  return(ksmax*pow(v,(n-1.0)/(2.0*n))*pow(pow((1.0-v),(n-1.0)/n)-1.0,2.0));
//   double ah = alpha*abs(psi);
//   double m = 1.0 - (1.0/n);
//   double den = pow(1.0+pow(ah, n),m/2.0);
//   double num = pow(1.0-pow(ah,n-1.0)*pow(1.0+pow(ah, n),-m),2.0);
//   return(ksmax*(num/den));
}
// [[Rcpp::export("hydraulics.E2psiVanGenuchten")]]
double E2psiVanGenuchten(double E, double psiSoil, double ksmax, double n, double alpha, double psiStep = -0.01, double psiMax = -10.0) {
  if(E<0.0) stop("E has to be positive");
  if(E==0) return(psiSoil);
  double psi = psiSoil;
  double psiPrev = psi;
  double vgPrev = vanGenuchtenConductance(psi, ksmax, n, alpha);
  double vg = vgPrev;
  double Eg = 0.0;
  while(Eg<E) {
    psiPrev = psi;
    vgPrev = vg;
    psi = psi + psiStep;
    vg = vanGenuchtenConductance(psi, ksmax, n, alpha);
    Eg = Eg + ((vg+vgPrev)/2.0)*abs(psiStep);
    if(psi<psiMax) return(NA_REAL);
  }
  return(psiPrev);
}

// [[Rcpp::export("hydraulics.E2psiTwoElements")]]
double E2psiTwoElements(double E, double psiSoil, double ksmax, double kxmax, double n, double alpha, double c, double d, double psiStep = -0.001, double psiMax = -10.0) {
  if(E<0.0) stop("E has to be positive");
  if(E==0) return(psiSoil);
  double psiRoot = E2psiVanGenuchten(E, psiSoil, ksmax, n, alpha, psiStep, psiMax);
  if(NumericVector::is_na(psiRoot)) return(NA_REAL);
  return(E2psiXylem(E, psiRoot, kxmax, c, d, psiStep, psiMax));
}

// [[Rcpp::export("hydraulics.supplyFunction")]]
List supplyFunction(double Emax, double psiSoil, double ksmax, double kxmax, double n, double alpha, double c, double d, double dE = 0.1, double psiMax = -10.0) {
  dE = std::min(dE,Emax/5.0);
  int maxNsteps = round(Emax/dE)+1;
  NumericVector supplyE(maxNsteps);
  NumericVector supplydEdp(maxNsteps);
  NumericVector supplyFittedE(maxNsteps);
  NumericVector supplyPsiRoot(maxNsteps);
  NumericVector supplyPsiPlant(maxNsteps);
  double Eg1 = 0.0;
  double Eg2 = 0.0;
  double psiStep1 = -0.1;
  double psiStep2 = -0.1;
  double psiRoot = psiSoil;
  double psiPlant = psiSoil;
  double vgPrev = vanGenuchtenConductance(psiRoot, ksmax, n, alpha);
  double vg = 0.0;
  double wPrev = xylemConductance(psiPlant, kxmax, c, d);
  double w = 0.0;
  double incr = 0.0;
  supplyPsiRoot[0] = psiSoil;
  supplyPsiPlant[0] = psiSoil;
  supplyE[0] = 0.0;
  double psiPrec = -0.000001;
  for(int i=1;i<maxNsteps;i++) {
    supplyE[i] = supplyE[i-1]+dE;
    psiStep1 = -0.01;
    psiRoot = supplyPsiRoot[i-1];
    vgPrev = vanGenuchtenConductance(psiRoot, ksmax, n, alpha);
    while((psiStep1<psiPrec) & (psiRoot>psiMax))  {
      vg = vanGenuchtenConductance(psiRoot+psiStep1, ksmax, n, alpha);
      incr = ((vg+vgPrev)/2.0)*abs(psiStep1);
      if((Eg1+incr)>supplyE[i]) {
        psiStep1 = psiStep1*0.5;
      } else {
        psiRoot = psiRoot + psiStep1;
        Eg1 = Eg1+incr;
        vgPrev = vg;
      }
    }
    supplyPsiRoot[i] = psiRoot;
    if(supplyPsiRoot[i]<psiMax) supplyPsiRoot[i] = psiMax;
    
    psiStep2 = -0.01;
    Eg2 = 0.0;
    psiPlant = psiRoot;
    wPrev = xylemConductance(psiPlant, kxmax, c, d);
    while((psiStep2<psiPrec) & (psiPlant>psiMax))  {
      w = xylemConductance(psiPlant+psiStep2, kxmax, c, d);
      incr = ((w+wPrev)/2.0)*abs(psiStep2);
      if((Eg2+incr)>supplyE[i]) {
        psiStep2 = psiStep2*0.5;
      } else {
        psiPlant = psiPlant + psiStep2;
        Eg2 = Eg2+incr;
        wPrev = w;
      }
    }
    supplyPsiPlant[i] = psiPlant;
    if(supplyPsiPlant[i]<psiMax) supplyPsiPlant[i] = psiMax;
    supplyFittedE[i] = std::max(supplyFittedE[i-1], Eg2); //Ensure non-decreasing function
    // supplyFittedE[i] = Eg2;
    if(i==1) {
      supplydEdp[0] = (supplyFittedE[1]-supplyFittedE[0])/(supplyPsiPlant[0]-supplyPsiPlant[1]); 
      if((supplyPsiPlant[0]-supplyPsiPlant[1])==0.0) supplydEdp[0] = 0.0;
    }
    else if(i>1) {
      supplydEdp[i-1] = 0.5*(supplyFittedE[i-1]-supplyFittedE[i-2])/(supplyPsiPlant[i-2]-supplyPsiPlant[i-1])+0.5*(supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiPlant[i-1]-supplyPsiPlant[i]); 
      if((supplyPsiPlant[i-2]-supplyPsiPlant[i-1])==0.0) supplydEdp[i-1] = 0.0;
      else if((supplyPsiPlant[i-1]-supplyPsiPlant[i])==0.0) supplydEdp[i-1] = 0.0;
    }
    if(i==(maxNsteps-1)) {
      supplydEdp[i] = (supplyFittedE[i]-supplyFittedE[i-1])/(supplyPsiPlant[i-1]-supplyPsiPlant[i]); 
      if((supplyPsiPlant[i-1]-supplyPsiPlant[i])==0.0) supplydEdp[i] = 0.0;
    }
  }
  return(List::create(Named("E") = supplyE,
                      Named("FittedE") = supplyFittedE,
                      Named("PsiRoot")=supplyPsiRoot, 
                      Named("PsiPlant")=supplyPsiPlant,
                      Named("dEdP")=supplydEdp));
}

// [[Rcpp::export("hydraulics.regulatedPsiTwoElements")]]
NumericVector regulatedPsiTwoElements(double Emax, double psiSoil, double ksmax, double kxmax, double n, double alpha, double c, double d, double dE = 0.1, double psiMax = -10.0) {
  List s = supplyFunction(Emax, psiSoil, ksmax, kxmax, n, alpha, c, d, dE,psiMax);
  NumericVector supplyPsi = s["PsiPlant"];
  NumericVector Efitted = s["FittedE"];
  NumericVector dEdP = s["dEdP"];
  int maxNsteps = Efitted.size();
  double deltaPsiRegulated = 0.0;
  double deltaPsiRegulatedi=0.0;
  double dEdP0 = dEdP[0];
  for(int i=1;i<maxNsteps;i++) {
    if(supplyPsi[i]>psiMax) {
      deltaPsiRegulatedi = (supplyPsi[i] - psiSoil)*std::min(1.0, dEdP[i]/dEdP0);
      // Rcout<<supplydEdp <<" "<<deltaPsiRegulatedi<<"\n";
      if(deltaPsiRegulatedi < deltaPsiRegulated) {
        deltaPsiRegulated = deltaPsiRegulatedi;
      }
    }
  }
  //Regulated potential
  double psiRegulated = psiSoil + deltaPsiRegulated;
  //Find transpiration corresponding to regulated potential
  double ERegulated = 0.0, dEdPRegulated = 0.0;
  for(int i=1;i<maxNsteps;i++) {
    if((supplyPsi[i-1] >= psiRegulated) & (supplyPsi[i]<psiRegulated)) {
      ERegulated = Efitted[i]*abs((supplyPsi[i-1]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i])) + Efitted[i-1]*abs((supplyPsi[i]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i]));
      ERegulated = std::min(ERegulated, Emax);
      psiRegulated = supplyPsi[i]*abs((supplyPsi[i-1]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i])) + supplyPsi[i-1]*abs((supplyPsi[i]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i]));
      dEdPRegulated = dEdP[i]*abs((supplyPsi[i-1]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i])) + dEdP[i-1]*abs((supplyPsi[i]-psiRegulated)/(supplyPsi[i-1]-supplyPsi[i]));
      if((supplyPsi[i-1]-supplyPsi[i])==0.0) dEdPRegulated = 0.0;
      // Rcout<<dEdP[i]<< " "<<dEdP[i-1]<< " "<<dEdPRegulated<<"\n";
      break;
    }
  }
  return(NumericVector::create(supplyPsi[maxNsteps-1], psiRegulated, Efitted[maxNsteps-1], ERegulated, dEdPRegulated));
}
