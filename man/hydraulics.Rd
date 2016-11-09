\encoding{UTF-8}
\name{hydraulics}
\alias{hydraulics.E2psiXylem}
\alias{hydraulics.E2psiTwoElements}
\alias{hydraulics.E2psiVanGenuchten}
\alias{hydraulics.Ecrit}
\alias{hydraulics.EXylem}
\alias{hydraulics.psiCrit}
\alias{hydraulics.regulatedPsiXylem}
\alias{hydraulics.regulatedPsiTwoElements}
\alias{hydraulics.supplyFunction}
\alias{hydraulics.vanGenuchtenConductance}
\alias{hydraulics.xylemConductance}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Hydraulic model functions
}
\description{
Set of functions used in the implementation of Sperry and Love's (2015) hydraulic model.
}
\usage{
hydraulics.E2psiXylem(E, psiUpstream, kxmax, c, d, psiStep = -0.01, 
                 psiMax = -10.0)
hydraulics.E2psiVanGenuchten(E, psiSoil, ksmax, n, alpha, 
                             psiStep = -0.01, psiMax = -10.0)
hydraulics.E2psiTwoElements(E, psiSoil, ksmax, kxmax, n, alpha, c, d, 
                            psiStep = -0.001, psiMax = -10.0)
hydraulics.Ecrit(psiUpstream, kxmax, c, d)
hydraulics.EXylem(psiPlant, psiUpstream, kxmax, c, d)
hydraulics.supplyFunction(Emax, psiSoil, ksmax, kxmax, n, alpha, c, d, 
                          dE = 0.1, psiMax = -10.0)
hydraulics.psiCrit(c, d)
hydraulics.regulatedPsiXylem(E, psiUpstream, kxmax, c, d, psiStep = -0.01)
hydraulics.regulatedPsiTwoElements(Emax, psiSoil, ksmax, kxmax, n, alpha, 
                                   c, d, dE = 0.1, psiMax = -10.0)
hydraulics.vanGenuchtenConductance(psi, ksmax, n, alpha)
hydraulics.xylemConductance(psi, kxmax, c, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{E}{Flow per surface unit.}
  \item{Emax}{Maximum flow per surface unit.}
  \item{psi}{Water potential (in MPa).}
  \item{psiUpstream}{Water potential upstream (in MPa). In a one-component model corresponds to soil potential. In a two-component model corresponds to the potential inside the roots.}
  \item{psiPlant}{Plant water potential (in MPa).}
  \item{psiSoil}{Soil water potential (in MPa).}
  \item{psiStep}{Water potential precision (in MPa).}
  \item{psiMax}{Minimum (maximum in absolute value) water potential to be considered (in MPa).}
  \item{dE}{Flow precision.}
  \item{kxmax}{Maximum xylem hydraulic conductance (defined as flow per surface unit and per pressure drop).}
  \item{ksmax}{Maximum rhizosphere hydraulic conductance (defined as flow per surface unit and per pressure drop).}
  \item{c, d}{Parameters of the Weibull function (plant vulnerability curve).}
  \item{n, alpha}{Parameters of the Van Genuchten function (rhizosphere vulnerability curve).}
}
\details{
Details of the hydraulic model are given in a vignette.
}
\value{
Values returned for each function are:
\itemize{
\item{\code{hydraulics.E2psiXylem}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the xylem supply function and given an upstream (soil or root) water potential.}
\item{\code{hydraulics.E2psiVanGenuchten}: The root water potential (in MPa) corresponding to the input flow, according to the rhizosphere supply function and given a soil water potential.}
\item{\code{hydraulics.E2psiTwoElements}: The plant (leaf) water potential (in MPa) corresponding to the input flow, according to the rhizosphere and plant supply functions and given an input soil water potential.}
\item{\code{hydraulics.Ecrit}: The critical flow according to the xylem supply function and given an input soil water potential.}
\item{\code{hydraulics.EXylem}: The flow (integral of the vulnerability curve) according to the xylem supply function and given an input drop in water potential (rhizosphere and plant).}
\item{\code{hydraulics.supplyFunction}: A list with different numeric vectors with information of the two-element supply function:
\itemize{
\item{\code{E}: Flow values (supply values).}
\item{\code{FittedE}: Fitted flow values.}
\item{\code{PsiRoot}: Water potential inside the root values.}
\item{\code{PsiPlant}: Plant water potential values.}
\item{\code{dEdP}: Derivatives of the supply function.}
}
}
\item{\code{hydraulics.psiCrit}: Critical water potential (in MPa) according to the xylem vulnerability curve. }
\item{\code{hydraulics.regulatedPsiXylem}: Plant water potential after regulation (one-element loss function) given an input water potential.}
\item{\code{hydraulics.regulatedPsiTwoElements}: Plant water potential after regulation (two-element loss function) given an input soil water potential.}
\item{\code{hydraulics.vanGenuchtenConductance}: Rhizosphere conductance corresponding to an input water potential (soil vulnerability curve).}
\item{\code{hydraulics.xylemConductance}: Xylem conductance (flow rate per pressure drop) corresponding to an input water potential (plant vulnerability curve).}
}
}
\references{
Sperry, J. S., and D. M. Love. 2015. What plant hydraulics can tell us about responses to climate-change droughts. New Phytologist 207:14–27.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{swb}}
}
