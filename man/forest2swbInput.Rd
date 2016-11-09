\encoding{UTF-8}
\name{forest2swbInput}
\alias{swbInput}
\alias{forest2swbInput}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Input variables for soil water balance
}
\description{
Function \code{forest2swbInput} takes an object of class \code{\link{forest}} to calculate input variables for  \code{\link{swb}}. Function \code{swbInput} does the same from input vectors.
}
\usage{
forest2swbInput(x, SpParams, d, gdd = NA, petMode = "Input", hydraulicMode = "Simple")
swbInput(SP, LAI, H, CR,  V, SpParams,  petMode = "Input", hydraulicMode = "Simple")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{\link{forest}}.}
  \item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}} and \code{\link{SpParamsMED}}).}
  \item{d}{A vector containing the depth of each soil layer.}
  \item{gdd}{Growth degree days to account for leaf phenology effects (in Celsius).}
  \item{petMode}{Potential evapotranspiration mode (either 'Input' or 'PenmanMonteith').}
  \item{hydraulicMode}{Hydraulic model (either 'Simple' or 'Sperry').}
  \item{SP}{An integer vector of plant cohort species identity (first species is 0).}
  \item{LAI}{A numeric vector (same length as \code{SP}) with plant cohort leaf area index (one-side leaf area relative to plot area).}
  \item{H}{A numeric vector (same length as \code{SP}) with plant cohort total heights (in cm).}
  \item{CR}{A numeric vector (same length as \code{SP}) with plant cohort crown ratio values (between 0 and 1).}
  \item{V}{A numeric matrix (with three columns and as many rows as the length as \code{SP}) containing the proportion of roots of each plant in each soil layer.}
  
}
\details{
Function \code{forest2swbInput} extracts height and species identity from plant cohorts of \code{x}. For each plant cohort also calculates leaf area index and the distribution of fine roots across soil. Both \code{forest2swbInput} and \code{swbInput} find parameter values for each plant cohort according to the parameters of its species as specified in \code{SpParams}.
}
\value{
A data frame with the following columns (additional columns will difer depending on the value of \code{hydraulicMode} and \code{petMode}):
\itemize{
  \item{\code{LAI}: Leaf area index (in m2/m2).}
  \item{\code{SP}: Species identity (an integer).}
  \item{\code{H}: Plant height (in cm).}
  \item{\code{V.1}, \code{V.2}, \code{V.3}: The proportion of fine roots in each layer.}
  \item{\code{RC_min}: Minimum canopy (stomatal) resistance (in s/m) (for \code{petMode = "PenmanMonteith"}).}
  \item{\code{psiExtr}: Water potential corresponding to 50\% relative conductance (in kPa) (for \code{hydraulicMode = "Simple"}).}
  \item{\code{VC_kxmax}: Maximum hydraulic conductance (in L/m2/day) (for \code{hydraulicMode = "Sperry"}).}
  \item{\code{VC_c}, \code{VC_d}: Parameters of the vulnerability curve (for \code{hydraulicMode = "Sperry"}).}
  \item{\code{k}: PAR extinction coefficient.}
  \item{\code{g}: Canopy water retention capacity per LAI unit (mm/LAI).}
  \item{\code{CR}: Crown ratio (crown length to total height).}
  \item{\code{Sgdd}: Growth degree days needed for leaf budburst (for winter decideous species).}
  \item{\code{Transpiration}: Plant cohort transpiration (filled with zeroes before simulations).}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{\code{\link{swb}}, \code{\link{soil}},  \code{\link{forest}}, \code{\link{SpParamsMED}}, \code{\link{defaultSoilParams}}}

\examples{
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params
s = soil(defaultSoilParams())

#Prepare input
swbInput = forest2swbInput(exampleforest,SpParamsMED, s$dVec)

}
