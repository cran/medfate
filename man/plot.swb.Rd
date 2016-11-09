\encoding{UTF-8}
\name{plot.swb}
\alias{plot.swb}
\alias{summary.swb}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Displays soil water balance results}
\description{
Function \code{plot} plots the results of the soil water balance model (see \code{\link{swb}}), whereas function \code{summary} summarizes the model's output in different temporal steps (i.e. weekly, annual, ...).
}
\usage{
\method{plot}{swb}(x, yearAxis=FALSE, type="PET_Precipitation",
    xlim = NULL, ylim=NULL, xlab=NULL, ylab=NULL,...)
\method{summary}{swb}(object, freq="years", output="DailyBalance", FUN=sum, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, object}{An object of class \code{swb}}
  \item{yearAxis}{A boolean to indicate whether the units of the x-axis are years (by default they are days).}
  \item{type}{The information to be plotted:
  \itemize{
    \item{\code{"PET_NetPrec"}:}{ Potential evapotranspiration and Net Precipitation.}
    \item{\code{"Export"}:}{ Water exported through deep drainage and surface runoff.}
    \item{\code{"ET"}:}{ Plant eranspiration and soil evaporation.}
    \item{\code{"Psi"}:}{ Soil water potential.}
    \item{\code{"Theta"}:}{ Soil relative water content.}
    \item{\code{"Vol"}:}{ Soil water volumetric content.}
    \item{\code{"PlantStress"}:}{ Plant cohort average daily drought stress.}
  }
  }
  \item{xlim}{Range of values for x.}
  \item{ylim}{Range of values for y.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{freq}{Frequency of summary statistics (see \code{\link{cut.Date}}).}
  \item{output}{The data table to be summarized. Accepted values are \code{"DailyBalance"}, \code{"PlantStress"}, \code{"PlantTranspiration"} and \code{"SoilWaterBalance"}.}
  \item{FUN}{The function to summarize results (e.g., \code{sum}, \code{mean}, ...)}
  \item{...}{Additional parameters for functions \code{plot} or \code{summary}.}
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (submitted) Coupling a water balance model with forest inventory data to evaluate plant drought stress at the regional level. Agricultural and Forest Meteorology.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swb}}, \code{\link{swbpoints}}
}
\examples{
#Load example daily meteorological data
data(examplemeteo)

#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params
examplesoil = soil(defaultSoilParams())

#Prepare input
x = forest2swbInput(exampleforest,SpParamsMED, examplesoil$dVec)

#Call simulation function
S<-swb(x, examplesoil, examplemeteo)
        
#Plot results
plot(S)

#Monthly summary (averages) of soil water balance
summary(S, freq="months",FUN=mean, output="SoilWaterBalance")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.