% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{species_values}
\alias{species_values}
\alias{species_basalArea}
\alias{species_cover}
\alias{species_density}
\alias{species_foliarBiomass}
\alias{species_fuelLoading}
\alias{species_LAI}
\alias{species_characterParameter}
\alias{species_parameter}
\title{Species description functions}
\usage{
species_basalArea(x, SpParams)

species_cover(x, SpParams)

species_density(x, SpParams)

species_foliarBiomass(x, SpParams, gdd = NA_real_)

species_fuelLoading(x, SpParams, gdd = NA_real_, includeDead = TRUE)

species_LAI(x, SpParams, gdd = NA_real_, bounded = TRUE)

species_characterParameter(species, SpParams, parName)

species_parameter(
  species,
  SpParams,
  parName,
  fillMissing = TRUE,
  fillWithGenus = TRUE
)
}
\arguments{
\item{x}{An object of class \code{\link{forest}}.}

\item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}

\item{gdd}{Growth degree days (to account for leaf phenology effects).}

\item{includeDead}{A flag to indicate that standing dead fuels (dead branches) are included.}

\item{bounded}{A boolean flag to indicate that extreme values should be prevented (maximum tree LAI = 7 and maximum shrub LAI = 3)}

\item{species}{A character vector of species names.}

\item{parName}{A string with a parameter name.}

\item{fillMissing}{A boolean flag to try imputation on missing values.}

\item{fillWithGenus}{A boolean flag to try imputation of missing values using genus values.}
}
\value{
A vector with values for each species in \code{SpParams}:
\itemize{
\item{\code{species_basalArea}: Species basal area (m2/ha).}
\item{\code{species_cover}: Shrub cover (in percent).}
\item{\code{species_density}: Plant density (ind/ha). Tree density is directly taken from the forest object, while the shrub density is estimated from cover and height by calculating the area of a single individual.}
\item{\code{species_foliarBiomass}: Standing biomass of leaves (in kg/m2).}
\item{\code{species_fuel}: Fine fuel load (in kg/m2).}
\item{\code{species_LAI}: Leaf area index (m2/m2).}
\item{\code{species_phytovolume}: Shrub phytovolume (m3/m2).}
\item{\code{species_parameter}: A numeric vector with the parameter values of each input species.}
\item{\code{species_characterParameter}: A character vector with the parameter values of each input species.}
}
}
\description{
Functions to calculate attributes of a \code{\link{forest}} object by species or to extract species parameters from a species parameter table (\code{\link{SpParamsMED}}).
}
\examples{
# Default species parameterization
data(SpParamsMED)

# Load example plot
data(exampleforest)

# Species basal area in the forest plot
species_basalArea(exampleforest, SpParamsMED)
  
# Value of parameter "Psi_Extract" for two species
species_parameter(c("Pinus halepensis", "Quercus ilex"), SpParamsMED, "Psi_Extract")
    
}
\seealso{
\code{\link{spwb}}, \code{\link{forest}}, \code{\link{plant_basalArea}}, \code{\link{summary.forest}}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
\keyword{internal}
