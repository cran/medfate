% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forest.R
\name{forest}
\alias{forest}
\alias{emptyforest}
\alias{summary.forest}
\alias{print.summary.forest}
\title{Forest description}
\usage{
emptyforest(ntree = 0, nshrub = 0, nseed = 0)

\method{summary}{forest}(object, SpParams, ...)

\method{print}{summary.forest}(x, digits = getOption("digits"), ...)
}
\arguments{
\item{ntree, nshrub}{Number of tree and shrub cohorts, respectively.}

\item{nseed}{Number of species in the seed bank.}

\item{object}{An object of class \code{forest} has the following structure (see details):
\itemize{
  \item{\code{treeData}: A data frame of tree cohorts (in rows) and the following columns:
      \itemize{
        \item{\code{Species}: String with species (taxon) name or a non-negative integer for tree species identity (i.e., 0,1,2,...) matching SpParams.}
        \item{\code{Height}: Total tree height (in cm).}
        \item{\code{DBH}: Tree diameter at breast height (in cm).}
        \item{\code{N}: Density (number of individuals/hectare) that the measured tree represents.}
        \item{\code{Z50}: Depth (in mm) corresponding to 50\% of fine roots.}
        \item{\code{Z95}: Depth (in mm) corresponding to 95\% of fine roots.}
     }
  }
  \item{\code{shrubData}: A data frame of shrub cohorts (in rows) and the following columns:
      \itemize{
        \item{\code{Species}: String with species (taxon) name or a non-negative integer for shrub species identity (i.e., 0,1,2,...) matching SpParams.}
        \item{\code{Height}: Average total height of plants (in cm).}
        \item{\code{Cover}: Percent cover.}
        \item{\code{Z50}: Depth (in mm) corresponding to 50\% of fine roots.}
        \item{\code{Z95}: Depth (in mm) corresponding to 95\% of fine roots.}
      }
  }
  \item{\code{herbCover}: Percent cover of the herb layer (optional).}
  \item{\code{herbHeight}: Mean height (in cm) of the herb layer (optional).}
}}

\item{SpParams}{A data frame with species parameters (see \code{\link{SpParamsMED}}).}

\item{...}{Additional parameters for functions \code{\link{summary}} and \code{\link{print}}.}

\item{x}{The object returned by \code{summary.forest}.}

\item{digits}{Minimal number of significant digits.}
}
\value{
Function \code{summary.forest} returns a list with several structural attributes, such as the basal area and LAI of the forest. 
Function \code{emptyforest} returns an empty \code{forest} object.
}
\description{
Description of a forest stand
}
\details{
Function \code{summary.forest} can be used to summarize a \code{forest} object in the console. 
Function \code{emptyforest} creates an empty \code{forest} object.

The structure presented above for \code{forest} objects corresponds to the required data elements. 
A \code{forest} object can contain additional information when this is available. Data frames \code{treeData} 
and \code{shrubData} can contain additional columns:
\itemize{
  \item{\code{LAI}: Leaf area index (m2/m2)}
  \item{\code{FoliarBiomass}: Standing dry biomass of leaves (kg/m2)}
  \item{\code{FuelLoading}: Fine fuel loading (kg/m2)}
  \item{\code{CrownRatio}: The ratio between crown length and total height (between 0 and 1)}
}
Similarly, one can define \code{forest} list elements \code{herbLAI}, \code{herbFoliarBiomass} or \code{herbFuelLoading}.
All these values are used to override allometry-based estimates of those variables when initializing
inputs for functions \code{\link{spwb}} or \code{\link{spwb_day}}. Note that leaf area index, foliar biomass and
fuel loading are related entities, and they are treated as such in medfate. Therefore, users are expected to supply 
one or the other, and not all of them at the same time.
}
\examples{
data(exampleforestMED)
data(SpParamsMED)

# Example forest
exampleforestMED

# Summary of example forest
summary(exampleforestMED, SpParamsMED)

}
\seealso{
\code{\link{exampleforestMED}}, \code{\link{forest_mergeTrees}},  \code{\link{plot.forest}}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, CREAF
}
