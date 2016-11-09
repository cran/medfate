\encoding{UTF-8}
\name{defaultSoilParams}
\alias{defaultSoilParams}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Default soil parameters}
\description{
Creates a list with default soil parameters for model functions.
}
\usage{
defaultSoilParams() 
}
\details{
The function returns a list with default parameters. Users can change those that need to be set to other values and use the list as input for function \code{\link{soil}}.
}
\value{
A list with the following elements (and default values):
\itemize{
\item{\code{SoilDepth (= 1000)}: Soil depth (accounting for the sum of topsoil and subsoil).}
\item{\code{RockLayerDepth (= 4000)}: Depth of the rock layer (i.e. total depth).}
\item{\code{TS_clay (= 25)}: Topsoil clay percentage.}
\item{\code{TS_silt (= 50)}: Topsoil silt percentage.}
\item{\code{TS_sand (= 25)}: Topsoil sand percentage.}
\item{\code{TS_macro (= 0.1)}: Topsoil macroporosity.}
\item{\code{TS_rfc (= 20)}: Topsoil percentage of rock fragment content.}
\item{\code{SS_clay (= 25)}: Subsoil clay percentage.}
\item{\code{SS_silt (= 50)}: Subsoil silt percentage.}
\item{\code{SS_sand (= 25)}: Subsoil sand percentage.}
\item{\code{SS_macro (= 0.1)}: Subsoil macroporosity.}
\item{\code{SS_rfc (= 40)}: Subsoil percentage of rock fragment content.}
\item{\code{RL_clay (= 25)}: Rock layer clay percentage.}
\item{\code{RL_silt (= 50)}: Rock layer silt percentage.}
\item{\code{RL_sand (= 25)}: Rock layer sand percentage.}
\item{\code{RL_macro (= 0.1)}: Rock layer macroporosity.}
\item{\code{RL_rfc (= 85)}: Rock layer percentage of rock fragment content.}
\item{\code{Gsoil (= 0.5)}: Gamma parameter for infiltration.}
\item{\code{Ksoil (= 0.05)}: Kappa parameter for infiltration.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{soil}}, \code{\link{defaultControl}}, \code{\link{SpParamsMED}}
}
\examples{
defaultSoilParams()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.