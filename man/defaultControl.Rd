\encoding{UTF-8}
\name{defaultControl}
\alias{defaultControl}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Default control parameters for models}
\description{
Creates a list with global default parameters for simulation models.
}
\usage{
defaultControl() 
}
\details{
The function returns a list with default parameters. Users can change those defaults that need to be set to other values and use the list as input for model functions. The relevant parameters are different for each model function.
}
\value{
A list, with the following options:
\itemize{
\item{\code{verbose (=FALSE)}: Boolean flag to indicate console output during calculations.}
\item{\code{petMode (="Input")}: Potential evapotranspiration mode. See \code{\link{swb}}.}
\item{\code{hydraulicMode (="Input")}: Hydraulic model. See \code{\link{swb}}.}
}
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{defaultSoilParams}}, \code{\link{SpParamsMED}}, \code{\link{swb}}, \code{\link{swbpoints}}
}