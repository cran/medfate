\encoding{UTF-8}
\name{root}
\alias{root.ldrRS}
\alias{root.ldrProfile}
\alias{root.conicRS}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Distribution of fine roots}
\description{
Functions to calculate the distribution of fine roots within the soil, given root system parameters and soil layer definition.
}
\usage{
root.conicRS(Z, d)
root.ldrRS(Z50, Z95, d)
root.ldrProfile(Z50, Z95, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z50}{A vector of depths (in mm) corresponding to 50\% of roots.}
  \item{Z95}{A vector of depths (in mm) corresponding to 95\% of roots.}
  \item{Z}{A vector of depths (in mm) corresponding to the root cone tip.}
  \item{d}{The width (in mm) corresponding to each soil layer.}
}
\details{
  Function \code{conicRS} assumes a conic distribution of fine roots, whereas function \code{ldrRS} uses the linear dose response model of Schenck & Jackson (2002).
}
\value{
  Functions \code{conicRS} and \code{ldrRS} return a matrix with as many rows as elements in \code{Z} and three columns (corresponding to soil layers). Function \code{ldrProfile} returns a matrix with as many rows as elements in \code{Z90} and as many columns as elements in \code{d}. Values in all cases correspond to the proportion of fine roots in each soil layer.
}
\references{
Schenk, H., Jackson, R., 2002. The global biogeography of roots. Ecol. Monogr. 72, 311–328.
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
 \code{\link{swb}}, \code{\link{forest2swbInput}}, \code{\link{soil}}
}
\examples{
#Load example plot plant data
data(exampleforest)

#Default species parameterization
data(SpParamsMED)

#Initialize soil with default soil params
S = soil(defaultSoilParams())

#Calculate conic root system for trees
V1 = root.conicRS(Z=rep(2000,nrow(exampleforest$treeData)), S$dVec)            
print(V1)

#Calculate LDR root system for trees (Schenck & Jackson 2002)
V2 = root.ldrRS(Z50 = rep(200,nrow(exampleforest$treeData)), 
          Z95 = rep(1000,nrow(exampleforest$treeData)), S$dVec)
print(V2)     

#Plot LDR root systems (31 layers of 5 mm each)
P = root.ldrProfile(c(100,500), c(200,1500), rep(50, 31))
matplot(x=t(P), y=seq(0,1500, by=50), type="l", xlim=c(0,0.5), ylim=c(1500,0), ylab = "Depth (mm)", 
        xlab="Root density")
legend("bottomright", legend=c("Z50 = 100, Z95 = 200", 
       "Z50 = 200, Z95 = 1500"), col=c("black","red"), lty=1:2, bty="n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.