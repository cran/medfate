\encoding{UTF-8}
\name{swb}
\alias{swb}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Soil water balance model}
\description{
Function \code{swb} is a soil water balance model that determines changes in soil moisture, soil water potentials and plant drought stress at daily steps for a given forest stand during a period specified in the input climatic data. 
}
\usage{
swb(x, soil, meteo, petMode = "Input", hydraulicMode = "Simple", 
    latitude = NA, elevation = NA, verbose=TRUE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A data frame with parameters for each plant cohort (see \code{\link{forest2swbInput}}).}
  \item{soil}{A list containing the description of the soil (see \code{\link{soil}}).}
  \item{meteo}{A data frame with daily meteorological data series. When  \code{petMode = "Input"} the following columns are required:
   \itemize{
     \item{\code{DOY}: Day of the year (Julian day).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{MeanTemperature}: Mean temperature (in degrees Celsius).}
     \item{\code{PET}: Potential evapotranspiration (in mm).}
    }
  When \code{petMode = "PenmanMonteith"} the following columns are required:
   \itemize{
     \item{\code{DOY}: Day of the year (Julian day).}
     \item{\code{Precipitation}: Precipitation (in mm).}
     \item{\code{MeanTemperature}: Mean temperature(in degrees Celsius).}
     \item{\code{MinTemperature}: Minimum temperature (in degrees Celsius).}
     \item{\code{MaxTemperature}: Maximum temperature (in degrees Celsius).}
     \item{\code{MinRelativeHumidity}: Minimum relative humidity (in percent).}
     \item{\code{MaxRelativeHumidity}: Maximum relative humidity (in percent).}
     \item{\code{Radiation}: Solar radiation (in MJ/m2/day).}
     \item{\code{WindSpeed}: Wind speed (in m/s). If not available, this column can be left with \code{NA} values.}
    }
  }
  \item{petMode}{Potential evapotranspiration mode (either \code{"Input"} or \code{"PenmanMonteith"}). If \code{petMode = "Input"} the function takes column 'PET' in \code{meteo} as potential evapotranspiration. If \code{petMode = "PenmanMonteith"} the function calculates 'PET' internally, taking into account the current leaf area index of the stand.}
  \item{hydraulicMode}{Hydraulic model (either 'Simple' or 'Sperry').}
  \item{latitude}{Latitude (in radians).}
  \item{elevation}{Elevation above sea level (in m).}
  \item{verbose}{Boolean flag to indicate console output during calculations.}
}
\details{
  Detailed model description is available in the vignettes section. The model using simple hydraulic model is described in De Caceres et al. (2015). WARNING: objects \code{x} and \code{soil} are modified during the simulation.
}
\value{
  A list of class 'swb' with the following elements:
  \itemize{
  \item{\code{"petMode"}: Potential evapotranspiration mode (either \code{"Input"} or \code{"PenmanMonteith"}).}
  \item{\code{"hydraulicMode"}: Hydraulic model (either 'Simple' or 'Sperry').}
  \item{\code{"DailyBalance"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
  \itemize{
  \item{\code{"LAIcell"}: The LAI of the stand (accounting for leaf phenology) (in m2/m2).}
  \item{\code{"Cm"}: The water retention capacity of the stand (in mm) (accounting for leaf phenology).}
  \item{\code{"Lground"}: The proportion of PAR that reaches the ground (accounting for leaf phenology).}
  \item{\code{"PET"}: Potential evapotranspiration (in mm).}
  \item{\code{"Rainfall"}: Input precipitation (in mm).}
  \item{\code{"NetPrec"}: Net precipitation, after accounting for interception (in mm).}
  \item{\code{"Infiltration"}: The amount of water infiltrating into the soil (in mm).}
  \item{\code{"Runoff"}: The amount of water exported via surface runoff (in mm).}
  \item{\code{"DeepDrainage"}: The amount of water exported via deep drainage (in mm).}
  \item{\code{"Etot"}: Evapotranspiration (in mm).}
  \item{\code{"Esoil"}: Bare soil evaporation (in mm).}
  \item{\code{"Eplanttot"}: Plant transpiration (considering all soil layers) (in mm).}
  \item{\code{"Eplant1"}: Plant transpiration from soil layer 1 (in mm).}
  \item{\code{"Eplant2"}: Plant transpiration from soil layer 2 (in mm).}
  \item{\code{"Eplant3"}: Plant transpiration from soil layer 3 (in mm).}
  }
  \item{\code{"SoilWaterBalance"}: A data frame where different variables (in columns) are given for each simulated day (in rows):}
  \itemize{
  \item{\code{"W1"}: Relative soil moisture content (relative to field capacity) in layer 1.}
  \item{\code{"W2"}: Relative soil moisture content (relative to field capacity) in layer 2.}
  \item{\code{"W3"}: Relative soil moisture content (relative to field capacity) in layer 3.}
  \item{\code{"ML1"}: Soil water volume in layer 1 (in L/m2).}
  \item{\code{"ML2"}: Soil water volume in layer 2 (in L/m2).}
  \item{\code{"ML3"}: Soil water volume in layer 3 (in L/m2).}
  \item{\code{"MLTot"}: Total soil water volume (in L/m2).}
  \item{\code{"psi1"}: Soil water potential in layer 1 (in kPa).}
  \item{\code{"psi2"}: Soil water potential in layer 2 (in kPa).}
  \item{\code{"psi3"}: Soil water potential in layer 3 (in kPa).}
  }
  \item{\code{"PlantTranspiration"}: A data frame with the amount of daily transpiration (in mm) for each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantStress"}: A data frame with the amount of daily stress suffered by each plant cohort (relative whole-plant conductance). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantTmax"}: A data frame with the amount of daily maximum transpiration (in mm) for each plant cohort. Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  \item{\code{"PlantPsi"} (only for \code{hydraulicMode = "Sperry"}): A data frame with the average daily water potential each plant (in kPa). Days are in rows and plant cohorts are in columns. Columns in this data frame correspond to the elements in 'SP'.}
  }
}
\references{
De \enc{Cáceres}{Caceres} M, \enc{Martínez}{Martinez}-Vilalta J, Coll L, Llorens P, Casals P, Poyatos R, Pausas JG, Brotons L. (2015) Coupling a water balance model with forest inventory data to predict drought stress: the role of forest structural changes vs. climate changes. Agricultural and Forest Meteorology (doi:10.1016/j.agrformet.2015.06.012).
}
\author{
Miquel De \enc{Cáceres}{Caceres} Ainsa, Centre Tecnologic Forestal de Catalunya
}
\seealso{
\code{\link{plot.swb}},  \code{\link{swbpoints}}, \code{\link{swbgrid}},\code{\link{forest}}
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

#Call simulation function
S1<-swb(forest2swbInput(exampleforest,SpParamsMED, examplesoil$dVec), 
        examplesoil, examplemeteo)
        
#Plot results
plot(S1)

#Monthly summary (averages) of soil water balance
summary(S1, freq="months",FUN=mean, output="SoilWaterBalance")


#Simulation with PET calculated with Penman-Monteith
examplesoil = soil(defaultSoilParams())
S2<-swb(forest2swbInput(exampleforest,SpParamsMED, examplesoil$dVec, petMode="PenmanMonteith"), 
        examplesoil, examplemeteo, latitude = 0.73, 
        elevation = 100, petMode="PenmanMonteith")
        
#Simulation with transpiration calculated using Sperry & Love
examplesoil = soil(defaultSoilParams())
S3<-swb(forest2swbInput(exampleforest,SpParamsMED, examplesoil$dVec, hydraulicMode="Sperry"), 
        examplesoil, examplemeteo, hydraulicMode="Sperry")

#Simulation with PET calculated using PenmanMonteiht
#and transpiration calculated using Sperry & Love
examplesoil = soil(defaultSoilParams())
S4<-swb(forest2swbInput(exampleforest,SpParamsMED, examplesoil$dVec, 
                        petMode = "PenmanMonteith", hydraulicMode="Sperry"), 
        examplesoil, examplemeteo, latitude = 0.73, 
        elevation = 100, petMode = "PenmanMonteith", hydraulicMode="Sperry")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.