\name{CNModel}
\alias{CNModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{"CNModel"
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CNModel(fnc_CNModel, CNavg = 75, IaFrac = 0.05, fnc_slope = 0, fnc_aspect = 0, func_DAWC = 0.3, func_z = 1000, fnc_fcres = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnc_CNModel}{
%%     ~~Describe \code{fnc_CNModel} here~~
}
  \item{CNavg}{
%%     ~~Describe \code{CNavg} here~~
}
  \item{IaFrac}{
%%     ~~Describe \code{IaFrac} here~~
}
  \item{fnc_slope}{
%%     ~~Describe \code{fnc_slope} here~~
}
  \item{fnc_aspect}{
%%     ~~Describe \code{fnc_aspect} here~~
}
  \item{func_DAWC}{
%%     ~~Describe \code{func_DAWC} here~~
}
  \item{func_z}{
%%     ~~Describe \code{func_z} here~~
}
  \item{fnc_fcres}{
%%     ~~Describe \code{fnc_fcres} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (fnc_CNModel, CNavg = 75, IaFrac = 0.05, fnc_slope = 0, 
    fnc_aspect = 0, func_DAWC = 0.3, func_z = 1000, fnc_fcres = 0.3) 
{
    attach(fnc_CNModel)
    SNO_Energy = SnowMelt(date, P, MaxTemp - 3, MinTemp - 3, 
        myflowgage$declat, slope = fnc_slope, aspect = fnc_aspect, 
        tempHt = 1, windHt = 2, groundAlbedo = 0.25, SurfEmissiv = 0.95, 
        windSp = 2, forest = 0, startingSnowDepth_m = 0, startingSnowDensity_kg_m3 = 450)
    detach(fnc_CNModel)
    fnc_CNModel$SNO = SNO_Energy$SnowWaterEq_mm
    fnc_CNModel$SNOmlt = SNO_Energy$SnowMelt_mm
    fnc_CNModel$SnowfallWatEq_mm = SNO_Energy$SnowfallWatEq_mm
    fnc_CNModel$SnowMelt_mm = SNO_Energy$SnowMelt_mm
    attach(fnc_CNModel)
    fnc_CNModel$Albedo = 0.23
    fnc_CNModel$Albedo[fnc_CNModel$SNO > 0] = 0.95
    PET = PET_fromTemp(Jday = (1 + as.POSIXlt(date)$yday), Tmax_C = MaxTemp, 
        Tmin_C = MinTemp, lat_radians = myflowgage$declat * pi/180) * 
        1000
    fnc_CNModel$PET = PET
    detach(fnc_CNModel)
    rm(list = "PET")
    fnc_CNModel$AWC = func_DAWC * func_z
    fnc_CNModel$dP = 0
    fnc_CNModel$ET = 0
    fnc_CNModel$AW = 0
    fnc_CNModel$Excess = 0
    fnc_CNModel$S = 0
    fnc_CNModel$Qpred = 0
    attach(fnc_CNModel)
    SSCNavg = (1000/CNavg - 10) * 25.4
    SSCN = SoilStorage(S_avg = SSCNavg, field_capacity = func_DAWC * 
        0.9, soil_water_content = 0.1 * func_DAWC, porosity = func_DAWC)
    Ia_init = IaFrac * SSCN
    fnc_CNModel$CNavg = CNavg
    fnc_CNModel$SSCNavg = SSCNavg
    fnc_CNModel$SSCN = SSCN
    detach(fnc_CNModel)
    rm(list = c("CNavg", "SSCN", "SSCNavg"))
    fnc_CNModel$Ia = Ia_init
    attach(fnc_CNModel)
    for (t in 2:length(AW)) {
        ET[t] = AW[t - 1]/AWC[t - 1] * PET[t]
        dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + SNO_Energy$SnowMelt_mm[t] + 
            HillslopeAboveExcess[t]
        if (AW[t - 1] + dP[t] >= AWC[t]) {
            Excess[t] = AW[t - 1] + dP[t] - AWC[t]
            AW[t] = AWC[t]
        }
        else if (dP[t] <= Ia[t]) {
            Excess[t] = 0
            AW[t] = AW[t - 1] + dP[t]
        }
        else {
            Excess[t] = (dP[t] - Ia[t])^2/(dP[t] - Ia[t] + SSCN[t])
            AW[t] = AW[t - 1] + dP[t] - Excess[t]
        }
        S[t] = S[t - 1] + Excess[t]
        Qpred[t] = fnc_fcres * S[t]
        S[t] = S[t] - Qpred[t]
    }
    fnc_CNModel$ET = ET
    fnc_CNModel$dP = dP
    fnc_CNModel$AW = AW
    fnc_CNModel$Excess = Excess
    fnc_CNModel$S = S
    fnc_CNModel$Qpred = Qpred
    rm(list = c("AW", "dP", "ET", "Excess", "Qpred", "S"))
    detach(fnc_CNModel)
    return(fnc_CNModel)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
