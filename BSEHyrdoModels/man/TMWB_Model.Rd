\name{TMWB_Model}
\alias{TMWB_Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{"TMWB_Model"
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
TMWB_Model(fnc_TMWB, fnc_slope = 0, fnc_aspect = 0, func_DAWC = 0.3, func_z = 1000, fnc_fcres = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fnc_TMWB}{
%%     ~~Describe \code{fnc_TMWB} here~~
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
function (fnc_TMWB, fnc_slope = 0, fnc_aspect = 0, func_DAWC = 0.3, 
    func_z = 1000, fnc_fcres = 0.3) 
{
    attach(fnc_TMWB)
    SNO_Energy = SnowMelt(date, P, MaxTemp - 3, MinTemp - 3, 
        myflowgage$declat, slope = fnc_slope, aspect = fnc_aspect, 
        tempHt = 1, windHt = 2, groundAlbedo = 0.25, SurfEmissiv = 0.95, 
        windSp = 2, forest = 0, startingSnowDepth_m = 0, startingSnowDensity_kg_m3 = 450)
    detach(fnc_TMWB)
    fnc_TMWB$SNO = SNO_Energy$SnowWaterEq_mm
    fnc_TMWB$SNOmlt = SNO_Energy$SnowMelt_mm
    attach(fnc_TMWB)
    fnc_TMWB$Albedo = 0.23
    fnc_TMWB$Albedo[fnc_TMWB$SNO > 0] = 0.95
    PET = PET_fromTemp(Jday = (1 + as.POSIXlt(date)$yday), Tmax_C = MaxTemp, 
        Tmin_C = MinTemp, lat_radians = myflowgage$declat * pi/180) * 
        1000
    fnc_TMWB$PET = PET
    detach(fnc_TMWB)
    rm(list = "PET")
    fnc_TMWB$AWC = func_DAWC * func_z
    fnc_TMWB$dP = 0
    fnc_TMWB$ET = 0
    fnc_TMWB$AW = 0
    fnc_TMWB$Excess = 0
    attach(fnc_TMWB)
    for (t in 2:length(AW)) {
        ET[t] = AW[t - 1]/AWC[t - 1] * PET[t]
        dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + SNO_Energy$SnowMelt_mm[t] + 
            HillslopeAboveExcess[t]
        if (dP[t] <= 0) {
            values <- soildrying(AW[t - 1], dP[t], AWC[t])
        }
        else if ((dP[t] > 0) & (AW[t - 1] + dP[t]) <= AWC[t]) {
            values <- soilwetting(AW[t - 1], dP[t], AWC[t])
        }
        else {
            values <- soil_wetting_above_capacity(AW[t - 1], 
                dP[t], AWC[t])
        }
        AW[t] <- values[1]
        Excess[t] <- values[2]
    }
    fnc_TMWB$AW = AW
    fnc_TMWB$Excess = Excess
    fnc_TMWB$dP = dP
    fnc_TMWB$ET = ET
    detach(fnc_TMWB)
    rm(list = c("AW", "dP", "ET", "Excess"))
    fnc_TMWB$Qpred = NA
    fnc_TMWB$Qpred[1] = 0
    fnc_TMWB$S = NA
    fnc_TMWB$S[1] = 0
    fcres = fnc_fcres
    attach(fnc_TMWB)
    for (t in 2:length(Qpred)) {
        S[t] = S[t - 1] + Excess[t]
        Qpred[t] = fcres * S[t]
        S[t] = S[t] - Qpred[t]
    }
    fnc_TMWB$S = S
    fnc_TMWB$Qpred = Qpred
    detach(fnc_TMWB)
    rm(list = c("Qpred", "S"))
    return(fnc_TMWB)
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
