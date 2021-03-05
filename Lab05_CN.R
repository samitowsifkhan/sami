# Cleaning up
objects()
rm(list=objects())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(elevatr,raster,soilDB,rgdal)
pacman::p_load(EcoHydRology,curl,httr,rnoaa)
myflowgage_id="0205551460"  ### For LICKRUN
myflowgage=get_usgs_gage(myflowgage_id,begin_date = "2015-01-01",end_date = "2021-03-01")

myflowgage$flowdata$Qmm = myflowgage$flowdata$flow/myflowgage$area/10^3
stns=meteo_distance(
  station_data=ghcnd_stations(),
  lat=myflowgage$declat,
  long=myflowgage$declon,
  units = "deg",
  radius = 30,
  limit = NULL
)

WXStn=stns[stns$element=="TMAX"&stns$last_year>=2020,]$id[2]
WXData=meteo_pull_monitors(
  monitors=WXStn,
  keep_flags = FALSE,
  date_min = "2016-01-01",
  date_max = NULL,
  var = c("TMAX","TMIN","PRCP") 
)
summary(WXData)


modeldata=merge(WXData,myflowgage$flowdata,by.x="date",by.y="mdate")
modeldata$MaxTemp=modeldata$tmax/10 # Converting to C
modeldata$MinTemp=modeldata$tmin/10 # Converting to C
modeldata$P=modeldata$prcp/10 # Converting to mm

mean(modeldata$Qmm)
mean(modeldata$P)
modeldata$P[is.na(modeldata$P)]=0
modeldata$MinTemp[is.na(modeldata$MinTemp)]=0
modeldata$MaxTemp[is.na(modeldata$MaxTemp)]=
  modeldata$MinTemp[is.na(modeldata$MaxTemp)] +1
modeldata$MaxTemp[modeldata$MaxTemp<=modeldata$MinTemp]=
  modeldata$MinTemp[modeldata$MaxTemp<=modeldata$MinTemp]+1
modeldata$AvgTemp=(modeldata$MaxTemp+modeldata$MinTemp)/2.0


soildrying<-function(AWprev,dP,AWC){
  AW<-AWprev*exp(dP/AWC)
  excess<-0.0
  c(AW,excess)
}
soil_wetting_above_capacity<-function(AWprev,dP,AWC){
  AW<-AWC
  excess<-AWprev+dP-AWC
  c(AW,excess)
}
soilwetting<-function(AWprev,dP,AWC){
  AW<-AWprev+dP
  excess<-0.0
  c(AW,excess)
}

# Start of Model
TMWB_Model=function(fnc_TMWB,fnc_slope=0, 
                    fnc_aspect=0,func_DAWC=.3,
                    func_z=1000,fnc_fcres=.3) {
  
  attach(fnc_TMWB)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1, 
                      windHt = 2, groundAlbedo = 0.25,
                      SurfEmissiv = 0.95, windSp = 2, forest = 0, 
                      startingSnowDepth_m = 0, startingSnowDensity_kg_m3=450)
  
  detach(fnc_TMWB)
  fnc_TMWB$SNO=SNO_Energy$SnowWaterEq_mm
  fnc_TMWB$SNOmlt=SNO_Energy$SnowMelt_mm
  attach(fnc_TMWB)
  fnc_TMWB$Albedo=.23
  fnc_TMWB$Albedo[fnc_TMWB$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),Tmax_C = MaxTemp,Tmin_C = MinTemp,lat_radians = myflowgage$declat*pi/180) * 1000
  fnc_TMWB$PET=PET
  detach(fnc_TMWB)
  rm(list="PET")
  
  fnc_TMWB$AWC=func_DAWC*func_z
  
  fnc_TMWB$dP = 0 # Initializing Net Precipitation
  fnc_TMWB$ET = 0 # Initializing ET
  fnc_TMWB$AW = 0 # Initializing AW
  fnc_TMWB$Excess = 0 # Initializing Excess
  
  
  
  attach(fnc_TMWB)
  for (t in 2:length(AW)){
    
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + 
      SNO_Energy$SnowMelt_mm[t] + HillslopeAboveExcess[t]
    
    
    if (dP[t]<=0) {
      values<-soildrying(AW[t-1],dP[t],AWC[t])
    } else if((dP[t]>0) & (AW[t-1]+dP[t])<=AWC[t]) {
      values<-soilwetting(AW[t-1],dP[t],AWC[t])
    } else {
      values<-soil_wetting_above_capacity(AW[t-1],dP[t],AWC[t])
    }
    AW[t]<-values[1]
    Excess[t]<-values[2]
  }
  fnc_TMWB$AW=AW
  fnc_TMWB$Excess=Excess
  fnc_TMWB$dP=dP
  fnc_TMWB$ET=ET
  detach(fnc_TMWB) # IMPORTANT TO DETACH
  rm(list=c("AW", "dP", "ET", "Excess"))
  
  fnc_TMWB$Qpred=NA
  fnc_TMWB$Qpred[1]=0
  fnc_TMWB$S=NA
  fnc_TMWB$S[1]=0
  
  fcres=fnc_fcres
  attach(fnc_TMWB)
  for (t in 2:length(Qpred)){
    S[t]=S[t-1]+Excess[t]     
    Qpred[t]=fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  fnc_TMWB$S=S
  fnc_TMWB$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  detach(fnc_TMWB) # IMPORTANT TO DETACH
  rm(list=c("Qpred", "S"))
  return(fnc_TMWB)
}

#Initializing the model and making a column to put data
modeldata$HillslopeAboveExcess=0
summary(modeldata)
TopSlope=modeldata
MidSlope=modeldata
BotSlope=modeldata

# Running the model with the three HRU dataframes
# Low slope but highest ksat
# These 3 function calls are what you will vary for the Lab 04 homework 
TopSlope = TMWB_Model(fnc_TMWB = TopSlope,fnc_slope=0, 
                      fnc_aspect=0,func_DAWC=.3,
                      func_z=500,fnc_fcres=.3)
MidSlope$HillslopeAboveExcess=TopSlope$Excess
# Higher slope, medium ksat, fcres=0.5 
MidSlope = TMWB_Model(fnc_TMWB = MidSlope,fnc_slope=0, 
                      fnc_aspect=0,func_DAWC=.3,
                      func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlope$HillslopeAboveExcess=MidSlope$Excess
BotSlope = TMWB_Model(fnc_TMWB = BotSlope,fnc_slope=0, 
                      fnc_aspect=0,func_DAWC=.3,
                      func_z=1000,fnc_fcres=.2)

#AW Plots HW1
plot(TopSlope$date,TopSlope$AW,type="l",ylim=c(0,400),col="red",xlab="Date",ylab="AW")
lines(MidSlope$date,MidSlope$AW,col="blue")
lines(BotSlope$date,BotSlope$AW,col="green")
legend(17150,400,legend=c("TopSlope","MidSlope","BottomSlope"),lty = c(1,1,1),lwd=c(1,1,1),col= c("red", "blue","green"),
       ncol=1)

# Excess Plots HW1
plot(TopSlope$date,TopSlope$Excess,type="l",ylim=c(0,250),col="red",xlab="Date",ylab="Excess")
lines(MidSlope$date,MidSlope$Excess,col="blue")
lines(BotSlope$date,BotSlope$Excess,col="green")
legend(17150,250,legend=c("TopSlope","MidSlope","BottomSlope"),lty = c(1,1,1),lwd=c(1,1,1),col= c("red", "blue","green"),
       ncol=1)

plot(BotSlope$date,BotSlope$PET,type="l",col=1,xlab="Date",ylab="(P)ET (mm)")
lines(BotSlope$date,BotSlope$ET,type="l",col=2)
lines(MidSlope$date,MidSlope$ET,type="l",col=3)
lines(TopSlope$date,TopSlope$ET,type="l",col=4)


plot(BotSlope$date,cumsum(BotSlope$Qpred),type="l",
     xlab="Date",ylab="Flow Q Cumulative Summary (mm)")
lines(MidSlope$date,cumsum(MidSlope$Qpred),col="red")
lines(TopSlope$date,cumsum(TopSlope$Qpred),col="green")

plot(BotSlope$date,BotSlope$Qpred,type="l")
NSeff(BotSlope$Qmm,BotSlope$Qpred)




# 
#CNModel
# 
CNModel<-function(fnc_CNModel, CNavg = 75,IaFrac = 0.05,fnc_slope=0, 
                  fnc_aspect=0,func_DAWC=.3,func_z=1000,fnc_fcres=.3) {
  
  # Energy Balance based Snow Accumulation 
  # and Melt model from the EcoHydRology package.
  attach(fnc_CNModel)
  SNO_Energy=SnowMelt(date, P, MaxTemp-3, MinTemp-3, myflowgage$declat, 
                      slope = fnc_slope, aspect = fnc_aspect, tempHt = 1, 
                      windHt = 2, groundAlbedo = 0.25,SurfEmissiv = 0.95, windSp = 2, 
                      forest = 0, startingSnowDepth_m = 0,startingSnowDensity_kg_m3=450)
  # We will update the -3 in the above to be a lapse rate adjustment
  detach(fnc_CNModel)
  fnc_CNModel$SNO=SNO_Energy$SnowWaterEq_mm
  fnc_CNModel$SNOmlt=SNO_Energy$SnowMelt_mm
  fnc_CNModel$SnowfallWatEq_mm=SNO_Energy$SnowfallWatEq_mm
  fnc_CNModel$SnowMelt_mm=SNO_Energy$SnowMelt_mm
  attach(fnc_CNModel)
  fnc_CNModel$Albedo=.23
  fnc_CNModel$Albedo[fnc_CNModel$SNO>0]=.95
  PET=PET_fromTemp(Jday=(1+as.POSIXlt(date)$yday),
                   Tmax_C = MaxTemp,Tmin_C = MinTemp,
                   lat_radians = myflowgage$declat*pi/180) * 1000
  fnc_CNModel$PET=PET
  detach(fnc_CNModel)
  rm(list="PET")
  
  fnc_CNModel$AWC=func_DAWC*func_z
  # Oh, this we want to vary some of these around our watershed!
  fnc_CNModel$dP = 0 # Initializing Net Precipitation
  fnc_CNModel$ET = 0 # Initializing ET
  fnc_CNModel$AW = 0 # Initializing AW
  fnc_CNModel$Excess = 0 # Initializing Excess
  fnc_CNModel$S =0 # Initializing S
  fnc_CNModel$Qpred=0 # Initializing Qpred
  
  
  attach(fnc_CNModel)
  
  
  SSCNavg=(1000/CNavg-10)*25.4
  SSCN=SoilStorage(S_avg=SSCNavg, field_capacity=func_DAWC*.9,
                   soil_water_content=0.1*func_DAWC, porosity=func_DAWC)
  Ia_init=IaFrac*SSCN   
  fnc_CNModel$CNavg = CNavg
  fnc_CNModel$SSCNavg = SSCNavg
  fnc_CNModel$SSCN = SSCN
  detach(fnc_CNModel)
  rm(list=c("CNavg", "SSCN", "SSCNavg"))
  fnc_CNModel$Ia = Ia_init
  attach(fnc_CNModel)
  # Those processes that are dependant on prior days conditions, we run as a 
  # loop through each of the days.
  for (t in 2:length(AW)){
    ET[t] = AW[t-1]/AWC[t-1]*PET[t]
    # Calculating Net Precipitation which adds in slope above's Excess
    dP[t] = SNO_Energy$Rain_mm[t] - ET[t] + 
      SNO_Energy$SnowMelt_mm[t] + HillslopeAboveExcess[t]    # CN Solution
    # Is the soil saturated, and thus can't take more dP? 
    if (AW[t-1] + dP[t]>=AWC[t]){
      Excess[t]=AW[t-1] + dP[t] -AWC[t]
      AW[t]=AWC[t]
      # Otherwise, if dP is less than the initial abstraction? 
      # https://en.wikipedia.org/wiki/Runoff_curve_number#Definition
    } else if (dP[t]<=Ia[t]) {
      Excess[t]=0.0
      AW[t]=AW[t-1] + dP[t]
    } else {
      Excess[t]=(dP[t]-Ia[t])^2/(dP[t]-Ia[t]+SSCN[t])
      AW[t]=AW[t-1] + dP[t] -Excess[t]
    }
    S[t]=S[t-1]+Excess[t]
    Qpred[t]=fnc_fcres*S[t]
    S[t]=S[t]-Qpred[t]
  }
  fnc_CNModel$ET=ET
  fnc_CNModel$dP=dP
  fnc_CNModel$AW=AW
  fnc_CNModel$Excess=Excess
  fnc_CNModel$S=S
  fnc_CNModel$Qpred=Qpred # UPDATE vector BEFORE DETACHING
  rm(list=c("AW", "dP", "ET", "Excess", "Qpred", "S"))
  detach(fnc_CNModel)
  return(fnc_CNModel)
}
modeldata$HillslopeAboveExcess=0
TopSlopeCN=modeldata
MidSlopeCN=modeldata
BotSlopeCN=modeldata
# Call the new CNModel() function with Top,Mid,BotSlope HRU objects,
# passing the Qpred into the lower HRUs HillslopeAboveExcess (as area scaled flow)
TopSlopeCN=CNModel(TopSlopeCN, CNavg = 60)
TopSlopeCN = CNModel(fnc_CNModel = TopSlopeCN, CNavg = 60,fnc_slope=0,
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=500,fnc_fcres=.3)
MidSlopeCN$HillslopeAboveExcess=TopSlopeCN$Excess
# Higher slope, medium ksat, fcres=0.5 
MidSlopeCN = CNModel(fnc_CNModel = MidSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=750,fnc_fcres=.5)
# Low Slope and lowest ksat, $fcres=0.2
BotSlopeCN$HillslopeAboveExcess=MidSlopeCN$Excess
BotSlopeCN = CNModel(fnc_CNModel = BotSlopeCN, CNavg = 60,fnc_slope=0, 
                     fnc_aspect=0,func_DAWC=.3,
                     func_z=1000,fnc_fcres=.2)


#AW Plots HW1
plot(BotSlopeCN$date,BotSlopeCN$AW,type="l",col=1,xlab="Date",ylab="AW (mm)")
lines(MidSlopeCN$date,MidSlopeCN$AW,type="l",col=2)
lines(TopSlopeCN$date,TopSlopeCN$AW,type="l",col=3)
# Excess Plots HW1
plot(BotSlopeCN$date,BotSlopeCN$Excess,type="l",col=1,xlab="Date",ylab="Excess (mm)")
lines(MidSlopeCN$date,MidSlopeCN$Excess,type="l",col=2)
lines(TopSlopeCN$date,TopSlopeCN$Excess,type="l",col=3)

# PET and ET HW2
plot(BotSlopeCN$date,BotSlopeCN$PET,type="l",col=1,xlab="Date",ylab="(P)ET (mm)")
lines(BotSlopeCN$date,BotSlopeCN$ET,type="l",col=2)
lines(MidSlopeCN$date,MidSlopeCN$ET,type="l",col=3)
lines(TopSlopeCN$date,TopSlopeCN$ET,type="l",col=4)
# or as cumulative summations
plot(TopSlopeCN$date,cumsum(BotSlopeCN$PET),type="l",
     xlab="Date",ylab="(P)ET")
lines(TopSlopeCN$date,cumsum(TopSlopeCN$ET),col="red")
lines(MidSlopeCN$date,cumsum(MidSlopeCN$ET),col="green")
lines(BotSlopeCN$date,cumsum(BotSlopeCN$ET),col="blue")


# Cumulative Summary of QPred is very informative
plot(BotSlopeCN$date,cumsum(BotSlopeCN$Qpred),type="l",
     xlab="Date",ylab="Flow Q Cumulative Summary (mm)")
lines(MidSlopeCN$date,cumsum(MidSlopeCN$Qpred),col="red")
lines(TopSlopeCN$date,cumsum(TopSlopeCN$Qpred),col="green")

# Model Performance 
plot(BotSlopeCN$date,BotSlopeCN$Qpred,type="l")
NSeff(BotSlopeCN$Qmm,BotSlopeCN$Qpred)

# finish building all the hillslope HRUs..
