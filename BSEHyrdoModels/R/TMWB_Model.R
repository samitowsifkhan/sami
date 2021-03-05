TMWB_Model <-
function(fnc_TMWB,fnc_slope=0, 
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
