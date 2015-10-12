
#Vadose recharge time series generator.
# This function reads in a timeseries of surface recharge for zones in a catchment and estimates the timeseries of vadose recharge.
# Each zone has recharge for dryland, irrigated land and from the river.
# The dryland and irrigated land can be recharging to different aquifers which contribute different amounts to the groundwater that is being modelled.
# The vadose recharge is assumed to be well estimated by an exponential weighted moving average.
#This is all taken from Bidwell and Burbery 2011 Groundwater Data Analysis - quantifying aquifer dynamics Lincoln Ventures Report No 4110/1

# Check that the input values make sense

vadose.recharge <- function(ZoneStorageTime = c(30,30,30,20),
                            AquiferZoneDryFractions = list(c(0.3,0.3,0.3,0.3),
                                                              c(0,0,0.1,0.1),
                                                              c(0.991,0.9398,0.7,0.7)),
                            AquiferZoneIrrigFractions = list(c(0,0,0,0),
                                                              c(0,0,0,0),
                                                              c(0.009,0.0602,0.0069,0.1737)),
                            RiverRechargeFractions = c(0,0,0,0),
                            RechargeFileName ="Fishcreek recharge dataV2.csv")

#***************************************************************
#  **************Description of the input paramaters***********
#ZoneStorageTime is a vector of times (in days)for each zone for how long it takes for the water to pass through the vadose zone
#
#AquiferZoneDryFractions is a list of vectors of the proportion of each aquifer's dryland recharge that is contributing to each zones vadose recharge
#AquiferZoneDryFractions <- list(c(Aquifer1Zone1Dry,Aquifer1Zone2Dry,...,Aquifer1ZonemDry),
#                          c(Aquifer2Zone1Dry,Aquifer2Zone2Dry,...,Aquifer2ZonemDry),
#                           ...
#                          c(AquifernZone1Dry,AquifernZone2Dry,AquifernZone3Dry,AquifernZonemDry))
#
#AquiferZoneDryFractions is a list of vectors of the proportion of each aquifer's irrigation recharge that is contributing to each zones vadose recharge
#AquiferZoneIrrigFractions <- list(c(Aquifer1Zone1Irrig,Aquifer1Zone2Irrig,...,Aquifer1ZonemIrrig),
#                          c(Aquifer2Zone1Irrig,Aquifer2Zone2Irrig,...,Aquifer2ZonemIrrig),
#                           ...
#                          c(AquifernZone1Irrig,AquifernZone2Irrig,AquifernZone3Irrig,AquifernZonemIrrig))
#
# Recharge filename is a csv file of the land surface recharge data combined with observed data, riverflow and pump data
# First column is date,Observed piezometric head,Observed discharge,
#             Zone1Aquifer1dryland,Zone1Aquifer1irrigated,Zone1Aquifer2Dryland,Zone1Aquifer2Irrigated,...,Zone1AquifernDryland,Zone1Aquifernirrigated,
#             Z2A1Dry,Z2A1Irr,Z2A2Dry,Z2A2Irr,...,Z2AnDry,Z2AnIrr,
#             ....
#             ZnA1Dry,ZnA1Irr,ZnA2Dry,ZnA2Irr,...,ZnAnDry,ZnAnIrr,
#             Z1River,Z2River,...,ZnRiver,
#             Z1Pumped,Z2Pumped,...,ZnPumped
#********************************************************************
  
{  #Start of the function
  
#Load libraries
library(TTR)          #This library includes the EMA exponentially weighted moving average

#read in the data
ZoneTimeseries        <-  read.csv(RechargeFileName)                              #this is the daily timeseries of observed discharge and groundwater level, vadose recharge for each zone, river recharge to groundwater and groundwater pumping in mm

#Calculate the number of zones and aquifers based on the 
NumberOfZones         <- length(AquiferZoneDryFractions[[1]])
NumberOfAquifers      <- length(AquiferZoneDryFractions)

#Multiply the aquifer recharge fractions by the recharge timeseries to calculate a total recharge for each zone
ZoneFractions  <- c()                                               #Initialise ZoneFractions to an empty set
for (ZoneNo in 1:NumberOfZones){
  for (AquiferNo in 1:NumberOfAquifers){
    ZoneFractions <- c(ZoneFractions,AquiferZoneDryFractions[[AquiferNo]][ZoneNo],AquiferZoneIrrigFractions[[AquiferNo]][ZoneNo])
  }
}

#Multiply the recharge timeseries by their respective fractions
ZoneSurfaceRecharge  <- sweep(ZoneTimeseries[4:(ncol(ZoneTimeseries)-2*NumberOfZones)],MARGIN=2,ZoneFractions,'*')
ZoneRiverRecharge <- sweep(ZoneTimeseries[(ncol(ZoneTimeseries)-2*NumberOfZones+1):(ncol(ZoneTimeseries)-NumberOfZones)],MARGIN=2,RiverRechargeFractions,'*')

#Get the sums for each zone so that we are left with just one timeseries for each zone
ZoneSurfaceRechargeTotals    <- c()     #Initialise
for (ZoneNo in 1:(NumberOfZones)){
  ZoneSurfaceRechargeTotals <- cbind(ZoneSurfaceRechargeTotals,rowSums(ZoneSurfaceRecharge[,c(((ZoneNo-1)*(2*NumberOfAquifers)+1):(ZoneNo*2*NumberOfAquifers))]))
}

#Add in the river recharge
ZoneRechargeTotals <- ZoneSurfaceRechargeTotals + ZoneRiverRecharge

#build a list of lists of the recharge data with the weighting time constant, needed to apply the exp weighted moving average
RechargeList <- list()
for (ZoneNo in 1:(NumberOfZones)){
  RechargeList[[length(RechargeList)+1]]<-list(series=ZoneRechargeTotals[,ZoneNo],timeConstant=ZoneStorageTime[ZoneNo])
}

#Apply the exponential weighted moving average
ZoneVadoseRecharge <- sapply(RechargeList, function(x) EMA(x$series,n=1,ratio=1-exp(-1/max(x$timeConstant,0.000001))))
                          
#Take out any pumping
ZoneVadoseRecharge <- ZoneVadoseRecharge - ZoneTimeseries[,(4+NumberOfZones*(2*NumberOfAquifers+1)):(3+NumberOfZones*(2*NumberOfAquifers+2))] 

#Convert from millimetres to metres
ZoneVadoseRecharge<- ZoneVadoseRecharge / 1000

#Put the dates back on
rownames(ZoneVadoseRecharge) <- ZoneTimeseries[,1]  

#Put the observed data back on
ZoneVadoseRecharge <- cbind(ZoneTimeseries[,2:3],ZoneVadoseRecharge) 

return(ZoneVadoseRecharge)
} #End of function