#This script prepares estimates of groundwater discharge for Takaka streams using the Groundwater Data Analysis tools functions

library(zoo)

#Ideally, in the future, the next two lines could be replaced with the call to a package 
source("//aqualinc-sbs/data/ARL Projects/Other/C15066_WaterWheel2/Modeling/Groundwater/R eigenmodels/BoussEigen.R")
source("//aqualinc-sbs/data/ARL Projects/Other/C15066_WaterWheel2/Modeling/Groundwater/R eigenmodels/Vadose.R")

#Get the observational data
Observeddata <- read.csv("Observeddata.csv")
ObserveddataZoo <- read.zoo(Observeddata,format="%d/%m/%Y")



fishRecharge<-vadose.recharge(ZoneStorageTime=c(4,3,3,3),AquiferZoneDryFractions = list(c(0.8,0.8,0.8,0.8),c(0,0,0,0),c(0,0,0.193,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(1,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
fishDischarge <- bouss.eigen(WellDistance=66000,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.01,Transmisivity=313632,DischargeScaleFactor=950,RechargeData=fishRecharge,GWBypassFlow=1)
bouss.plot(fishDischarge[,4],ObserveddataZoo$FishCreekDischarge)

PupuMainSpringRecharge<-vadose.recharge(ZoneStorageTime=c(3,2,2,2),AquiferZoneDryFractions = list(c(1,1,1,1),c(0,0,0,0),c(0,0,0.143,0)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.0)),RiverRechargeFractions=c(1,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
PupuMainSpringdischarge <- bouss.eigen(WellDistance=66000,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.05,Transmisivity=228690,DischargeScaleFactor=3250,RechargeData=PupuMainSpringRecharge,GWBypassFlow=1.9)
bouss.plot(PupuMainSpringdischarge[,4],ObserveddataZoo$PupuMainSpringDischarge)

SpringRiverRecharge<-vadose.recharge(ZoneStorageTime=c(4,3,3,0.7),AquiferZoneDryFractions = list(c(1,1,1,1),c(0,0,0,0),c(0,0,0.193,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(0,0,1,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
SpringRiverdischarge <- bouss.eigen(WellDistance=66000,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.01,Transmisivity=108900,DischargeScaleFactor=2350,RechargeData=SpringRiverRecharge,GWBypassFlow=0.5)
bouss.plot(SpringRiverdischarge[,4],ObserveddataZoo$SpringRiver)

MotupipiRecharge<-vadose.recharge(ZoneStorageTime=c(1,0,0,0),AquiferZoneDryFractions = list(c(0,0,0,0),c(0,0,1,1),c(0,0,0.08,0.36)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(0,0,0.1,0.1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
Motupipidischarge <- bouss.eigen(WellDistance=66000,ZoneLengths=c(0,0,2508,40920),Storativity=0.005,Transmisivity=54450,DischargeScaleFactor=100,RechargeData=MotupipiRecharge,GWBypassFlow=-0.26)
bouss.plot(Motupipidischarge[,4],ObserveddataZoo$MotupipiDischarge,startDate="2005-01-01")

PaynesFordRecharge<-vadose.recharge(ZoneStorageTime=c(0,0,0,0),AquiferZoneDryFractions = list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.033,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0.0602,0.0069,0.1737)),RiverRechargeFractions=c(0,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
PaynesForddischarge <- bouss.eigen(WellDistance=66000,ZoneLengths=c(0,23100,19800,23100),Storativity=0.001,Transmisivity=653400,DischargeScaleFactor=1450,RechargeData=PaynesFordRecharge,GWBypassFlow=2.45)
bouss.plot(PaynesForddischarge[,4],ObserveddataZoo$PaynesFord)

BallRecharge<-vadose.recharge(ZoneStorageTime=c(4,3,2,1),AquiferZoneDryFractions = list(c(1,1,1,1),c(0,0,0,0),c(0,0,0.193,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(0,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
BallGWSim <- bouss.eigen(WellDistance=64680,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.047,Transmisivity=22500,DischargeScaleFactor=1450,RechargeData=BallRecharge,GWBypassFlow=-0,InitialEigenState=24.2491)
bouss.plot(BallGWSim$estimated_groundwater_depth,ObserveddataZoo$Ball, startDate="1994-01-01",endDate="2005-01-01")

BennettRecharge<-vadose.recharge(ZoneStorageTime=c(4,0,0,0),AquiferZoneDryFractions = list(c(1,1,1,1),c(0,0,0,0),c(0.9910,0.9398,0.9931,0)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0.009,0.0602,0.0069,0)),RiverRechargeFractions=c(0,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
BennettGWSim <- bouss.eigen(WellDistance=44220,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.523,Transmisivity=487961,DischargeScaleFactor=1450,RechargeData=BennettRecharge,GWBypassFlow=-0,InitialEigenState=30.0072)
bouss.plot(BennettGWSim$estimated_groundwater_depth,ObserveddataZoo$Bennett, startDate="1996-01-01",endDate="2000-01-01")

HamamaRecharge<-vadose.recharge(ZoneStorageTime=c(4,3,2,1),AquiferZoneDryFractions = list(c(1,1,1,1),c(0,0,0,0),c(0.9910,0.9398,0.7931,0)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0.009,0.0602,0.0069,0)),RiverRechargeFractions=c(0,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
HamamaGWSim <- bouss.eigen(WellDistance=49500,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.17,Transmisivity=292505,DischargeScaleFactor=1,RechargeData=HamamaRecharge,GWBypassFlow=-0,InitialEigenState=11.9916)
bouss.plot(HamamaGWSim$estimated_groundwater_depth,ObserveddataZoo$Hamama, startDate="1988-01-01",endDate="2006-01-01")

PupuMainSpringRecharge<-vadose.recharge(ZoneStorageTime=c(2,1,0,0),AquiferZoneDryFractions = list(c(1,1,1,1),c(0,0,0,0),c(0,0,0.143,0)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0)),RiverRechargeFractions=c(0,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
PupuMainSpringGWSim <- bouss.eigen(WellDistance=48180,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.5,Transmisivity=258093,DischargeScaleFactor=1,RechargeData=PupuMainSpringRecharge,GWBypassFlow=-0,InitialEigenState=4.6020)
bouss.plot(PupuMainSpringGWSim$estimated_groundwater_depth,ObserveddataZoo$PupuMainSpring, startDate="1999-01-01")

SavageRecharge<-vadose.recharge(ZoneStorageTime=c(2,1,0,0),AquiferZoneDryFractions = list(c(1,1,1,1),c(0,0,0,0),c(0.991,0.9398,0.2931,0)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0.009,0.0602,0.0069,0)),RiverRechargeFractions=c(0,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
SavageGWSim <- bouss.eigen(WellDistance=48180,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.05,Transmisivity=210177,DischargeScaleFactor=1,RechargeData=SavageRecharge,GWBypassFlow=-0,InitialEigenState=4.2364)
bouss.plot(SavageGWSim$estimated_groundwater_depth,ObserveddataZoo$Savage, startDate="2002-01-01")

SowmanRecharge<-vadose.recharge(ZoneStorageTime=c(2,1,0,0),AquiferZoneDryFractions = list(c(1,1,1,1),c(0,0,0,0),c(0.08,0.36,0.2931,0)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0.009,0.0602,0.0069,0)),RiverRechargeFractions=c(0,1,0.958,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
SowmanGWSim <- bouss.eigen(WellDistance=46200,ZoneLengths=c(35600,5700,8900,15800),Storativity=0.04,Transmisivity=217800,DischargeScaleFactor=1,RechargeData=SowmanRecharge,GWBypassFlow=-0,InitialEigenState=3.7282)
bouss.plot(SowmanGWSim$estimated_groundwater_depth,ObserveddataZoo$Sowman, startDate="1999-01-01")

CserneyRecharge<-vadose.recharge(ZoneStorageTime=c(0,0,4,3),AquiferZoneDryFractions = list(c(0,0,0,0),c(0,0,1,1),c(0,0,0.08,0.36)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(0,0,0.1,0.1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
CserneyGWSim <- bouss.eigen(WellDistance=56100,ZoneLengths=c(0,0,25255.10202,40744.89797),Storativity=0.033,Transmisivity=201247,DischargeScaleFactor=1,RechargeData=CserneyRecharge,GWBypassFlow=-0,InitialEigenState=1)
bouss.plot(CserneyGWSim$estimated_groundwater_depth,ObserveddataZoo$Cserney, startDate="1987-01-01",endDate="2007-01-01")

GroveOrchardRecharge<-vadose.recharge(ZoneStorageTime=c(0,0,6,5),AquiferZoneDryFractions = list(c(0,0,0,0),c(0,0,1,1),c(0,0,0.08,0.36)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.0069,0.1737)),RiverRechargeFractions=c(0,0,0,0),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
GroveOrchardGWSim <- bouss.eigen(WellDistance=56100,ZoneLengths=c(0,0,25255.10202,40744.89797),Storativity=0.075,Transmisivity=69260,DischargeScaleFactor=1,RechargeData=GroveOrchardRecharge,GWBypassFlow=-0,InitialEigenState=10.5927)
bouss.plot(GroveOrchardGWSim$estimated_groundwater_depth,ObserveddataZoo$GroveOrchard, startDate="1987-01-01",endDate="2007-01-01")

MotupipiSubstationRecharge<-vadose.recharge(ZoneStorageTime=c(0,0,4,3),AquiferZoneDryFractions = list(c(0,0,0,0),c(0,0,1,1),c(0,0,0.08,0.36)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0.0602,0.0069,0.1737)),RiverRechargeFractions=c(0,0,0.1,0.1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
MotupipiSubstationGWSim <- bouss.eigen(WellDistance=52800,ZoneLengths=c(0,0,25255.10202,40744.89797),Storativity=0.009,Transmisivity=686070,DischargeScaleFactor=1,RechargeData=MotupipiSubstationRecharge,GWBypassFlow=-0,InitialEigenState=0.095)
bouss.plot(MotupipiSubstationGWSim$estimated_groundwater_depth,ObserveddataZoo$MotupipiSubstation, startDate="1981-01-01",endDate="1989-01-01")

FireStationRecharge<-vadose.recharge(ZoneStorageTime=c(0,0,0,0),AquiferZoneDryFractions = list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.033,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0.0602,0.0069,0.1737)),RiverRechargeFractions=c(0,0,0.96,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
FireStationGWSim <- bouss.eigen(WellDistance=64350,ZoneLengths=c(0,23175.57,19396.938,23427.492),Storativity=0.15,Transmisivity=79715,DischargeScaleFactor=1,RechargeData=FireStationRecharge,GWBypassFlow=-0,InitialEigenState=30)
bouss.plot(FireStationGWSim$estimated_groundwater_depth,ObserveddataZoo$FireStation, startDate="2004-01-01")

JeffersonRecharge<-vadose.recharge(ZoneStorageTime=c(0,0,0,0),AquiferZoneDryFractions = list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.033,0.8263)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0.0602,0.0069,0.1737)),RiverRechargeFractions=c(0,0,0.96,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
JeffersonGWSim <- bouss.eigen(WellDistance=33660,ZoneLengths=c(0,23175.57,19396.938,23427.492),Storativity=0.13,Transmisivity=302960,DischargeScaleFactor=1,RechargeData=JeffersonRecharge,GWBypassFlow=-0,InitialEigenState=5.7669)
bouss.plot(JeffersonGWSim$estimated_groundwater_depth,ObserveddataZoo$Jefferson, startDate="1998-01-01",endDate="2003-01-01")

TDCOfficeRecharge<-vadose.recharge(ZoneStorageTime=c(0,0,0,0),AquiferZoneDryFractions = list(c(0,0,0,0),c(0,0,0,0),c(0,0,0.033,0.82632168)),AquiferZoneIrrigFractions=list(c(0,0,0,0),c(0,0,0,0),c(0,0.0602,0.0069,0.17367832)),RiverRechargeFractions=c(0,0,0.96,1),RechargeFileName="GoldenBayLandSurfaceRechargeData.csv")
TDCOfficeGWSim <- bouss.eigen(WellDistance=64680,ZoneLengths=c(0,23175.57,19396.938,23427.492),Storativity=0.16,Transmisivity=65514,DischargeScaleFactor=1,RechargeData=TDCOfficeRecharge,GWBypassFlow=-0,InitialEigenState=20)
bouss.plot(TDCOfficeGWSim$estimated_groundwater_depth,ObserveddataZoo$TDCOffice, startDate="1999-01-01")
