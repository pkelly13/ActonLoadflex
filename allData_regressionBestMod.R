#Try to fit model using calibration data for just one year and sequentially build from there -- compare to the model using all 20 years as a comparison to see whic hworks best

library(knitr)
opts_chunk$set(echo=TRUE, message=FALSE)
library(loadflex)
library(rloadest)

#NEED TO COLLECT CORRECTED FRACTION DATA AND MAYBE RESIDUALS AS WELL
#load Acton concentration data
setwd('~/Documents/Miami U/Acton load data/R')
fm<-read.csv('AllConcs9414.csv') #four mile creek

#load discharge data
q.daily<-read.csv('AllHDMQ9414_daily.csv')
q.hour<-read.csv('AllHDMQ9414_hourly.csv')
q.monthly<-read.csv('AllHDMQ9414_monthly.csv')

#fix date so its in POSIX format
fm$FM.Date_time<-ymd_hm(fm$FM.Date_time)
q.hour$Date.Time<-ymd_hm(q.hour$Date.Time)

#Start with four mile creek
#merge concentrations and solute data
fm.data<-merge(fm,q.hour,by.x='FM.Date_time',by.y='Date.Time')
#remove blank days
#fm.data<-fm.data[fm.data$dateTime!='',]

#remove NA columns and unused MB and LFM columns
fm.data<-fm.data[,-c(6,7,8,9,10,11,12,13,14,15,17,18,19,20,21,22)]

#rename columns
colnames(fm.data)<-c('dateTime','NH4.ug.L','NO3.ug.L','SRP.ug.L','SS.mg.L','Q.m3.sec')

#remove NAs
fm.data<-fm.data[!is.na(fm.data$dateTime),]

#use every seventh day
regr<-TRUE
for(i in 2:nrow(fm.data)){
	if(as.Date(fm.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(fm.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr<-c(regr,rep(FALSE,8))
regr[is.na(regr)]=FALSE
fm.data$REGR=regr

#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
#fm.data$REGR[fm.data$NH4.ug.L>344]=TRUE

#add year to data frame
fm.data$year<-format(ymd_hms(fm.data$dateTime),'%Y')

#Trim down to just work with NH4 first
fm.nh4<-fm.data[,c(1,2,6,7,8)]
#remove NAs
fm.nh4<-fm.nh4[!is.na(fm.nh4$NH4.ug.L),]
fm.nh4<-fm.nh4[!is.na(fm.nh4$Q.m3.sec),]

#use only positive numbers
fm.nh4<-fm.nh4[fm.nh4$NH4.ug.L>0,]
fm.nh4<-fm.nh4[fm.nh4$Q.m3.sec>0,]

#subset data for interpolation purposes
fm.reg.nh4<-subset(fm.nh4,REGR)
fm.discharge<-data.frame(dateTime=q.hour$Date.Time,Q.m3.sec=q.hour$FMQ.m3.sec)
fm.discharge<-fm.discharge[!is.na(fm.discharge$Q.m3.sec),]
fm.discharge$Q.m3.sec[fm.discharge$Q.m3.sec<=0]=0.0001
#add year to discharge
fm.discharge$year<-format(ymd_hms(fm.discharge$dateTime),'%Y')

#remove NAs
fm.reg.nh4<-fm.reg.nh4[!is.na(fm.reg.nh4$NH4.ug.L),]
fm.reg.nh4<-fm.reg.nh4[!is.na(fm.reg.nh4$Q.m3.sec),]
#remove anything with zeros
fm.reg.nh4<-fm.reg.nh4[fm.reg.nh4$Q.m3.sec>0,]
fm.reg.nh4<-fm.reg.nh4[fm.reg.nh4$NH4.ug.L>0,]

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='NH4.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')


mods<-list()
lr.1<-loadReg2(loadReg(NH4.ug.L~model(1),data=fm.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NH4.ug.L~model(2),data=fm.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NH4.ug.L~model(3),data=fm.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NH4.ug.L~model(4),data=fm.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NH4.ug.L~model(5),data=fm.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NH4.ug.L~model(6),data=fm.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NH4.ug.L~model(7),data=fm.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NH4.ug.L~model(8),data=fm.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NH4.ug.L~model(9),data=fm.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

fm.lr.nh4.flux<-predictSolute(mods[[which(aics==min(aics))]],'conc',fm.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Now do the same for NO3
#need new regression dataset
#use every seventh day
regr<-TRUE
for(i in 2:nrow(fm.data)){
	if(as.Date(fm.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(fm.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr<-c(regr,rep(FALSE,8))
regr[is.na(regr)]=FALSE
fm.data$REGR=regr

#mean(fm.data$NO3.ug.L,na.rm=T)+(sd(fm.data$NO3.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
fm.data$REGR[fm.data$NO3.ug.L>=12350]=TRUE

#Trim down to just work with NO3 
fm.no3<-fm.data[,c(1,3,6,7,8)]
#remove NAs
fm.no3<-fm.no3[!is.na(fm.no3$NO3.ug.L),]
fm.no3<-fm.no3[!is.na(fm.no3$Q.m3.sec),]

#use only positive numbers
fm.no3<-fm.no3[fm.no3$NO3.ug.L>0,]
fm.no3<-fm.no3[fm.no3$Q.m3.sec>0,]

#subset data for interpolation purposes
fm.reg.no3<-subset(fm.no3,REGR)

#create a metadata description of the dataset and desired output
meta.no3<-metadata(constituent='NO3.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

mods<-list()
lr.1<-loadReg2(loadReg(NO3.ug.L~model(1),data=fm.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NO3.ug.L~model(2),data=fm.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NO3.ug.L~model(3),data=fm.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NO3.ug.L~model(4),data=fm.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NO3.ug.L~model(5),data=fm.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NO3.ug.L~model(6),data=fm.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NO3.ug.L~model(7),data=fm.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NO3.ug.L~model(8),data=fm.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NO3.ug.L~model(9),data=fm.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

fm.lr.no3.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',fm.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Now do the same for SRP
#need new regression dataset
#use every seventh day
regr<-TRUE
for(i in 2:nrow(fm.data)){
	if(as.Date(fm.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(fm.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr<-c(regr,rep(FALSE,8))
regr[is.na(regr)]=FALSE
fm.data$REGR=regr

mean(fm.data$SRP.ug.L,na.rm=T)+(sd(fm.data$SRP.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
fm.data$REGR[fm.data$SRP.ug.L>=182]=TRUE

#Trim down to just work with NO3 
fm.srp<-fm.data[,c(1,4,6,7,8)]
#remove NAs
fm.srp<-fm.srp[!is.na(fm.srp$SRP.ug.L),]
fm.srp<-fm.srp[!is.na(fm.srp$Q.m3.sec),]

#use only positive numbers
fm.srp<-fm.srp[fm.srp$SRP.ug.L>0,]
fm.srp<-fm.srp[fm.srp$Q.m3.sec>0,]

#subset data for interpolation purposes
fm.reg.srp<-subset(fm.srp,REGR)

#create a metadata description of the dataset and desired output
meta.srp<-metadata(constituent='SRP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')
mods<-list()
lr.1<-loadReg2(loadReg(SRP.ug.L~model(1),data=fm.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SRP.ug.L~model(2),data=fm.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SRP.ug.L~model(3),data=fm.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SRP.ug.L~model(4),data=fm.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SRP.ug.L~model(5),data=fm.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SRP.ug.L~model(6),data=fm.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SRP.ug.L~model(7),data=fm.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SRP.ug.L~model(8),data=fm.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SRP.ug.L~model(9),data=fm.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

fm.lr.srp.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',fm.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Now do the same for SS
#need new regression dataset
#use every seventh day
regr<-TRUE
for(i in 2:nrow(fm.data)){
	if(as.Date(fm.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(fm.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=8){
		regr[i]=TRUE
	}
}
regr<-c(regr,rep(FALSE,6))
regr[is.na(regr)]=FALSE
fm.data$REGR=regr

mean(fm.data$SS.mg.L,na.rm=T)+(sd(fm.data$SS.mg.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
#fm.data$REGR[fm.data$SS.mg.L>=637]=TRUE

#Trim down to just work with SS
fm.ss<-fm.data[,c(1,5,6,7,8)]
#remove NAs
fm.ss<-fm.ss[!is.na(fm.ss$SS.mg.L),]
fm.ss<-fm.ss[!is.na(fm.ss$Q.m3.sec),]

#use only positive numbers
fm.ss<-fm.ss[fm.ss$SS.mg.L>0,]
fm.ss<-fm.ss[fm.ss$Q.m3.sec>0,]

fm.ss<-fm.ss[-grep('1997-02-03',fm.ss$dateTime),]

#subset data for interpolation purposes
fm.reg.ss<-subset(fm.ss,REGR)
fm.reg.ss<-fm.reg.ss[-grep('1997-02-11',fm.reg.ss$dateTime),]

#create a metadata description of the dataset and desired output
meta.ss<-metadata(constituent='SS.mg.L',flow='Q.m3.sec',dates='dateTime',conc.units='mg L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')
mods<-list()
lr.1<-loadReg2(loadReg(SS.mg.L~model(1),data=fm.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SS.mg.L~model(2),data=fm.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SS.mg.L~model(3),data=fm.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SS.mg.L~model(4),data=fm.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SS.mg.L~model(5),data=fm.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SS.mg.L~model(6),data=fm.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SS.mg.L~model(7),data=fm.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SS.mg.L~model(8),data=fm.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SS.mg.L~model(9),data=fm.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

fm.lr.ss.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',fm.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Little Four Mile Creek
setwd('~/Documents/Miami U/Acton load data/R')
lf<-read.csv('AllConcs9414_LF.csv') #four mile creek

#fix date so its in POSIX format
lf$LF.Date_time<-ymd_hm(lf$LF.Date_time)

#Start with four mile creek
#merge concentrations and solute data
lf.data<-merge(lf,q.hour,by.x='LF.Date_time',by.y='Date.Time')
#remove blank days
#fm.data<-fm.data[fm.data$dateTime!='',]

#remove NA columns and unused MB and LFM columns
lf.data<-lf.data[,-c(6,7,8,9,11,12,13,14)]

#rename columns
colnames(lf.data)<-c('dateTime','NH4.ug.L','NO3.ug.L','SRP.ug.L','SS.mg.L','Q.m3.sec')

#remove NAs
lf.data<-lf.data[!is.na(lf.data$dateTime),]

#use every seventh day
regr<-TRUE
for(i in 2:nrow(lf.data)){
	if(as.Date(lf.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(lf.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr[is.na(regr)]=FALSE
lf.data$REGR=regr

#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
#lf.data$REGR[lf.data$NH4.ug.L>344]=TRUE

#add year to data frame
lf.data$year<-format(ymd_hms(lf.data$dateTime),'%Y')

#Trim down to just work with NH4 first
lf.nh4<-lf.data[,c(1,2,6,7,8)]
#remove NAs
lf.nh4<-lf.nh4[!is.na(lf.nh4$NH4.ug.L),]
lf.nh4<-lf.nh4[!is.na(lf.nh4$Q.m3.sec),]

#use only positive numbers
lf.nh4<-lf.nh4[lf.nh4$NH4.ug.L>0,]
lf.nh4<-lf.nh4[lf.nh4$Q.m3.sec>0,]

#subset data for interpolation purposes
lf.reg.nh4<-subset(lf.nh4,REGR)
lf.discharge<-data.frame(dateTime=q.hour$Date.Time,Q.m3.sec=q.hour$LFQ.m3.sec)
lf.discharge<-lf.discharge[!is.na(lf.discharge$Q.m3.sec),]
lf.discharge$Q.m3.sec[lf.discharge$Q.m3.sec<=0]=0.0001
#add year to discharge
lf.discharge$year<-format(ymd_hms(lf.discharge$dateTime),'%Y')

#remove NAs
lf.reg.nh4<-lf.reg.nh4[!is.na(lf.reg.nh4$NH4.ug.L),]
lf.reg.nh4<-lf.reg.nh4[!is.na(lf.reg.nh4$Q.m3.sec),]
#remove anything with zeros
lf.reg.nh4<-lf.reg.nh4[lf.reg.nh4$Q.m3.sec>0,]
lf.reg.nh4<-lf.reg.nh4[lf.reg.nh4$NH4.ug.L>0,]

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='NH4.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')


mods<-list()
lr.1<-loadReg2(loadReg(NH4.ug.L~model(1),data=lf.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NH4.ug.L~model(2),data=lf.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NH4.ug.L~model(3),data=lf.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NH4.ug.L~model(4),data=lf.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NH4.ug.L~model(5),data=lf.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NH4.ug.L~model(6),data=lf.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NH4.ug.L~model(7),data=lf.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NH4.ug.L~model(8),data=lf.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NH4.ug.L~model(9),data=lf.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

lf.lr.nh4.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',lf.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Now do the same for NO3
#need new regression dataset
#use every seventh day

#mean(lf.data$NO3.ug.L,na.rm=T)+(sd(lf.data$NO3.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
lf.data$REGR[lf.data$NO3.ug.L>=12350]=TRUE

#Trim down to just work with NO3 
lf.no3<-lf.data[,c(1,3,6,7,8)]
#remove NAs
lf.no3<-lf.no3[!is.na(lf.no3$NO3.ug.L),]
lf.no3<-lf.no3[!is.na(lf.no3$Q.m3.sec),]

#use only positive numbers
lf.no3<-lf.no3[lf.no3$NO3.ug.L>0,]
lf.no3<-lf.no3[lf.no3$Q.m3.sec>0,]

#subset data for interpolation purposes
lf.reg.no3<-subset(lf.no3,REGR)

#create a metadata description of the dataset and desired output
meta.no3<-metadata(constituent='NO3.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

mods<-list()
lr.1<-loadReg2(loadReg(NO3.ug.L~model(1),data=lf.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NO3.ug.L~model(2),data=lf.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NO3.ug.L~model(3),data=lf.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NO3.ug.L~model(4),data=lf.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NO3.ug.L~model(5),data=lf.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NO3.ug.L~model(6),data=lf.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NO3.ug.L~model(7),data=lf.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NO3.ug.L~model(8),data=lf.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NO3.ug.L~model(9),data=lf.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

lf.lr.no3.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',lf.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Now do the same for SRP
#need new regression dataset
#use every seventh day

mean(lf.data$SRP.ug.L,na.rm=T)+(sd(lf.data$SRP.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
lf.data$REGR[lf.data$SRP.ug.L>=182]=TRUE

#Trim down to just work with NO3 
lf.srp<-lf.data[,c(1,4,6,7,8)]
#remove NAs
lf.srp<-lf.srp[!is.na(lf.srp$SRP.ug.L),]
lf.srp<-lf.srp[!is.na(lf.srp$Q.m3.sec),]

#use only positive numbers
lf.srp<-lf.srp[lf.srp$SRP.ug.L>0,]
lf.srp<-lf.srp[lf.srp$Q.m3.sec>0,]

#subset data for interpolation purposes
lf.reg.srp<-subset(lf.srp,REGR)

#create a metadata description of the dataset and desired output
meta.srp<-metadata(constituent='SRP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')
mods<-list()
lr.1<-loadReg2(loadReg(SRP.ug.L~model(1),data=lf.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SRP.ug.L~model(2),data=lf.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SRP.ug.L~model(3),data=lf.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SRP.ug.L~model(4),data=lf.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SRP.ug.L~model(5),data=lf.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SRP.ug.L~model(6),data=lf.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SRP.ug.L~model(7),data=lf.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SRP.ug.L~model(8),data=lf.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SRP.ug.L~model(9),data=lf.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

lf.lr.srp.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',lf.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Now do the same for SS
#need new regression dataset
#use every seventh day

mean(lf.data$SS.mg.L,na.rm=T)+(sd(lf.data$SS.mg.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
#lf.data$REGR[lf.data$SS.mg.L>=637]=TRUE

#Trim down to just work with SS
lf.ss<-lf.data[,c(1,5,6,7,8)]
#remove NAs
lf.ss<-lf.ss[!is.na(lf.ss$SS.mg.L),]
lf.ss<-lf.ss[!is.na(lf.ss$Q.m3.sec),]

#use only positive numbers
lf.ss<-lf.ss[lf.ss$SS.mg.L>0,]
lf.ss<-lf.ss[lf.ss$Q.m3.sec>0,]

lf.ss<-lf.ss[-grep('1997-02-03',lf.ss$dateTime),]

#subset data for interpolation purposes
lf.reg.ss<-subset(lf.ss,REGR)
#lf.reg.ss<-lf.reg.ss[-grep('1997-02-11',lf.reg.ss$dateTime),]

#create a metadata description of the dataset and desired output
meta.ss<-metadata(constituent='SS.mg.L',flow='Q.m3.sec',dates='dateTime',conc.units='mg L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')
mods<-list()
lr.1<-loadReg2(loadReg(SS.mg.L~model(1),data=lf.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SS.mg.L~model(2),data=lf.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SS.mg.L~model(3),data=lf.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SS.mg.L~model(4),data=lf.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SS.mg.L~model(5),data=lf.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SS.mg.L~model(6),data=lf.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SS.mg.L~model(7),data=lf.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SS.mg.L~model(8),data=lf.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SS.mg.L~model(9),data=lf.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

lf.lr.ss.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',lf.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Marshalls Branch
setwd('~/Documents/Miami U/Acton load data/R')
mb<-read.csv('AllConcs9414_MB.csv') #four mile creek

#fix date so its in POSIX format
mb$MB.Date_time<-ymd_hm(mb$MB.Date_time)

#Start with four mile creek
#merge concentrations and solute data
mb.data<-merge(mb,q.hour,by.x='MB.Date_time',by.y='Date.Time')
#remove blank days

#remove NA columns and unused MB and LFM columns
mb.data<-mb.data[,-c(6,7,8,9,10,11,13,14,15)]

#rename columns
colnames(mb.data)<-c('dateTime','NH4.ug.L','NO3.ug.L','SRP.ug.L','SS.mg.L','Q.m3.sec')

#remove NAs
mb.data<-mb.data[!is.na(mb.data$dateTime),]

#use every seventh day
regr<-TRUE
for(i in 2:nrow(mb.data)){
	if(as.Date(mb.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(mb.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr[is.na(regr)]=FALSE
mb.data$REGR=regr

#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
#mb.data$REGR[mb.data$NH4.ug.L>344]=TRUE

#add year to data frame
mb.data$year<-format(ymd_hms(mb.data$dateTime),'%Y')

#Trim down to just work with NH4 first
mb.nh4<-mb.data[,c(1,2,6,7,8)]
#remove NAs
mb.nh4<-mb.nh4[!is.na(mb.nh4$NH4.ug.L),]
mb.nh4<-mb.nh4[!is.na(mb.nh4$Q.m3.sec),]

#use only positive numbers
mb.nh4<-mb.nh4[mb.nh4$NH4.ug.L>0,]
mb.nh4<-mb.nh4[mb.nh4$Q.m3.sec>0,]

#subset data for interpolation purposes
mb.reg.nh4<-subset(mb.nh4,REGR)
mb.discharge<-data.frame(dateTime=q.hour$Date.Time,Q.m3.sec=q.hour$MBQ.m3.sec)
mb.discharge<-mb.discharge[!is.na(mb.discharge$Q.m3.sec),]
mb.discharge$Q.m3.sec[mb.discharge$Q.m3.sec<=0]=0.0001
#add year to discharge
mb.discharge$year<-format(ymd_hms(mb.discharge$dateTime),'%Y')

#remove NAs
mb.reg.nh4<-mb.reg.nh4[!is.na(mb.reg.nh4$NH4.ug.L),]
mb.reg.nh4<-mb.reg.nh4[!is.na(mb.reg.nh4$Q.m3.sec),]
#remove anything with zeros
mb.reg.nh4<-mb.reg.nh4[mb.reg.nh4$Q.m3.sec>0,]
mb.reg.nh4<-mb.reg.nh4[mb.reg.nh4$NH4.ug.L>0,]

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='NH4.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')


mods<-list()
lr.1<-loadReg2(loadReg(NH4.ug.L~model(1),data=mb.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NH4.ug.L~model(2),data=mb.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NH4.ug.L~model(3),data=mb.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NH4.ug.L~model(4),data=mb.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NH4.ug.L~model(5),data=mb.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NH4.ug.L~model(6),data=mb.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NH4.ug.L~model(7),data=mb.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NH4.ug.L~model(8),data=mb.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NH4.ug.L~model(9),data=mb.reg.nh4,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

mb.lr.nh4.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',mb.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Now do the same for NO3
#need new regression dataset


#mean(mb.data$NO3.ug.L,na.rm=T)+(sd(mb.data$NO3.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
mb.data$REGR[mb.data$NO3.ug.L>=12350]=TRUE

#Trim down to just work with NO3 
mb.no3<-mb.data[,c(1,3,6,7,8)]
#remove NAs
mb.no3<-mb.no3[!is.na(mb.no3$NO3.ug.L),]
mb.no3<-mb.no3[!is.na(mb.no3$Q.m3.sec),]

#use only positive numbers
mb.no3<-mb.no3[mb.no3$NO3.ug.L>0,]
mb.no3<-mb.no3[mb.no3$Q.m3.sec>0,]

#subset data for interpolation purposes
mb.reg.no3<-subset(mb.no3,REGR)

#create a metadata description of the dataset and desired output
meta.no3<-metadata(constituent='NO3.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

mods<-list()
lr.1<-loadReg2(loadReg(NO3.ug.L~model(1),data=mb.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NO3.ug.L~model(2),data=mb.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NO3.ug.L~model(3),data=mb.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NO3.ug.L~model(4),data=mb.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NO3.ug.L~model(5),data=mb.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NO3.ug.L~model(6),data=mb.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NO3.ug.L~model(7),data=mb.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NO3.ug.L~model(8),data=mb.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NO3.ug.L~model(9),data=mb.reg.no3,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

mb.lr.no3.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',mb.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Now do the same for SRP

mean(mb.data$SRP.ug.L,na.rm=T)+(sd(mb.data$SRP.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
mb.data$REGR[mb.data$SRP.ug.L>=182]=TRUE

#Trim down to just work with NO3 
mb.srp<-mb.data[,c(1,4,6,7,8)]
#remove NAs
mb.srp<-mb.srp[!is.na(mb.srp$SRP.ug.L),]
mb.srp<-mb.srp[!is.na(mb.srp$Q.m3.sec),]

#use only positive numbers
mb.srp<-mb.srp[mb.srp$SRP.ug.L>0,]
mb.srp<-mb.srp[mb.srp$Q.m3.sec>0,]

#subset data for interpolation purposes
mb.reg.srp<-subset(mb.srp,REGR)

#create a metadata description of the dataset and desired output
meta.srp<-metadata(constituent='SRP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')
mods<-list()
lr.1<-loadReg2(loadReg(SRP.ug.L~model(1),data=mb.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SRP.ug.L~model(2),data=mb.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SRP.ug.L~model(3),data=mb.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SRP.ug.L~model(4),data=mb.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SRP.ug.L~model(5),data=mb.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SRP.ug.L~model(6),data=mb.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SRP.ug.L~model(7),data=mb.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SRP.ug.L~model(8),data=mb.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SRP.ug.L~model(9),data=mb.reg.srp,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

mb.lr.srp.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',mb.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#Now do the same for SS


mean(mb.data$SS.mg.L,na.rm=T)+(sd(mb.data$SS.mg.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
#mb.data$REGR[mb.data$SS.mg.L>=637]=TRUE

#Trim down to just work with SS
mb.ss<-mb.data[,c(1,5,6,7,8)]
#remove NAs
mb.ss<-mb.ss[!is.na(mb.ss$SS.mg.L),]
mb.ss<-mb.ss[!is.na(mb.ss$Q.m3.sec),]

#use only positive numbers
mb.ss<-mb.ss[mb.ss$SS.mg.L>0,]
mb.ss<-mb.ss[mb.ss$Q.m3.sec>0,]

#mb.ss<-mb.ss[-grep('1997-02-03',mb.ss$dateTime),]

#subset data for interpolation purposes
mb.reg.ss<-subset(mb.ss,REGR)
#mb.reg.ss<-mb.reg.ss[-grep('1997-02-11',mb.reg.ss$dateTime),]

#create a metadata description of the dataset and desired output
meta.ss<-metadata(constituent='SS.mg.L',flow='Q.m3.sec',dates='dateTime',conc.units='mg L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')
mods<-list()
lr.1<-loadReg2(loadReg(SS.mg.L~model(1),data=mb.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SS.mg.L~model(2),data=mb.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SS.mg.L~model(3),data=mb.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SS.mg.L~model(4),data=mb.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SS.mg.L~model(5),data=mb.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SS.mg.L~model(6),data=mb.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SS.mg.L~model(7),data=mb.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SS.mg.L~model(8),data=mb.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SS.mg.L~model(9),data=mb.reg.ss,flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))

mb.lr.ss.flux<-predictSolute(mods[[which(aics==min(aics))]],'flux',mb.discharge,se.pred=T,date=T)
getFittedModel(mods[[which(aics==min(aics))]])

#setwd('~/Documents/Miami U/loadflex/best model_data')
#save.image('allData_bestModelRegression.RData')