#Try to fit model using calibration data for just one year and sequentially build from there -- compare to the model using all 20 years as a comparison to see whic hworks best

rm(list=ls())
setwd('~')
source('.Rprofile')

library(knitr)
opts_chunk$set(echo=TRUE, message=FALSE)
library(loadflex)
library(rloadest)

#load Acton concentration data
setwd('~/Documents/Miami U/Acton load data/R')
mb<-read.csv('AllConcs9414_MB.csv') #Marshall's Branch

#load discharge data
q.daily<-read.csv('AllHDMQ9414_daily.csv')
q.hour<-read.csv('AllHDMQ9414_hourly.csv')
q.monthly<-read.csv('AllHDMQ9414_monthly.csv')

#fix date so its in POSIX format
mb$MB.Date_time<-ymd_hm(mb$MB.Date_time)
q.hour$Date.Time<-ymd_hm(q.hour$Date.Time)

mb<-mb[,c(1,2,3,4,5)]
q.hour<-q.hour[c(1,6)]

#Start with four mile creek
#merge concentrations and solute data
mb.data<-merge(mb,q.hour,by.x='MB.Date_time',by.y='Date.Time')
#remove blank days
#mb.data<-mb.data[mb.data$dateTime!='',]

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
mean(mb.data$NH4.ug.L,na.rm=T)+(2*sd(mb.data$NH4.ug.L,na.rm=T))
mb.data$REGR[mb.data$NH4.ug.L>=390]=TRUE

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
mb.discharge<-mb.discharge[mb.discharge$Q.m3.sec>0,]
#add year to discharge
mb.discharge$year<-format(ymd_hms(mb.discharge$dateTime),'%Y')

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='NH4.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

years<-unique(mb.data$year)

li.data<-c()
lm.data<-c()
lr.data<-c()
lc.data<-c()
summaryStats.lm<-list()
summaryStats.lr<-list()
summaryStats.lc<-list()
summaryStats.lc.resid<-list()
summaryStats.lc.CF<-c()
niter=1
for(i in seq(1,length(years),3)){
print(niter)
#fit four models: interpolation, linear, rloadest, and composite. 
mb.nh4_li<-loadInterp(interp.format='conc',interp.fun=rectangularInterpolation,data=mb.nh4[mb.nh4$year>=years[i] & mb.nh4$year<=years[i+2],],metadata=meta) #interpolation model
li.conc.pred<-predictSolute(mb.nh4_li,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
li.flux.pred<-predictSolute(mb.nh4_li,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
li.pred<-data.frame(date=li.conc.pred$date,conc.pred=li.conc.pred$fit,conc.se=li.conc.pred$se.pred,flux.pred=li.flux.pred$fit,flux.se=li.flux.pred$se.pred)
li.data<-rbind(li.data,li.pred)

mb.nh4_lm<-loadLm(formula=log(NH4.ug.L)~log(Q.m3.sec),pred.format='conc',data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],metadata=meta,retrans=exp) #linear model
lm.conc.pred<-predictSolute(mb.nh4_lm,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.flux.pred<-predictSolute(mb.nh4_lm,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.pred<-data.frame(date=lm.conc.pred$date,conc.pred=lm.conc.pred$fit,conc.se=lm.conc.pred$se.pred,flux.pred=lm.flux.pred$fit,flux.se=lm.flux.pred$se.pred)
lm.data<-rbind(lm.data,lm.pred)
summaryStats.lm[[niter]]<-summary(getFittedModel(mb.nh4_lm))

mods<-list()
lr.1<-loadReg2(loadReg(NH4.ug.L~model(1),data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NH4.ug.L~model(2),data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NH4.ug.L~model(3),data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NH4.ug.L~model(4),data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NH4.ug.L~model(5),data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NH4.ug.L~model(6),data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NH4.ug.L~model(7),data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NH4.ug.L~model(8),data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NH4.ug.L~model(9),data=mb.reg.nh4[mb.reg.nh4$year>=years[i] & mb.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))
lr.conc.pred<-predictSolute(mods[[which(aics==min(aics))]],'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.flux.pred<-predictSolute(mods[[which(aics==min(aics))]],'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.pred<-data.frame(date=lr.conc.pred$date,conc.pred=lr.conc.pred$fit,conc.se=lr.conc.pred$se.pred,flux.pred=lr.flux.pred$fit,flux.se=lr.flux.pred$se.pred)
lr.data<-rbind(lr.data,lr.pred)
summaryStats.lr[[niter]]<-getFittedModel(mods[[which(aics==min(aics))]])


mb.nh4_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=mb.nh4[mb.nh4$year>=years[i] & mb.nh4$year<=years[i+2],])
lc.conc.pred<-predictSolute(mb.nh4_lc,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(mb.nh4_lc,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
lc.data<-rbind(lc.data,lc.pred)
summaryStats.lc[[niter]]<-getFittedModel(mb.nh4_lc)
summaryStats.lc.resid[[niter]]<-getResiduals(mb.nh4_lc,'flux')
summaryStats.lc.CF[niter]<-getCorrectionFraction(mb.nh4_lc)

print(i)
niter=niter+1
}


#Now do the same for NO3
#need new regression dataset
#use every seventh day
regr<-TRUE
for(i in 2:nrow(mb.data)){
	if(as.Date(mb.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(mb.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr[is.na(regr)]=FALSE
mb.data$REGR=regr

mean(mb.data$NO3.ug.L,na.rm=T)+(sd(mb.data$NO3.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
mb.data$REGR[mb.data$NO3.ug.L>=15420]=TRUE

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

years<-unique(mb.data$year)

li.no3.data<-c()
lm.no3.data<-c()
lr.no3.data<-c()
lc.no3.data<-c()
summaryStats.no3.lm<-list()
summaryStats.no3.lr<-list()
summaryStats.no3.lc<-list()
summaryStats.no3.lc.resid<-list()
summaryStats.no3.lc.CF<-c()
niter=1
for(i in seq(1,length(years),3)){

#fit four models: interpolation, linear, rloadest, and composite. 
mb.no3_li<-loadInterp(interp.format='conc',interp.fun=rectangularInterpolation,data=mb.no3[mb.no3$year>=years[i] & mb.no3$year<=years[i+2],],metadata=meta.no3) #interpolation model
li.conc.pred<-predictSolute(mb.no3_li,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
li.flux.pred<-predictSolute(mb.no3_li,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
li.pred<-data.frame(date=li.conc.pred$date,conc.pred=li.conc.pred$fit,conc.se=li.conc.pred$se.pred,flux.pred=li.flux.pred$fit,flux.se=li.flux.pred$se.pred)
li.no3.data<-rbind(li.no3.data,li.pred)

mb.no3_lm<-loadLm(formula=log(NO3.ug.L)~log(Q.m3.sec),pred.format='conc',data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],metadata=meta.no3,retrans=exp) #linear model
lm.conc.pred<-predictSolute(mb.no3_lm,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.flux.pred<-predictSolute(mb.no3_lm,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.pred<-data.frame(date=lm.conc.pred$date,conc.pred=lm.conc.pred$fit,conc.se=lm.conc.pred$se.pred,flux.pred=lm.flux.pred$fit,flux.se=lm.flux.pred$se.pred)
lm.no3.data<-rbind(lm.no3.data,lm.pred)
summaryStats.no3.lm[[niter]]<-summary(getFittedModel(mb.no3_lm))

mods<-list()
lr.1<-loadReg2(loadReg(NO3.ug.L~model(1),data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NO3.ug.L~model(2),data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NO3.ug.L~model(3),data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NO3.ug.L~model(4),data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NO3.ug.L~model(5),data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NO3.ug.L~model(6),data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NO3.ug.L~model(7),data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NO3.ug.L~model(8),data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NO3.ug.L~model(9),data=mb.reg.no3[mb.reg.no3$year>=years[i] & mb.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))
lr.conc.pred<-predictSolute(mods[[which(aics==min(aics))]],'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.flux.pred<-predictSolute(mods[[which(aics==min(aics))]],'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.pred<-data.frame(date=lr.conc.pred$date,conc.pred=lr.conc.pred$fit,conc.se=lr.conc.pred$se.pred,flux.pred=lr.flux.pred$fit,flux.se=lr.flux.pred$se.pred)
lr.no3.data<-rbind(lr.no3.data,lr.pred)
summaryStats.no3.lr[[niter]]<-getFittedModel(mods[[which(aics==min(aics))]])

mb.no3_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=mb.no3[mb.no3$year>=years[i] & mb.no3$year<=years[i+2],])
lc.conc.pred<-predictSolute(mb.no3_lc,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(mb.no3_lc,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
lc.no3.data<-rbind(lc.no3.data,lc.pred)
summaryStats.no3.lc[[niter]]<-getFittedModel(mb.no3_lc)
summaryStats.no3.lc.resid[[niter]]<-getResiduals(mb.no3_lc,'conc')
summaryStats.no3.lc.CF[niter]<-getCorrectionFraction(mb.no3_lc)

niter=niter+1
print(i)
}

#Now do the same for SRP
#need new regression dataset
#use every seventh day
regr<-TRUE
for(i in 2:nrow(mb.data)){
	if(as.Date(mb.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(mb.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr[is.na(regr)]=FALSE
mb.data$REGR=regr

mean(mb.data$SRP.ug.L,na.rm=T)+(sd(mb.data$SRP.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
mb.data$REGR[mb.data$SRP.ug.L>=178]=TRUE

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

li.srp.data<-c()
lm.srp.data<-c()
lr.srp.data<-c()
lc.srp.data<-c()
summaryStats.srp.lm<-list()
summaryStats.srp.lr<-list()
summaryStats.srp.lc<-list()
summaryStats.srp.lc.resid<-list()
summaryStats.srp.lc.CF<-c()
niter=1
for(i in seq(1,length(years),3)){

#fit four models: interpolation, linear, rloadest, and composite. 
mb.srp_li<-loadInterp(interp.format='conc',interp.fun=rectangularInterpolation,data=mb.srp[mb.srp$year>=years[i] & mb.srp$year<=years[i+2],],metadata=meta.srp) #interpolation model
li.conc.pred<-predictSolute(mb.srp_li,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
li.flux.pred<-predictSolute(mb.srp_li,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
li.pred<-data.frame(date=li.conc.pred$date,conc.pred=li.conc.pred$fit,conc.se=li.conc.pred$se.pred,flux.pred=li.flux.pred$fit,flux.se=li.flux.pred$se.pred)
li.srp.data<-rbind(li.srp.data,li.pred)

mb.srp_lm<-loadLm(formula=log(SRP.ug.L)~log(Q.m3.sec),pred.format='conc',data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],metadata=meta.srp,retrans=exp) #linear model
lm.conc.pred<-predictSolute(mb.srp_lm,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.flux.pred<-predictSolute(mb.srp_lm,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.pred<-data.frame(date=lm.conc.pred$date,conc.pred=lm.conc.pred$fit,conc.se=lm.conc.pred$se.pred,flux.pred=lm.flux.pred$fit,flux.se=lm.flux.pred$se.pred)
lm.srp.data<-rbind(lm.srp.data,lm.pred)
summaryStats.srp.lm[[niter]]<-summary(getFittedModel(mb.srp_lm))

mods<-list()
lr.1<-loadReg2(loadReg(SRP.ug.L~model(1),data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SRP.ug.L~model(2),data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SRP.ug.L~model(3),data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SRP.ug.L~model(4),data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SRP.ug.L~model(5),data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SRP.ug.L~model(6),data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SRP.ug.L~model(7),data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SRP.ug.L~model(8),data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SRP.ug.L~model(9),data=mb.reg.srp[mb.reg.srp$year>=years[i] & mb.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))
lr.conc.pred<-predictSolute(mods[[which(aics==min(aics))]],'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.flux.pred<-predictSolute(mods[[which(aics==min(aics))]],'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.pred<-data.frame(date=lr.conc.pred$date,conc.pred=lr.conc.pred$fit,conc.se=lr.conc.pred$se.pred,flux.pred=lr.flux.pred$fit,flux.se=lr.flux.pred$se.pred)
lr.srp.data<-rbind(lr.srp.data,lr.pred)
summaryStats.srp.lr[[niter]]<-getFittedModel(mods[[which(aics==min(aics))]])


mb.srp_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=mb.srp[mb.srp$year>=years[i] & mb.srp$year<=years[i+2],])
lc.conc.pred<-predictSolute(mb.srp_lc,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(mb.srp_lc,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
lc.srp.data<-rbind(lc.srp.data,lc.pred)
summaryStats.srp.lc[[niter]]<-getFittedModel(mb.srp_lc)
summaryStats.srp.lc.resid[[niter]]<-getResiduals(mb.srp_lc,'conc')
summaryStats.srp.lc.CF[niter]<-getCorrectionFraction(mb.srp_lc)

niter=niter+1
print(i)
}

#Now do the same for SS
#need new regression dataset
#use every seventh day
regr<-TRUE
for(i in 2:nrow(mb.data)){
	if(as.Date(mb.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(mb.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
#regr<-c(regr,rep(FALSE,10))
regr[is.na(regr)]=FALSE
mb.data$REGR=regr

mean(mb.data$SS.mg.L,na.rm=T)+(sd(mb.data$SS.mg.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
mb.data$REGR[mb.data$SS.mg.L>=389]=TRUE

#Trim down to just work with NO3 
mb.ss<-mb.data[,c(1,5,6,7,8)]
#remove NAs
mb.ss<-mb.ss[!is.na(mb.ss$SS.mg.L),]
mb.ss<-mb.ss[!is.na(mb.ss$Q.m3.sec),]

#use only positive numbers
mb.ss<-mb.ss[mb.ss$SS.mg.L>0,]
mb.ss<-mb.ss[mb.ss$Q.m3.sec>0,]

mb.ss$SS.ug.L<-mb.ss$SS.mg.L*1000

#subset data for interpolation purposes
mb.reg.ss<-subset(mb.ss,REGR)

#create a metadata description of the dataset and desired output
meta.ss<-metadata(constituent='SS.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

li.ss.data<-c()
lm.ss.data<-c()
lr.ss.data<-c()
lc.ss.data<-c()
summaryStats.ss.lm<-list()
summaryStats.ss.lr<-list()
summaryStats.ss.lc<-list()
summaryStats.ss.lc.resid<-list()
summaryStats.ss.lc.CF<-c()
niter=1
for(i in seq(1,length(years),3)){

#fit four models: interpolation, linear, rloadest, and composite. 
mb.ss_li<-loadInterp(interp.format='conc',interp.fun=rectangularInterpolation,data=mb.ss[mb.ss$year>=years[i] & mb.ss$year<=years[i+2],],metadata=meta.ss) #interpolation model
li.conc.pred<-predictSolute(mb.ss_li,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
li.flux.pred<-predictSolute(mb.ss_li,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
li.pred<-data.frame(date=li.conc.pred$date,conc.pred=li.conc.pred$fit,conc.se=li.conc.pred$se.pred,flux.pred=li.flux.pred$fit,flux.se=li.flux.pred$se.pred)
li.ss.data<-rbind(li.ss.data,li.pred)

mb.ss_lm<-loadLm(formula=log(SS.mg.L)~log(Q.m3.sec),pred.format='conc',data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],metadata=meta.ss,retrans=exp) #linear model
lm.conc.pred<-predictSolute(mb.ss_lm,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.flux.pred<-predictSolute(mb.ss_lm,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.pred<-data.frame(date=lm.conc.pred$date,conc.pred=lm.conc.pred$fit,conc.se=lm.conc.pred$se.pred,flux.pred=lm.flux.pred$fit,flux.se=lm.flux.pred$se.pred)
lm.ss.data<-rbind(lm.ss.data,lm.pred)
summaryStats.ss.lm[[niter]]<-summary(getFittedModel(mb.ss_lm))

mods<-list()
lr.1<-loadReg2(loadReg(SS.mg.L~model(1),data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SS.mg.L~model(2),data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SS.mg.L~model(3),data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SS.mg.L~model(4),data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SS.mg.L~model(5),data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SS.mg.L~model(6),data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SS.mg.L~model(7),data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SS.mg.L~model(8),data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SS.mg.L~model(9),data=mb.reg.ss[mb.reg.ss$year>=years[i] & mb.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))
lr.conc.pred<-predictSolute(mods[[which(aics==min(aics))]],'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.flux.pred<-predictSolute(mods[[which(aics==min(aics))]],'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.pred<-data.frame(date=lr.conc.pred$date,conc.pred=lr.conc.pred$fit,conc.se=lr.conc.pred$se.pred,flux.pred=lr.flux.pred$fit,flux.se=lr.flux.pred$se.pred)
lr.ss.data<-rbind(lr.ss.data,lr.pred)
summaryStats.ss.lr[[niter]]<-getFittedModel(mods[[which(aics==min(aics))]])

mb.ss_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=mb.ss[mb.ss$year>=years[i] & mb.ss$year<=years[i+2],])
lc.conc.pred<-predictSolute(mb.ss_lc,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(mb.ss_lc,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
lc.ss.data<-rbind(lc.ss.data,lc.pred)
summaryStats.ss.lc[[niter]]<-getFittedModel(mb.ss_lc)
summaryStats.ss.lc.resid[[niter]]<-getResiduals(mb.ss_lc,'conc')
summaryStats.ss.lc.CF[niter]<-getCorrectionFraction(mb.ss_lc)

niter=niter+1
print(i)
print(mods[[which(aics==min(aics))]])
}

setwd('~/Documents/Miami U/loadflex/best model_data')
save.image('TRIMMED.MB.loads.RData')