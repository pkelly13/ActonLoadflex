#Try to fit model using calibration data for just one year and sequentially build from there -- compare to the model using all 20 years as a comparison to see whic hworks best

library(knitr)
opts_chunk$set(echo=TRUE, message=FALSE)
library(loadflex)
library(rloadest)

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

fm<-fm[,c(1,2,3,4,5)]
q.hour<-q.hour[c(1,4)]

#Start with four mile creek
#merge concentrations and solute data
fm.data<-merge(fm,q.hour,by.x='FM.Date_time',by.y='Date.Time')
#remove blank days
#fm.data<-fm.data[fm.data$dateTime!='',]

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
mean(fm.data$NH4.ug.L,na.rm=T)+(2*sd(fm.data$NH4.ug.L,na.rm=T))
fm.data$REGR[fm.data$NH4.ug.L>=390]=TRUE

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
fm.discharge<-fm.discharge[fm.discharge$Q.m3.sec>0,]
#add year to discharge
fm.discharge$year<-format(ymd_hms(fm.discharge$dateTime),'%Y')

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='NH4.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

years<-unique(fm.data$year)

li.data<-c()
lm.data<-c()
lr.data<-c()
lc.data<-c()
fm.nh4.lr.ConcResid<-list()
fm.nh4.lr.FluxResid<-list()
summaryStats.lm<-list()
summaryStats.lr<-list()
summaryStats.lc<-list()
summaryStats.lc.resid<-list()
summaryStats.lc.CF<-c()
niter=1
for(i in seq(1,length(years),3)){
print(niter)
#fit four models: interpolation, linear, rloadest, and composite. 
fm.nh4_li<-loadInterp(interp.format='conc',interp.fun=rectangularInterpolation,data=fm.nh4[fm.nh4$year>=years[i] & fm.nh4$year<=years[i+2],],metadata=meta) #interpolation model
li.conc.pred<-predictSolute(fm.nh4_li,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
li.flux.pred<-predictSolute(fm.nh4_li,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
li.pred<-data.frame(date=li.conc.pred$date,conc.pred=li.conc.pred$fit,conc.se=li.conc.pred$se.pred,flux.pred=li.flux.pred$fit,flux.se=li.flux.pred$se.pred)
li.data<-rbind(li.data,li.pred)

fm.nh4_lm<-loadLm(formula=log(NH4.ug.L)~log(Q.m3.sec),pred.format='conc',data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],metadata=meta,retrans=exp) #linear model
lm.conc.pred<-predictSolute(fm.nh4_lm,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.flux.pred<-predictSolute(fm.nh4_lm,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.pred<-data.frame(date=lm.conc.pred$date,conc.pred=lm.conc.pred$fit,conc.se=lm.conc.pred$se.pred,flux.pred=lm.flux.pred$fit,flux.se=lm.flux.pred$se.pred)
lm.data<-rbind(lm.data,lm.pred)
summaryStats.lm[[niter]]<-summary(getFittedModel(fm.nh4_lm))

mods<-list()
lr.1<-loadReg2(loadReg(NH4.ug.L~model(1),data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NH4.ug.L~model(2),data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NH4.ug.L~model(3),data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NH4.ug.L~model(4),data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NH4.ug.L~model(5),data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NH4.ug.L~model(6),data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NH4.ug.L~model(7),data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NH4.ug.L~model(8),data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NH4.ug.L~model(9),data=fm.reg.nh4[fm.reg.nh4$year>=years[i] & fm.reg.nh4$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))
lr.conc.pred<-predictSolute(mods[[which(aics==min(aics))]],'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.flux.pred<-predictSolute(mods[[which(aics==min(aics))]],'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.pred<-data.frame(date=lr.conc.pred$date,conc.pred=lr.conc.pred$fit,conc.se=lr.conc.pred$se.pred,flux.pred=lr.flux.pred$fit,flux.se=lr.flux.pred$se.pred)
lr.data<-rbind(lr.data,lr.pred)
summaryStats.lr[[niter]]<-getFittedModel(mods[[which(aics==min(aics))]])
fm.nh4.lr.ConcResid[[niter]]<-getResiduals(mods[[which(aics==min(aics))]],'conc','absolute')
fm.nh4.lr.FluxResid[[niter]]<-getResiduals(mods[[which(aics==min(aics))]],'flux','absolute')

#fm.nh4_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=fm.nh4[fm.nh4$year>=years[i] & fm.nh4$year<=years[i+2],])
#lc.conc.pred<-predictSolute(fm.nh4_lc,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
#lc.flux.pred<-predictSolute(fm.nh4_lc,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
#lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
#lc.data<-rbind(lc.data,lc.pred)
#summaryStats.lc[[niter]]<-getFittedModel(fm.nh4_lc)
#summaryStats.lc.resid[[niter]]<-getResiduals(fm.nh4_lc,'flux')
#summaryStats.lc.CF[niter]<-getCorrectionFraction(fm.nh4_lc)

print(i)
niter=niter+1
}


#Now do the same for NO3
#need new regression dataset
#use every seventh day
regr<-TRUE
for(i in 2:nrow(fm.data)){
	if(as.Date(fm.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(fm.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr[is.na(regr)]=FALSE
fm.data$REGR=regr

mean(fm.data$NO3.ug.L,na.rm=T)+(sd(fm.data$NO3.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
fm.data$REGR[fm.data$NO3.ug.L>=15420]=TRUE

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

years<-unique(fm.data$year)

li.no3.data<-c()
lm.no3.data<-c()
lr.no3.data<-c()
lc.no3.data<-c()
fm.no3.lr.ConcResid<-list()
fm.no3.lr.FluxResid<-list()
summaryStats.no3.lm<-list()
summaryStats.no3.lr<-list()
summaryStats.no3.lc<-list()
summaryStats.no3.lc.resid<-list()
summaryStats.no3.lc.CF<-c()
niter=1
for(i in seq(1,length(years),3)){

#fit four models: interpolation, linear, rloadest, and composite. 
fm.no3_li<-loadInterp(interp.format='conc',interp.fun=rectangularInterpolation,data=fm.no3[fm.no3$year>=years[i] & fm.no3$year<=years[i+2],],metadata=meta.no3) #interpolation model
li.conc.pred<-predictSolute(fm.no3_li,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
li.flux.pred<-predictSolute(fm.no3_li,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
li.pred<-data.frame(date=li.conc.pred$date,conc.pred=li.conc.pred$fit,conc.se=li.conc.pred$se.pred,flux.pred=li.flux.pred$fit,flux.se=li.flux.pred$se.pred)
li.no3.data<-rbind(li.no3.data,li.pred)

fm.no3_lm<-loadLm(formula=log(NO3.ug.L)~log(Q.m3.sec),pred.format='conc',data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],metadata=meta.no3,retrans=exp) #linear model
lm.conc.pred<-predictSolute(fm.no3_lm,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.flux.pred<-predictSolute(fm.no3_lm,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.pred<-data.frame(date=lm.conc.pred$date,conc.pred=lm.conc.pred$fit,conc.se=lm.conc.pred$se.pred,flux.pred=lm.flux.pred$fit,flux.se=lm.flux.pred$se.pred)
lm.no3.data<-rbind(lm.no3.data,lm.pred)
summaryStats.no3.lm[[niter]]<-summary(getFittedModel(fm.no3_lm))

mods<-list()
lr.1<-loadReg2(loadReg(NO3.ug.L~model(1),data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(NO3.ug.L~model(2),data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(NO3.ug.L~model(3),data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(NO3.ug.L~model(4),data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(NO3.ug.L~model(5),data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(NO3.ug.L~model(6),data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(NO3.ug.L~model(7),data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(NO3.ug.L~model(8),data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(NO3.ug.L~model(9),data=fm.reg.no3[fm.reg.no3$year>=years[i] & fm.reg.no3$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))
lr.conc.pred<-predictSolute(mods[[which(aics==min(aics))]],'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.flux.pred<-predictSolute(mods[[which(aics==min(aics))]],'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.pred<-data.frame(date=lr.conc.pred$date,conc.pred=lr.conc.pred$fit,conc.se=lr.conc.pred$se.pred,flux.pred=lr.flux.pred$fit,flux.se=lr.flux.pred$se.pred)
lr.no3.data<-rbind(lr.no3.data,lr.pred)
summaryStats.no3.lr[[niter]]<-getFittedModel(mods[[which(aics==min(aics))]])
fm.no3.lr.ConcResid[[niter]]<-getResiduals(mods[[which(aics==min(aics))]],'conc','absolute')
fm.no3.lr.FluxResid[[niter]]<-getResiduals(mods[[which(aics==min(aics))]],'flux','absolute')

#fm.no3_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=fm.no3[fm.no3$year>=years[i] & fm.no3$year<=years[i+2],])
#lc.conc.pred<-predictSolute(fm.no3_lc,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
#lc.flux.pred<-predictSolute(fm.no3_lc,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
#lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
#lc.no3.data<-rbind(lc.no3.data,lc.pred)
#summaryStats.no3.lc[[niter]]<-getFittedModel(fm.no3_lc)
#summaryStats.no3.lc.resid[[niter]]<-getResiduals(fm.no3_lc,'conc')
#summaryStats.no3.lc.CF[niter]<-getCorrectionFraction(fm.no3_lc)

niter=niter+1
print(i)
}

#Now do the same for SRP
#need new regression dataset
#use every seventh day
regr<-TRUE
for(i in 2:nrow(fm.data)){
	if(as.Date(fm.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(fm.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr[is.na(regr)]=FALSE
fm.data$REGR=regr

mean(fm.data$SRP.ug.L,na.rm=T)+(sd(fm.data$SRP.ug.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
fm.data$REGR[fm.data$SRP.ug.L>=178]=TRUE

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

li.srp.data<-c()
lm.srp.data<-c()
lr.srp.data<-c()
lc.srp.data<-c()
fm.srp.lr.ConcResid<-list()
fm.srp.lr.FluxResid<-list()
summaryStats.srp.lm<-list()
summaryStats.srp.lr<-list()
summaryStats.srp.lc<-list()
summaryStats.srp.lc.resid<-list()
summaryStats.srp.lc.CF<-c()
niter=1
for(i in seq(1,length(years),3)){

#fit four models: interpolation, linear, rloadest, and composite. 
fm.srp_li<-loadInterp(interp.format='conc',interp.fun=rectangularInterpolation,data=fm.srp[fm.srp$year>=years[i] & fm.srp$year<=years[i+2],],metadata=meta.srp) #interpolation model
li.conc.pred<-predictSolute(fm.srp_li,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
li.flux.pred<-predictSolute(fm.srp_li,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
li.pred<-data.frame(date=li.conc.pred$date,conc.pred=li.conc.pred$fit,conc.se=li.conc.pred$se.pred,flux.pred=li.flux.pred$fit,flux.se=li.flux.pred$se.pred)
li.srp.data<-rbind(li.srp.data,li.pred)

fm.srp_lm<-loadLm(formula=log(SRP.ug.L)~log(Q.m3.sec),pred.format='conc',data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],metadata=meta.srp,retrans=exp) #linear model
lm.conc.pred<-predictSolute(fm.srp_lm,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.flux.pred<-predictSolute(fm.srp_lm,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.pred<-data.frame(date=lm.conc.pred$date,conc.pred=lm.conc.pred$fit,conc.se=lm.conc.pred$se.pred,flux.pred=lm.flux.pred$fit,flux.se=lm.flux.pred$se.pred)
lm.srp.data<-rbind(lm.srp.data,lm.pred)
summaryStats.srp.lm[[niter]]<-summary(getFittedModel(fm.srp_lm))

mods<-list()
lr.1<-loadReg2(loadReg(SRP.ug.L~model(1),data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SRP.ug.L~model(2),data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SRP.ug.L~model(3),data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SRP.ug.L~model(4),data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SRP.ug.L~model(5),data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SRP.ug.L~model(6),data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SRP.ug.L~model(7),data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SRP.ug.L~model(8),data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SRP.ug.L~model(9),data=fm.reg.srp[fm.reg.srp$year>=years[i] & fm.reg.srp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))
lr.conc.pred<-predictSolute(mods[[which(aics==min(aics))]],'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.flux.pred<-predictSolute(mods[[which(aics==min(aics))]],'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.pred<-data.frame(date=lr.conc.pred$date,conc.pred=lr.conc.pred$fit,conc.se=lr.conc.pred$se.pred,flux.pred=lr.flux.pred$fit,flux.se=lr.flux.pred$se.pred)
lr.srp.data<-rbind(lr.srp.data,lr.pred)
summaryStats.srp.lr[[niter]]<-getFittedModel(mods[[which(aics==min(aics))]])
fm.srp.lr.ConcResid[[niter]]<-getResiduals(mods[[which(aics==min(aics))]],'conc','absolute')
fm.srp.lr.FluxResid[niter]<-getResiduals(mods[[which(aics==min(aics))]],'flux','absolute')


#fm.srp_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=fm.srp[fm.srp$year>=years[i] & fm.srp$year<=years[i+2],])
#lc.conc.pred<-predictSolute(fm.srp_lc,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
#lc.flux.pred<-predictSolute(fm.srp_lc,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
#lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
#lc.srp.data<-rbind(lc.srp.data,lc.pred)
#summaryStats.srp.lc[[niter]]<-getFittedModel(fm.srp_lc)
#summaryStats.srp.lc.resid[[niter]]<-getResiduals(fm.srp_lc,'conc')
#summaryStats.srp.lc.CF[niter]<-getCorrectionFraction(fm.srp_lc)

niter=niter+1
print(i)
}

#Now do the same for SS
#need new regression dataset
#use every seventh day
regr<-TRUE
for(i in 2:nrow(fm.data)){
	if(as.Date(fm.data$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(fm.data$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
#regr<-c(regr,rep(FALSE,10))
regr[is.na(regr)]=FALSE
fm.data$REGR=regr

mean(fm.data$SS.mg.L,na.rm=T)+(sd(fm.data$SS.mg.L,na.rm=T)*2)
#make REGR where NO3 is greater than 2 sds of the mean also TRUE to make sure high concentrations are included in the regression
fm.data$REGR[fm.data$SS.mg.L>=389]=TRUE

#Trim down to just work with NO3 
fm.ss<-fm.data[,c(1,5,6,7,8)]
#remove NAs
fm.ss<-fm.ss[!is.na(fm.ss$SS.mg.L),]
fm.ss<-fm.ss[!is.na(fm.ss$Q.m3.sec),]

#use only positive numbers
fm.ss<-fm.ss[fm.ss$SS.mg.L>0,]
fm.ss<-fm.ss[fm.ss$Q.m3.sec>0,]

fm.ss$SS.ug.L<-fm.ss$SS.mg.L*1000

#subset data for interpolation purposes
fm.reg.ss<-subset(fm.ss,REGR)

#create a metadata description of the dataset and desired output
meta.ss<-metadata(constituent='SS.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

li.ss.data<-c()
lm.ss.data<-c()
lr.ss.data<-c()
lc.ss.data<-c()
fm.ss.lr.ConcResid<-list()
fm.ss.lr.FluxResid<-list()
summaryStats.ss.lm<-list()
summaryStats.ss.lr<-list()
summaryStats.ss.lc<-list()
summaryStats.ss.lc.resid<-list()
summaryStats.ss.lc.CF<-c()
niter=1
for(i in seq(1,length(years),3)){

#fit four models: interpolation, linear, rloadest, and composite. 
fm.ss_li<-loadInterp(interp.format='conc',interp.fun=rectangularInterpolation,data=fm.ss[fm.ss$year>=years[i] & fm.ss$year<=years[i+2],],metadata=meta.ss) #interpolation model
li.conc.pred<-predictSolute(fm.ss_li,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
li.flux.pred<-predictSolute(fm.ss_li,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
li.pred<-data.frame(date=li.conc.pred$date,conc.pred=li.conc.pred$fit,conc.se=li.conc.pred$se.pred,flux.pred=li.flux.pred$fit,flux.se=li.flux.pred$se.pred)
li.ss.data<-rbind(li.ss.data,li.pred)

fm.ss_lm<-loadLm(formula=log(SS.mg.L)~log(Q.m3.sec),pred.format='conc',data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],metadata=meta.ss,retrans=exp) #linear model
lm.conc.pred<-predictSolute(fm.ss_lm,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.flux.pred<-predictSolute(fm.ss_lm,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lm.pred<-data.frame(date=lm.conc.pred$date,conc.pred=lm.conc.pred$fit,conc.se=lm.conc.pred$se.pred,flux.pred=lm.flux.pred$fit,flux.se=lm.flux.pred$se.pred)
lm.ss.data<-rbind(lm.ss.data,lm.pred)
summaryStats.ss.lm[[niter]]<-summary(getFittedModel(fm.ss_lm))

mods<-list()
lr.1<-loadReg2(loadReg(SS.mg.L~model(1),data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(SS.mg.L~model(2),data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(SS.mg.L~model(3),data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(SS.mg.L~model(4),data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(SS.mg.L~model(5),data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(SS.mg.L~model(6),data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(SS.mg.L~model(7),data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(SS.mg.L~model(8),data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(SS.mg.L~model(9),data=fm.reg.ss[fm.reg.ss$year>=years[i] & fm.reg.ss$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='mg/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))
lr.conc.pred<-predictSolute(mods[[which(aics==min(aics))]],'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.flux.pred<-predictSolute(mods[[which(aics==min(aics))]],'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lr.pred<-data.frame(date=lr.conc.pred$date,conc.pred=lr.conc.pred$fit,conc.se=lr.conc.pred$se.pred,flux.pred=lr.flux.pred$fit,flux.se=lr.flux.pred$se.pred)
lr.ss.data<-rbind(lr.ss.data,lr.pred)
summaryStats.ss.lr[[niter]]<-getFittedModel(mods[[which(aics==min(aics))]])
fm.ss.lr.ConcResid[[niter]]<-getResiduals(mods[[which(aics==min(aics))]],'conc','absolute')
fm.ss.lr.FluxResid[[niter]]<-getResiduals(mods[[which(aics==min(aics))]],'flux','absolute')

#fm.ss_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=fm.ss[fm.ss$year>=years[i] & fm.ss$year<=years[i+2],])
#lc.conc.pred<-predictSolute(fm.ss_lc,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
#lc.flux.pred<-predictSolute(fm.ss_lc,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
#lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
#lc.ss.data<-rbind(lc.ss.data,lc.pred)
#summaryStats.ss.lc[[niter]]<-getFittedModel(fm.ss_lc)
#summaryStats.ss.lc.resid[[niter]]<-getResiduals(fm.ss_lc,'conc')
#summaryStats.ss.lc.CF[niter]<-getCorrectionFraction(fm.ss_lc)

niter=niter+1
print(i)
print(mods[[which(aics==min(aics))]])
}

setwd('~/Documents/Miami U/loadflex/best model_data')
save.image('Reg.Residuals.RData')