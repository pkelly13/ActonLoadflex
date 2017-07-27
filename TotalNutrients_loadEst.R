#Calculate TDN, TDP, PP, and PN for Mike using NO3, SRP, and SS data - this requires converting NO3, SRP, and SS into those constituent concentrations using established linear regressions. These are stream specific.

#Start with Four Mile Creek
#load concentration data data
setwd('~/Documents/Miami U/Acton load data/R')

#Concentration data
fm.data<-read.csv('AllConcs9414.csv')

#Need dischareg data as well
q<-read.csv('AllHDMQ9414_hourly.csv')

#fix dates to POSIX format
fm.data$FM.Date_time<-ymd_hm(fm.data$FM.Date_time)
q$Date.Time<-ymd_hm(q$Date.Time)

fm.discharge<-data.frame(dateTime=q$Date.Time,Q.m3.sec=q$FMQ.m3.sec)
fm.discharge<-fm.discharge[!is.na(fm.discharge$Q.m3.sec),]
fm.discharge<-fm.discharge[fm.discharge$Q.m3.sec>0,]
#add year to discharge
fm.discharge$year<-format(ymd_hms(fm.discharge$dateTime),'%Y')
fm.discharge$dateTime<-ymd_hms(fm.discharge$dateTime)


fm.q<-q[,c(1,4)]
colnames(fm.q)[1]<-'dateTime'

#get rid of NA columns
fm.data<-fm.data[,1:5]

#make different data frames for TDN, TDP, Particulate N, and Particulate P
fm.tdn<-fm.data$FM.NH4.ug.L*1.004+885.74 #R2  = 0.821
fm.tdn<-data.frame(dateTime=fm.data$FM.Date_time,TDN.ug.L=fm.tdn)
fm.tdn<-merge(fm.tdn,fm.q,by='dateTime')
colnames(fm.tdn)<-c('dateTime','TDN.ug.L','Q.m3.sec')
fm.tdn<-fm.tdn[!is.na(fm.tdn$dateTime),]
fm.tdn<-fm.tdn[fm.tdn$Q.m3.sec>0,]

fm.tdp<-fm.data$FM.SRP.ug.L*1.170+10.385 #R2 = 0.816
fm.tdp<-data.frame(dateTime=fm.data$FM.Date_time,TDP.ug.L=fm.tdp)
fm.tdp<-merge(fm.tdp,fm.q,by='dateTime')
fm.tdp<-fm.tdp[!is.na(fm.tdp$dateTime),]
colnames(fm.tdp)<-c('dateTime','TDP.ug.L','Q.m3.sec')
fm.tdp<-fm.tdp[!is.na(fm.tdp$dateTime),]
fm.tdp<-fm.tdp[fm.tdp$Q.m3.sec>0,]

fm.pn<-fm.data$FM.SS.mg.L*2.618+115.82 #R2 = 0.924
fm.pn<-data.frame(dateTime=fm.data$FM.Date_time,PN.ug.L=fm.pn)
fm.pn<-merge(fm.pn,fm.q,by='dateTime')
fm.pn<-fm.pn[!is.na(fm.pn$dateTime),]
colnames(fm.pn)<-c('dateTime','PN.ug.L','Q.m3.sec')
fm.pn<-fm.pn[!is.na(fm.pn$dateTime),]
fm.pn<-fm.pn[fm.pn$Q.m3.sec>0,]

fm.pp<-fm.data$FM.SS.mg.L*0.279+58.683 #R2 = 0.673
fm.pp<-data.frame(dateTime=fm.data$FM.Date_time,PP.ug.L=fm.pp)
fm.pp<-merge(fm.pp,fm.q,by='dateTime')
fm.pp<-fm.pp[!is.na(fm.pp$dateTime),]
colnames(fm.pp)<-c('dateTime','PP.ug.L','Q.m3.sec')
fm.pp<-fm.pp[!is.na(fm.pp$dateTime),]
fm.pp<-fm.pp[fm.pp$Q.m3.sec>0,]

#Start with TDN for Four Mile Creek
#use the composite method
#need to first trim the dataset to only use every seven days
regr<-TRUE
for(i in 2:nrow(fm.tdn)){
	if(as.Date(fm.tdn$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(fm.tdn$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr[is.na(regr)]=FALSE
fm.tdn$REGR=regr
fm.tdp$REGR<-regr
fm.pn$REGR<-regr
fm.pp$REGR<-regr

#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
fm.tdn$REGR[fm.tdn$TDN.ug.L>=(mean(fm.tdn$TDN.ug.L,na.rm=T)+(2*sd(fm.tdn$TDN.ug.L,na.rm=T)))]=TRUE

#add year to data frame
fm.tdn$year<-format(ymd_hms(fm.tdn$dateTime),'%Y')

#remove NAs
fm.tdn<-fm.tdn[!is.na(fm.tdn$TDN.ug.L),]
fm.tdn<-fm.tdn[!is.na(fm.tdn$Q.m3.sec),]

#subset data for interpolation purposes
fm.reg.tdn<-subset(fm.tdn,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='TDN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

years<-unique(fm.tdn$year)


fm.tdn.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(TDN.ug.L~model(1),data=fm.reg.tdn[fm.reg.tdn$year>=years[i] & fm.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(TDN.ug.L~model(2),data=fm.reg.tdn[fm.reg.tdn$year>=years[i] & fm.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(TDN.ug.L~model(3),data=fm.reg.tdn[fm.reg.tdn$year>=years[i] & fm.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(TDN.ug.L~model(4),data=fm.reg.tdn[fm.reg.tdn$year>=years[i] & fm.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(TDN.ug.L~model(5),data=fm.reg.tdn[fm.reg.tdn$year>=years[i] & fm.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(TDN.ug.L~model(6),data=fm.reg.tdn[fm.reg.tdn$year>=years[i] & fm.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(TDN.ug.L~model(7),data=fm.reg.tdn[fm.reg.tdn$year>=years[i] & fm.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(TDN.ug.L~model(8),data=fm.reg.tdn[fm.reg.tdn$year>=years[i] & fm.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(TDN.ug.L~model(9),data=fm.reg.tdn[fm.reg.tdn$year>=years[i] & fm.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


fm.tdn_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=fm.tdn[fm.tdn$year>=years[i] & fm.tdn$year<=years[i+2],])
lc.conc.pred<-predictSolute(fm.tdn_lc,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(fm.tdn_lc,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
fm.tdn.comp<-rbind(fm.tdn.comp,lc.pred)

print(i)
}

#Now do TDP____________________________________
#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
fm.tdp$REGR[fm.tdp$TDP.ug.L>=(mean(fm.tdp$TDP.ug.L,na.rm=T)+(2*sd(fm.tdp$TDP.ug.L,na.rm=T)))]=TRUE

#add year to data frame
fm.tdp$year<-format(ymd_hms(fm.tdp$dateTime),'%Y')

#remove NAs
fm.tdp<-fm.tdp[!is.na(fm.tdp$TDP.ug.L),]
fm.tdp<-fm.tdp[!is.na(fm.tdp$Q.m3.sec),]

#subset data for interpolation purposes
fm.reg.tdp<-subset(fm.tdp,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='TDP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')


fm.tdp.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(TDP.ug.L~model(1),data=fm.reg.tdp[fm.reg.tdp$year>=years[i] & fm.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(TDP.ug.L~model(2),data=fm.reg.tdp[fm.reg.tdp$year>=years[i] & fm.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(TDP.ug.L~model(3),data=fm.reg.tdp[fm.reg.tdp$year>=years[i] & fm.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(TDP.ug.L~model(4),data=fm.reg.tdp[fm.reg.tdp$year>=years[i] & fm.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(TDP.ug.L~model(5),data=fm.reg.tdp[fm.reg.tdp$year>=years[i] & fm.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(TDP.ug.L~model(6),data=fm.reg.tdp[fm.reg.tdp$year>=years[i] & fm.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(TDP.ug.L~model(7),data=fm.reg.tdp[fm.reg.tdp$year>=years[i] & fm.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(TDP.ug.L~model(8),data=fm.reg.tdp[fm.reg.tdp$year>=years[i] & fm.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(TDP.ug.L~model(9),data=fm.reg.tdp[fm.reg.tdp$year>=years[i] & fm.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


fm.tdp_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=fm.tdp[fm.tdp$year>=years[i] & fm.tdp$year<=years[i+2],])
lc.conc.pred<-predictSolute(fm.tdp_lc,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(fm.tdp_lc,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
fm.tdp.comp<-rbind(fm.tdp.comp,lc.pred)

print(i)
}


#Now do PN____________________________
#Now do PN____________________________________
#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
fm.pn$REGR[fm.pn$PN.ug.L>=(mean(fm.pn$PN.ug.L,na.rm=T)+(2*sd(fm.pn$PN.ug.L,na.rm=T)))]=TRUE

#add year to data frame
fm.pn$year<-format(ymd_hms(fm.pn$dateTime),'%Y')

#remove NAs
fm.pn<-fm.pn[!is.na(fm.pn$PN.ug.L),]
fm.pn<-fm.pn[!is.na(fm.pn$Q.m3.sec),]

#subset data for interpolation purposes
fm.reg.pn<-subset(fm.pn,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='PN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')


fm.pn.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(PN.ug.L~model(1),data=fm.reg.pn[fm.reg.pn$year>=years[i] & fm.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(PN.ug.L~model(2),data=fm.reg.pn[fm.reg.pn$year>=years[i] & fm.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(PN.ug.L~model(3),data=fm.reg.pn[fm.reg.pn$year>=years[i] & fm.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(PN.ug.L~model(4),data=fm.reg.pn[fm.reg.pn$year>=years[i] & fm.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(PN.ug.L~model(5),data=fm.reg.pn[fm.reg.pn$year>=years[i] & fm.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(PN.ug.L~model(6),data=fm.reg.pn[fm.reg.pn$year>=years[i] & fm.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(PN.ug.L~model(7),data=fm.reg.pn[fm.reg.pn$year>=years[i] & fm.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(PN.ug.L~model(8),data=fm.reg.pn[fm.reg.pn$year>=years[i] & fm.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(PN.ug.L~model(9),data=fm.reg.pn[fm.reg.pn$year>=years[i] & fm.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


fm.pn_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=fm.pn[fm.pn$year>=years[i] & fm.pn$year<=years[i+2],])
lc.conc.pred<-predictSolute(fm.pn_lc,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(fm.pn_lc,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
fm.pn.comp<-rbind(fm.pn.comp,lc.pred)

print(i)
}

#NOW do PP____________________________
#Now do PP____________________________________
#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
fm.pp$REGR[fm.pp$PP.ug.L>=(mean(fm.pp$PP.ug.L,na.rm=T)+(2*sd(fm.pp$PP.ug.L,na.rm=T)))]=TRUE

#add year to data frame
fm.pp$year<-format(ymd_hms(fm.pp$dateTime),'%Y')

#remove NAs
fm.pp<-fm.pp[!is.na(fm.pp$PP.ug.L),]
fm.pp<-fm.pp[!is.na(fm.pp$Q.m3.sec),]

#subset data for interpolation purposes
fm.reg.pp<-subset(fm.pp,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='PP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')


fm.pp.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(PP.ug.L~model(1),data=fm.reg.pp[fm.reg.pp$year>=years[i] & fm.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(PP.ug.L~model(2),data=fm.reg.pp[fm.reg.pp$year>=years[i] & fm.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(PP.ug.L~model(3),data=fm.reg.pp[fm.reg.pp$year>=years[i] & fm.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(PP.ug.L~model(4),data=fm.reg.pp[fm.reg.pp$year>=years[i] & fm.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(PP.ug.L~model(5),data=fm.reg.pp[fm.reg.pp$year>=years[i] & fm.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(PP.ug.L~model(6),data=fm.reg.pp[fm.reg.pp$year>=years[i] & fm.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(PP.ug.L~model(7),data=fm.reg.pp[fm.reg.pp$year>=years[i] & fm.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(PP.ug.L~model(8),data=fm.reg.pp[fm.reg.pp$year>=years[i] & fm.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(PP.ug.L~model(9),data=fm.reg.pp[fm.reg.pp$year>=years[i] & fm.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


fm.pp_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=fm.pp[fm.pp$year>=years[i] & fm.pp$year<=years[i+2],])
lc.conc.pred<-predictSolute(fm.pp_lc,'conc',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(fm.pp_lc,'flux',fm.discharge[fm.discharge$year>=years[i] & fm.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
fm.pp.comp<-rbind(fm.pp.comp,lc.pred)

print(i)
}

#Calculate TDN, TDP, PP, and PN for Mike using NO3, SRP, and SS data - this requires converting NO3, SRP, and SS into those constituent concentrations using established linear regressions. These are stream specific.

#Little Four Mile Creek
#load concentration data data
setwd('~/Documents/Miami U/Acton load data/R')

#Concentration data
lf.data<-read.csv('AllConcs9414_LF.csv')

#Need dischareg data as well
q<-read.csv('AllHDMQ9414_hourly.csv')

#fix dates to POSIX format
lf.data$LF.Date_time<-ymd_hm(lf.data$LF.Date_time)
q$Date.Time<-ymd_hm(q$Date.Time)

lf.discharge<-data.frame(dateTime=q$Date.Time,Q.m3.sec=q$LFQ.m3.sec)
lf.discharge<-lf.discharge[!is.na(lf.discharge$Q.m3.sec),]
lf.discharge<-lf.discharge[lf.discharge$Q.m3.sec>0,]
#add year to discharge
lf.discharge$year<-format(ymd_hms(lf.discharge$dateTime),'%Y')
lf.discharge$dateTime<-ymd_hms(lf.discharge$dateTime)


lf.q<-q[,c(1,4)]
colnames(lf.q)[1]<-'dateTime'

#get rid of NA columns
lf.data<-lf.data[,1:5]

#make different data frames for TDN, TDP, Particulate N, and Particulate P
lf.tdn<-lf.data$LF.NH4.ug.L*1.013+508.5 #R2  = 0.821
lf.tdn<-data.frame(dateTime=lf.data$LF.Date_time,TDN.ug.L=lf.tdn)
lf.tdn<-merge(lf.tdn,lf.q,by='dateTime')
colnames(lf.tdn)<-c('dateTime','TDN.ug.L','Q.m3.sec')
lf.tdn<-lf.tdn[!is.na(lf.tdn$dateTime),]
lf.tdn<-lf.tdn[lf.tdn$Q.m3.sec>0,]

lf.tdp<-lf.data$LF.SRP.ug.L*1.190+8.122 #R2 = 0.816
lf.tdp<-data.frame(dateTime=lf.data$LF.Date_time,TDP.ug.L=lf.tdp)
lf.tdp<-merge(lf.tdp,lf.q,by='dateTime')
lf.tdp<-lf.tdp[!is.na(lf.tdp$dateTime),]
colnames(lf.tdp)<-c('dateTime','TDP.ug.L','Q.m3.sec')
lf.tdp<-lf.tdp[!is.na(lf.tdp$dateTime),]
lf.tdp<-lf.tdp[lf.tdp$Q.m3.sec>0,]

lf.pn<-lf.data$LF.SS.mg.L*3.780+131.25 #R2 = 0.924
lf.pn<-data.frame(dateTime=lf.data$LF.Date_time,PN.ug.L=lf.pn)
lf.pn<-merge(lf.pn,lf.q,by='dateTime')
lf.pn<-lf.pn[!is.na(lf.pn$dateTime),]
colnames(lf.pn)<-c('dateTime','PN.ug.L','Q.m3.sec')
lf.pn<-lf.pn[!is.na(lf.pn$dateTime),]
lf.pn<-lf.pn[lf.pn$Q.m3.sec>0,]

lf.pp<-lf.data$LF.SS.mg.L*0.620+39.277 #R2 = 0.673
lf.pp<-data.frame(dateTime=lf.data$LF.Date_time,PP.ug.L=lf.pp)
lf.pp<-merge(lf.pp,lf.q,by='dateTime')
lf.pp<-lf.pp[!is.na(lf.pp$dateTime),]
colnames(lf.pp)<-c('dateTime','PP.ug.L','Q.m3.sec')
lf.pp<-lf.pp[!is.na(lf.pp$dateTime),]
lf.pp<-lf.pp[lf.pp$Q.m3.sec>0,]

#Start with TDN for Four Mile Creek
#use the composite method
#need to first trim the dataset to only use every seven days
regr<-TRUE
for(i in 2:nrow(lf.tdn)){
	if(as.Date(lf.tdn$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(lf.tdn$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr[is.na(regr)]=FALSE
lf.tdn$REGR=regr
lf.tdp$REGR<-regr
lf.pn$REGR<-regr
lf.pp$REGR<-regr

#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
lf.tdn$REGR[lf.tdn$TDN.ug.L>=(mean(lf.tdn$TDN.ug.L,na.rm=T)+(2*sd(lf.tdn$TDN.ug.L,na.rm=T)))]=TRUE

#add year to data frame
lf.tdn$year<-format(ymd_hms(lf.tdn$dateTime),'%Y')

#remove NAs
lf.tdn<-lf.tdn[!is.na(lf.tdn$TDN.ug.L),]
lf.tdn<-lf.tdn[!is.na(lf.tdn$Q.m3.sec),]

#subset data for interpolation purposes
lf.reg.tdn<-subset(lf.tdn,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='TDN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Little Four Mile Creek')

years<-unique(lf.tdn$year)


lf.tdn.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(TDN.ug.L~model(1),data=lf.reg.tdn[lf.reg.tdn$year>=years[i] & lf.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(TDN.ug.L~model(2),data=lf.reg.tdn[lf.reg.tdn$year>=years[i] & lf.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(TDN.ug.L~model(3),data=lf.reg.tdn[lf.reg.tdn$year>=years[i] & lf.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(TDN.ug.L~model(4),data=lf.reg.tdn[lf.reg.tdn$year>=years[i] & lf.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(TDN.ug.L~model(5),data=lf.reg.tdn[lf.reg.tdn$year>=years[i] & lf.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(TDN.ug.L~model(6),data=lf.reg.tdn[lf.reg.tdn$year>=years[i] & lf.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(TDN.ug.L~model(7),data=lf.reg.tdn[lf.reg.tdn$year>=years[i] & lf.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(TDN.ug.L~model(8),data=lf.reg.tdn[lf.reg.tdn$year>=years[i] & lf.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(TDN.ug.L~model(9),data=lf.reg.tdn[lf.reg.tdn$year>=years[i] & lf.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


lf.tdn_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=lf.tdn[lf.tdn$year>=years[i] & lf.tdn$year<=years[i+2],])
lc.conc.pred<-predictSolute(lf.tdn_lc,'conc',lf.discharge[lf.discharge$year>=years[i] & lf.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(lf.tdn_lc,'flux',lf.discharge[lf.discharge$year>=years[i] & lf.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
lf.tdn.comp<-rbind(lf.tdn.comp,lc.pred)

print(i)
}

#Now do TDP____________________________________
#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
lf.tdp$REGR[lf.tdp$TDP.ug.L>=(mean(lf.tdp$TDP.ug.L,na.rm=T)+(2*sd(lf.tdp$TDP.ug.L,na.rm=T)))]=TRUE

#add year to data frame
lf.tdp$year<-format(ymd_hms(lf.tdp$dateTime),'%Y')

#remove NAs
lf.tdp<-lf.tdp[!is.na(lf.tdp$TDP.ug.L),]
lf.tdp<-lf.tdp[!is.na(lf.tdp$Q.m3.sec),]

#subset data for interpolation purposes
lf.reg.tdp<-subset(lf.tdp,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='TDP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Little Four Mile Creek')


lf.tdp.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(TDP.ug.L~model(1),data=lf.reg.tdp[lf.reg.tdp$year>=years[i] & lf.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(TDP.ug.L~model(2),data=lf.reg.tdp[lf.reg.tdp$year>=years[i] & lf.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(TDP.ug.L~model(3),data=lf.reg.tdp[lf.reg.tdp$year>=years[i] & lf.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(TDP.ug.L~model(4),data=lf.reg.tdp[lf.reg.tdp$year>=years[i] & lf.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(TDP.ug.L~model(5),data=lf.reg.tdp[lf.reg.tdp$year>=years[i] & lf.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(TDP.ug.L~model(6),data=lf.reg.tdp[lf.reg.tdp$year>=years[i] & lf.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(TDP.ug.L~model(7),data=lf.reg.tdp[lf.reg.tdp$year>=years[i] & lf.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(TDP.ug.L~model(8),data=lf.reg.tdp[lf.reg.tdp$year>=years[i] & lf.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(TDP.ug.L~model(9),data=lf.reg.tdp[lf.reg.tdp$year>=years[i] & lf.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


lf.tdp_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=lf.tdp[lf.tdp$year>=years[i] & lf.tdp$year<=years[i+2],])
lc.conc.pred<-predictSolute(lf.tdp_lc,'conc',lf.discharge[lf.discharge$year>=years[i] & lf.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(lf.tdp_lc,'flux',lf.discharge[lf.discharge$year>=years[i] & lf.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
lf.tdp.comp<-rbind(lf.tdp.comp,lc.pred)

print(i)
}


#Now do PN____________________________
#Now do PN____________________________________
#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
lf.pn$REGR[lf.pn$PN.ug.L>=(mean(lf.pn$PN.ug.L,na.rm=T)+(2*sd(lf.pn$PN.ug.L,na.rm=T)))]=TRUE

#add year to data frame
lf.pn$year<-format(ymd_hms(lf.pn$dateTime),'%Y')

#remove NAs
lf.pn<-lf.pn[!is.na(lf.pn$PN.ug.L),]
lf.pn<-lf.pn[!is.na(lf.pn$Q.m3.sec),]

#subset data for interpolation purposes
lf.reg.pn<-subset(lf.pn,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='PN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Little Four Mile Creek')


lf.pn.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(PN.ug.L~model(1),data=lf.reg.pn[lf.reg.pn$year>=years[i] & lf.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(PN.ug.L~model(2),data=lf.reg.pn[lf.reg.pn$year>=years[i] & lf.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(PN.ug.L~model(3),data=lf.reg.pn[lf.reg.pn$year>=years[i] & lf.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(PN.ug.L~model(4),data=lf.reg.pn[lf.reg.pn$year>=years[i] & lf.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(PN.ug.L~model(5),data=lf.reg.pn[lf.reg.pn$year>=years[i] & lf.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(PN.ug.L~model(6),data=lf.reg.pn[lf.reg.pn$year>=years[i] & lf.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(PN.ug.L~model(7),data=lf.reg.pn[lf.reg.pn$year>=years[i] & lf.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(PN.ug.L~model(8),data=lf.reg.pn[lf.reg.pn$year>=years[i] & lf.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(PN.ug.L~model(9),data=lf.reg.pn[lf.reg.pn$year>=years[i] & lf.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


lf.pn_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=lf.pn[lf.pn$year>=years[i] & lf.pn$year<=years[i+2],])
lc.conc.pred<-predictSolute(lf.pn_lc,'conc',lf.discharge[lf.discharge$year>=years[i] & lf.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(lf.pn_lc,'flux',lf.discharge[lf.discharge$year>=years[i] & lf.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
lf.pn.comp<-rbind(lf.pn.comp,lc.pred)

print(i)
}

#NOW do PP____________________________
#Now do PP____________________________________
#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
lf.pp$REGR[lf.pp$PP.ug.L>=(mean(lf.pp$PP.ug.L,na.rm=T)+(2*sd(lf.pp$PP.ug.L,na.rm=T)))]=TRUE

#add year to data frame
lf.pp$year<-format(ymd_hms(lf.pp$dateTime),'%Y')

#remove NAs
lf.pp<-lf.pp[!is.na(lf.pp$PP.ug.L),]
lf.pp<-lf.pp[!is.na(lf.pp$Q.m3.sec),]

#subset data for interpolation purposes
lf.reg.pp<-subset(lf.pp,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='PP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Little Four Mile Creek')


lf.pp.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(PP.ug.L~model(1),data=lf.reg.pp[lf.reg.pp$year>=years[i] & lf.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(PP.ug.L~model(2),data=lf.reg.pp[lf.reg.pp$year>=years[i] & lf.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(PP.ug.L~model(3),data=lf.reg.pp[lf.reg.pp$year>=years[i] & lf.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(PP.ug.L~model(4),data=lf.reg.pp[lf.reg.pp$year>=years[i] & lf.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(PP.ug.L~model(5),data=lf.reg.pp[lf.reg.pp$year>=years[i] & lf.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(PP.ug.L~model(6),data=lf.reg.pp[lf.reg.pp$year>=years[i] & lf.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(PP.ug.L~model(7),data=lf.reg.pp[lf.reg.pp$year>=years[i] & lf.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(PP.ug.L~model(8),data=lf.reg.pp[lf.reg.pp$year>=years[i] & lf.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(PP.ug.L~model(9),data=lf.reg.pp[lf.reg.pp$year>=years[i] & lf.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


lf.pp_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=lf.pp[lf.pp$year>=years[i] & lf.pp$year<=years[i+2],])
lc.conc.pred<-predictSolute(lf.pp_lc,'conc',lf.discharge[lf.discharge$year>=years[i] & lf.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(lf.pp_lc,'flux',lf.discharge[lf.discharge$year>=years[i] & lf.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
lf.pp.comp<-rbind(lf.pp.comp,lc.pred)

print(i)
}

#Calculate TDN, TDP, PP, and PN for Mike using NO3, SRP, and SS data - this requires converting NO3, SRP, and SS into those constituent concentrations using established linear regressions. These are stream specific.

#Marshall's Branch
#load concentration data data
setwd('~/Documents/Miami U/Acton load data/R')

#Concentration data
mb.data<-read.csv('AllConcs9414_MB.csv')

#Need dischareg data as well
q<-read.csv('AllHDMQ9414_hourly.csv')

#fix dates to POSIX format
mb.data$MB.Date_time<-ymd_hm(mb.data$MB.Date_time)
q$Date.Time<-ymd_hm(q$Date.Time)

mb.discharge<-data.frame(dateTime=q$Date.Time,Q.m3.sec=q$MBQ.m3.sec)
mb.discharge<-mb.discharge[!is.na(mb.discharge$Q.m3.sec),]
mb.discharge<-mb.discharge[mb.discharge$Q.m3.sec>0,]
#add year to discharge
mb.discharge$year<-format(ymd_hms(mb.discharge$dateTime),'%Y')
mb.discharge$dateTime<-ymd_hms(mb.discharge$dateTime)


mb.q<-q[,c(1,4)]
colnames(mb.q)[1]<-'dateTime'

#get rid of NA columns
mb.data<-mb.data[,1:5]

#make different data frames for TDN, TDP, Particulate N, and Particulate P
mb.tdn<-mb.data$MB.NH4.ug.L*1.040+542.7 #R2  = 0.821
mb.tdn<-data.frame(dateTime=mb.data$MB.Date_time,TDN.ug.L=mb.tdn)
mb.tdn<-merge(mb.tdn,mb.q,by='dateTime')
colnames(mb.tdn)<-c('dateTime','TDN.ug.L','Q.m3.sec')
mb.tdn<-mb.tdn[!is.na(mb.tdn$dateTime),]
mb.tdn<-mb.tdn[mb.tdn$Q.m3.sec>0,]

mb.tdp<-mb.data$MB.SRP.ug.L*1.068+16.089 #R2 = 0.816
mb.tdp<-data.frame(dateTime=mb.data$MB.Date_time,TDP.ug.L=mb.tdp)
mb.tdp<-merge(mb.tdp,mb.q,by='dateTime')
mb.tdp<-mb.tdp[!is.na(mb.tdp$dateTime),]
colnames(mb.tdp)<-c('dateTime','TDP.ug.L','Q.m3.sec')
mb.tdp<-mb.tdp[!is.na(mb.tdp$dateTime),]
mb.tdp<-mb.tdp[mb.tdp$Q.m3.sec>0,]

mb.pn<-mb.data$MB.SS.mg.L*3.108+142.11 #R2 = 0.924
mb.pn<-data.frame(dateTime=mb.data$MB.Date_time,PN.ug.L=mb.pn)
mb.pn<-merge(mb.pn,mb.q,by='dateTime')
mb.pn<-mb.pn[!is.na(mb.pn$dateTime),]
colnames(mb.pn)<-c('dateTime','PN.ug.L','Q.m3.sec')
mb.pn<-mb.pn[!is.na(mb.pn$dateTime),]
mb.pn<-mb.pn[mb.pn$Q.m3.sec>0,]

mb.pp<-mb.data$MB.SS.mg.L*0.682+35.681 #R2 = 0.673
mb.pp<-data.frame(dateTime=mb.data$MB.Date_time,PP.ug.L=mb.pp)
mb.pp<-merge(mb.pp,mb.q,by='dateTime')
mb.pp<-mb.pp[!is.na(mb.pp$dateTime),]
colnames(mb.pp)<-c('dateTime','PP.ug.L','Q.m3.sec')
mb.pp<-mb.pp[!is.na(mb.pp$dateTime),]
mb.pp<-mb.pp[mb.pp$Q.m3.sec>0,]

#Start with TDN for Four Mile Creek
#use the composite method
#need to first trim the dataset to only use every seven days
regr<-TRUE
for(i in 2:nrow(mb.tdn)){
	if(as.Date(mb.tdn$dateTime[i],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')-as.Date(mb.tdn$dateTime[regr==TRUE][length(regr[regr==TRUE])],'%Y-%m-%d %H:%M:%S',tz='America/Jamaica')>=7){
		regr[i]=TRUE
	}
}
regr[is.na(regr)]=FALSE
mb.tdn$REGR=regr
mb.tdp$REGR<-regr
mb.pn$REGR<-regr
mb.pp$REGR<-regr

#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
mb.tdn$REGR[mb.tdn$TDN.ug.L>=(mean(mb.tdn$TDN.ug.L,na.rm=T)+(2*sd(mb.tdn$TDN.ug.L,na.rm=T)))]=TRUE

#add year to data frame
mb.tdn$year<-format(ymd_hms(mb.tdn$dateTime),'%Y')

#remove NAs
mb.tdn<-mb.tdn[!is.na(mb.tdn$TDN.ug.L),]
mb.tdn<-mb.tdn[!is.na(mb.tdn$Q.m3.sec),]

#subset data for interpolation purposes
mb.reg.tdn<-subset(mb.tdn,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='TDN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Marshalls Branch')

years<-unique(mb.tdn$year)


mb.tdn.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(TDN.ug.L~model(1),data=mb.reg.tdn[mb.reg.tdn$year>=years[i] & mb.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(TDN.ug.L~model(2),data=mb.reg.tdn[mb.reg.tdn$year>=years[i] & mb.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(TDN.ug.L~model(3),data=mb.reg.tdn[mb.reg.tdn$year>=years[i] & mb.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(TDN.ug.L~model(4),data=mb.reg.tdn[mb.reg.tdn$year>=years[i] & mb.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(TDN.ug.L~model(5),data=mb.reg.tdn[mb.reg.tdn$year>=years[i] & mb.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(TDN.ug.L~model(6),data=mb.reg.tdn[mb.reg.tdn$year>=years[i] & mb.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(TDN.ug.L~model(7),data=mb.reg.tdn[mb.reg.tdn$year>=years[i] & mb.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(TDN.ug.L~model(8),data=mb.reg.tdn[mb.reg.tdn$year>=years[i] & mb.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(TDN.ug.L~model(9),data=mb.reg.tdn[mb.reg.tdn$year>=years[i] & mb.reg.tdn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


mb.tdn_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=mb.tdn[mb.tdn$year>=years[i] & mb.tdn$year<=years[i+2],])
lc.conc.pred<-predictSolute(mb.tdn_lc,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(mb.tdn_lc,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
mb.tdn.comp<-rbind(mb.tdn.comp,lc.pred)

print(i)
}

#Now do TDP____________________________________
#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
mb.tdp$REGR[mb.tdp$TDP.ug.L>=(mean(mb.tdp$TDP.ug.L,na.rm=T)+(2*sd(mb.tdp$TDP.ug.L,na.rm=T)))]=TRUE

#add year to data frame
mb.tdp$year<-format(ymd_hms(mb.tdp$dateTime),'%Y')

#remove NAs
mb.tdp<-mb.tdp[!is.na(mb.tdp$TDP.ug.L),]
mb.tdp<-mb.tdp[!is.na(mb.tdp$Q.m3.sec),]

#subset data for interpolation purposes
mb.reg.tdp<-subset(mb.tdp,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='TDP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Marshalls Branch')


mb.tdp.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(TDP.ug.L~model(1),data=mb.reg.tdp[mb.reg.tdp$year>=years[i] & mb.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(TDP.ug.L~model(2),data=mb.reg.tdp[mb.reg.tdp$year>=years[i] & mb.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(TDP.ug.L~model(3),data=mb.reg.tdp[mb.reg.tdp$year>=years[i] & mb.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(TDP.ug.L~model(4),data=mb.reg.tdp[mb.reg.tdp$year>=years[i] & mb.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(TDP.ug.L~model(5),data=mb.reg.tdp[mb.reg.tdp$year>=years[i] & mb.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(TDP.ug.L~model(6),data=mb.reg.tdp[mb.reg.tdp$year>=years[i] & mb.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(TDP.ug.L~model(7),data=mb.reg.tdp[mb.reg.tdp$year>=years[i] & mb.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(TDP.ug.L~model(8),data=mb.reg.tdp[mb.reg.tdp$year>=years[i] & mb.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(TDP.ug.L~model(9),data=mb.reg.tdp[mb.reg.tdp$year>=years[i] & mb.reg.tdp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


mb.tdp_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=mb.tdp[mb.tdp$year>=years[i] & mb.tdp$year<=years[i+2],])
lc.conc.pred<-predictSolute(mb.tdp_lc,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(mb.tdp_lc,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
mb.tdp.comp<-rbind(mb.tdp.comp,lc.pred)

print(i)
}


#Now do PN____________________________
#Now do PN____________________________________
#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
mb.pn$REGR[mb.pn$PN.ug.L>=(mean(mb.pn$PN.ug.L,na.rm=T)+(2*sd(mb.pn$PN.ug.L,na.rm=T)))]=TRUE

#add year to data frame
mb.pn$year<-format(ymd_hms(mb.pn$dateTime),'%Y')

#remove NAs
mb.pn<-mb.pn[!is.na(mb.pn$PN.ug.L),]
mb.pn<-mb.pn[!is.na(mb.pn$Q.m3.sec),]

#subset data for interpolation purposes
mb.reg.pn<-subset(mb.pn,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='PN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Marshalls Branch')


mb.pn.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(PN.ug.L~model(1),data=mb.reg.pn[mb.reg.pn$year>=years[i] & mb.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(PN.ug.L~model(2),data=mb.reg.pn[mb.reg.pn$year>=years[i] & mb.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(PN.ug.L~model(3),data=mb.reg.pn[mb.reg.pn$year>=years[i] & mb.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(PN.ug.L~model(4),data=mb.reg.pn[mb.reg.pn$year>=years[i] & mb.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(PN.ug.L~model(5),data=mb.reg.pn[mb.reg.pn$year>=years[i] & mb.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(PN.ug.L~model(6),data=mb.reg.pn[mb.reg.pn$year>=years[i] & mb.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(PN.ug.L~model(7),data=mb.reg.pn[mb.reg.pn$year>=years[i] & mb.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(PN.ug.L~model(8),data=mb.reg.pn[mb.reg.pn$year>=years[i] & mb.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(PN.ug.L~model(9),data=mb.reg.pn[mb.reg.pn$year>=years[i] & mb.reg.pn$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


mb.pn_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=mb.pn[mb.pn$year>=years[i] & mb.pn$year<=years[i+2],])
lc.conc.pred<-predictSolute(mb.pn_lc,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(mb.pn_lc,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
mb.pn.comp<-rbind(mb.pn.comp,lc.pred)

print(i)
}

#NOW do PP____________________________
#Now do PP____________________________________
#include all samples where the concentration is greater than 2 sd's to get the high peaks in the regression models
mb.pp$REGR[mb.pp$PP.ug.L>=(mean(mb.pp$PP.ug.L,na.rm=T)+(2*sd(mb.pp$PP.ug.L,na.rm=T)))]=TRUE

#add year to data frame
mb.pp$year<-format(ymd_hms(mb.pp$dateTime),'%Y')

#remove NAs
mb.pp<-mb.pp[!is.na(mb.pp$PP.ug.L),]
mb.pp<-mb.pp[!is.na(mb.pp$Q.m3.sec),]

#subset data for interpolation purposes
mb.reg.pp<-subset(mb.pp,REGR)

#create a metadata description of the dataset and desired output
meta<-metadata(constituent='PP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Marshalls Branch')


mb.pp.comp<-c()
for(i in seq(1,length(years),3)){
#fit four models: interpolation, linear, rloadest, and composite. 
mods<-list()
lr.1<-loadReg2(loadReg(PP.ug.L~model(1),data=mb.reg.pp[mb.reg.pp$year>=years[i] & mb.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[1]]<-lr.1
lr.2<-loadReg2(loadReg(PP.ug.L~model(2),data=mb.reg.pp[mb.reg.pp$year>=years[i] & mb.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[2]]<-lr.2
lr.3<-loadReg2(loadReg(PP.ug.L~model(3),data=mb.reg.pp[mb.reg.pp$year>=years[i] & mb.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[3]]<-lr.3
lr.4<-loadReg2(loadReg(PP.ug.L~model(4),data=mb.reg.pp[mb.reg.pp$year>=years[i] & mb.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[4]]<-lr.4
lr.5<-loadReg2(loadReg(PP.ug.L~model(5),data=mb.reg.pp[mb.reg.pp$year>=years[i] & mb.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[5]]<-lr.5
lr.6<-loadReg2(loadReg(PP.ug.L~model(6),data=mb.reg.pp[mb.reg.pp$year>=years[i] & mb.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[6]]<-lr.6
lr.7<-loadReg2(loadReg(PP.ug.L~model(7),data=mb.reg.pp[mb.reg.pp$year>=years[i] & mb.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[7]]<-lr.7
lr.8<-loadReg2(loadReg(PP.ug.L~model(8),data=mb.reg.pp[mb.reg.pp$year>=years[i] & mb.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[8]]<-lr.8
lr.9<-loadReg2(loadReg(PP.ug.L~model(9),data=mb.reg.pp[mb.reg.pp$year>=years[i] & mb.reg.pp$year<=years[i+2],],flow='Q.m3.sec',dates='dateTime',time.step='instantaneous',flow.units='cms',conc.units='ug/L',load.units='kg')) #rloadest model
mods[[9]]<-lr.9
aics<-c(as.numeric(getFittedModel(lr.1)$cfit[25]),as.numeric(getFittedModel(lr.2)$cfit[25]),as.numeric(getFittedModel(lr.3)$cfit[25]),as.numeric(getFittedModel(lr.4)$cfit[25]),as.numeric(getFittedModel(lr.5)$cfit[25]),as.numeric(getFittedModel(lr.6)$cfit[25]),as.numeric(getFittedModel(lr.7)$cfit[25]),as.numeric(getFittedModel(lr.8)$cfit[25]),as.numeric(getFittedModel(lr.9)$cfit[25]))


mb.pp_lc<-loadComp(reg.model=mods[[which(aics==min(aics))]],interp.format='conc',interp.data=mb.pp[mb.pp$year>=years[i] & mb.pp$year<=years[i+2],])
lc.conc.pred<-predictSolute(mb.pp_lc,'conc',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.flux.pred<-predictSolute(mb.pp_lc,'flux',mb.discharge[mb.discharge$year>=years[i] & mb.discharge$year<=years[i+2],],se.pred=T,date=T)
lc.pred<-data.frame(date=lc.conc.pred$date,conc.pred=lc.conc.pred$fit,flux.pred=lc.flux.pred$fit,flux.se=lc.flux.pred$se.pred)
mb.pp.comp<-rbind(mb.pp.comp,lc.pred)

print(i)
}

setwd('~/Desktop')
save.image('totalNutrients.RData')