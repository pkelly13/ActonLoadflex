#This script looks at uncertainty meausrements for each stream and constituent
#Patrick Kelly 20 April 2016

rm(list=ls())
setwd('~')
source('.Rprofile')

library(loadflex)
library(rloadest)

#Start with Four Mile Creek
#load data
setwd('~/Documents/Miami U/loadflex/best model_data')
load('TRIMMED.LF.loads.RData')

#start by estimating daily loads 
#need to change the column headings and get rid of concentration data
meta.nh4<-metadata(constituent='NH4.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')


lc.nh4.data<-lc.data[,c(1,3,4)]
colnames(lc.nh4.data)<-c('date','fit','se.pred')
daily.nh4<-aggregateSolute(lc.nh4.data,meta.nh4,'flux total','day')
monthly.nh4<-aggregateSolute(lc.nh4.data,meta.nh4,'flux total','month')
yearly.nh4<-aggregateSolute(lc.nh4.data,meta.nh4,'flux total','calendar year')
wyearly.nh4<-aggregateSolute(lc.nh4.data,meta.nh4,'flux total','water year')

monthly.nh4.conc<-aggregateSolute(lc.nh4.data,meta.nh4,'conc','month')


lr.nh4.data<-lr.data[,c(1,4,5)]
colnames(lr.nh4.data)<-c('date','fit','se.pred')
monthly.lr.nh4<-aggregateSolute(lr.nh4.data,meta.nh4,'flux total','month')
daily.lr.nh4<-aggregateSolute(lr.nh4.data,meta.nh4,'flux total','day')
yearly.lr.nh4<-aggregateSolute(lr.nh4.data,meta.nh4,'flux total','calendar year')
wyearly.lr.nh4<-aggregateSolute(lr.nh4.data,meta.nh4,'flux total','water year')

li.nh4.data<-li.data[,c(1,4,5)]
colnames(li.nh4.data)<-c('date','fit','se.pred')
monthly.li.nh4<-aggregateSolute(li.nh4.data,meta.nh4,'flux total','month')
daily.li.nh4<-aggregateSolute(li.nh4.data,meta.nh4,'flux total','day')
yearly.li.nh4<-aggregateSolute(li.nh4.data,meta.nh4,'flux total','calendar year')
wyearly.li.nh4<-aggregateSolute(li.nh4.data,meta.nh4,'flux total','water year')

lr.nh4.reduced.data<-lr.data.reduced[,c(1,4,5)]
colnames(lr.nh4.reduced.data)<-c('date','fit','se.pred')
daily.lr.reduced.nh4<-aggregateSolute(lr.nh4.data,meta.nh4,'flux total','day')
monthly.lr.reduced.nh4<-aggregateSolute(lr.nh4.reduced.data,meta.nh4,'flux total','month')
yearly.lr.reduced.nh4<-aggregateSolute(lr.nh4.reduced.data,meta.nh4,'flux total','calendar year')
wyearly.lr.reduced.nh4<-aggregateSolute(lr.nh4.reduced.data,meta.nh4,'flux total','water year')


#NO3
lc.no3.data<-lc.no3.data[,c(1,3,4)]
lr.no3.data<-lr.no3.data[,c(1,4,5)]
colnames(lc.no3.data)<-c('date','fit','se.pred')
colnames(lr.no3.data)<-c('date','fit','se.pred')

daily.no3<-aggregateSolute(lc.no3.data,meta.no3,'flux total','day')
monthly.no3<-aggregateSolute(lc.no3.data,meta.no3,'flux total','month')
monthly.lr.no3<-aggregateSolute(lr.no3.data,meta.no3,'flux total','month')
yearly.no3<-aggregateSolute(lc.no3.data,meta.no3,'flux total','calendar year')
wyearly.no3<-aggregateSolute(lc.no3.data,meta.no3,'flux total','water year')

monthly.lr.no3<-aggregateSolute(lr.no3.data,meta.no3,'flux total','month')
daily.lr.no3<-aggregateSolute(lr.no3.data,meta.no3,'flux total','day')
yearly.lr.no3<-aggregateSolute(lr.no3.data,meta.no3,'flux total','calendar year')
wyearly.lr.no3<-aggregateSolute(lr.no3.data,meta.no3,'flux total','water year')

li.no3.data<-li.data[,c(1,4,5)]
colnames(li.no3.data)<-c('date','fit','se.pred')
monthly.li.no3<-aggregateSolute(li.no3.data,meta.no3,'flux total','month')
daily.li.no3<-aggregateSolute(li.no3.data,meta.no3,'flux total','day')
yearly.li.no3<-aggregateSolute(li.no3.data,meta.no3,'flux total','calendar year')
wyearly.li.no3<-aggregateSolute(li.no3.data,meta.no3,'flux total','water year')

monthly.li.no3.conc<-aggregateSolute(li.no3.data,meta.no3,'conc','month')

#srp
lc.srp.data<-lc.srp.data[,c(1,3,4)]
lr.srp.data<-lr.srp.data[,c(1,4,5)]
colnames(lc.srp.data)<-c('date','fit','se.pred')
colnames(lr.srp.data)<-c('date','fit','se.pred')

daily.srp<-aggregateSolute(lc.srp.data,meta.srp,'flux total','day')
monthly.srp<-aggregateSolute(lc.srp.data,meta.srp,'flux total','month')
monthly.lr.srp<-aggregateSolute(lr.srp.data,meta.srp,'flux total','month')
yearly.srp<-aggregateSolute(lc.srp.data,meta.srp,'flux total','calendar year')
wyearly.srp<-aggregateSolute(lc.srp.data,meta.srp,'flux total','water year')

monthly.srp.conc<-aggregateSolute(lc.srp.data,meta.srp,'conc','month')

monthly.lr.srp<-aggregateSolute(lr.srp.data,meta.srp,'flux total','month')
daily.lr.srp<-aggregateSolute(lr.srp.data,meta.srp,'flux total','day')
yearly.lr.srp<-aggregateSolute(lr.srp.data,meta.srp,'flux total','calendar year')
wyearly.lr.srp<-aggregateSolute(lr.srp.data,meta.srp,'flux total','water year')

li.srp.data<-li.data[,c(1,4,5)]
colnames(li.srp.data)<-c('date','fit','se.pred')
monthly.li.srp<-aggregateSolute(li.srp.data,meta.srp,'flux total','month')
daily.li.srp<-aggregateSolute(li.srp.data,meta.srp,'flux total','day')
yearly.li.srp<-aggregateSolute(li.srp.data,meta.srp,'flux total','calendar year')
wyearly.li.srp<-aggregateSolute(li.srp.data,meta.srp,'flux total','water year')

#ss
lc.ss.data<-lc.ss.data[,c(1,3,4)]
lr.ss.data<-lr.ss.data[,c(1,4,5)]
li.ss.data<-li.ss.data[,c(1,4,5)]
colnames(lc.ss.data)<-c('date','fit','se.pred')
colnames(lr.ss.data)<-c('date','fit','se.pred')
colnames(li.ss.data)<-c('date','fit','se.pred')

daily.ss<-aggregateSolute(lc.ss.data,meta.ss,'flux total','day')
monthly.ss<-aggregateSolute(lc.ss.data,meta.ss,'flux total','month')
monthly.lr.ss<-aggregateSolute(lr.ss.data,meta.ss,'flux total','month')
yearly.ss<-aggregateSolute(lc.ss.data,meta.ss,'flux total','calendar year')
wyearly.ss<-aggregateSolute(lc.ss.data,meta.ss,'flux total','water year')
monthly.li.ss<-aggregateSolute(li.ss.data,meta.ss,'flux total','month')

monthly.li.ss.conc<-aggregateSolute(li.ss.data,meta.ss,'conc','month')

monthly.lr.ss<-aggregateSolute(lr.ss.data,meta.ss,'flux total','month')
daily.lr.ss<-aggregateSolute(lr.ss.data,meta.ss,'flux total','day')
yearly.lr.ss<-aggregateSolute(lr.ss.data,meta.ss,'flux total','calendar year')
wyearly.lr.ss<-aggregateSolute(lr.ss.data,meta.ss,'flux total','water year')

monthly.li.ss<-aggregateSolute(li.ss.data,meta.ss,'flux total','month')
daily.li.ss<-aggregateSolute(li.ss.data,meta.ss,'flux total','day')
yearly.li.ss<-aggregateSolute(li.ss.data,meta.ss,'flux total','calendar year')
wyearly.li.ss<-aggregateSolute(li.ss.data,meta.ss,'flux total','water year')

#save annual and monthly estimates for use in MAR and other data analyses
setwd('~/Documents/Miami U/acton data/monthly load estimates')
write.csv(monthly.nh4,'LF.monthlyNH4Load_comp.csv')
write.csv(monthly.no3,'LF.monthlyNO3Load_comp.csv')
write.csv(monthly.li.no3,'LF.monthlyNO3Load_interp.csv')
write.csv(monthly.srp,'LF.monthlySRPLoad_comp.csv')
write.csv(monthly.li.ss,'LF.monthlySSLoad_interp.csv')
write.csv(monthly.nh4.conc,'LF.monthlyNH4conc_comp.csv')
write.csv(monthly.li.no3.conc,'LF.monthlyNO3conc_interp.csv')
write.csv(monthly.srp.conc,'LF.monthlySRPconc_comp.csv')
write.csv(monthly.li.ss.conc,'LF.monthlySSconc_interp.csv')
write.csv(monthly.ss,'LF.monthlySSLoad_comp.csv')


setwd('~/Documents/Miami U/acton data/annual load estimates')
write.csv(yearly.nh4,'LF.annualNH4Load_comp.csv')
write.csv(yearly.no3,'LF.annualNO3Load_comp.csv')
write.csv(yearly.li.no3,'LF.annualNO3Load_interp.csv')
write.csv(yearly.srp,'LF.annualSRPLoad_comp.csv')
write.csv(yearly.li.ss,'LF.annualSSLoad_interp.csv')

write.csv(wyearly.nh4,'LF.annualNH4Load_comp.csv')
write.csv(wyearly.no3,'LF.annualNO3Load_comp.csv')
write.csv(wyearly.li.no3,'LF.annualNO3Load_interp.csv')
write.csv(wyearly.srp,'LF.annualSRPLoad_comp.csv')
write.csv(wyearly.li.ss,'LF.annualSSLoad_interp.csv')


setwd('~/Documents/Miami U/loadflex/aggregated data')
save.image('LF_aggregatedFlux.RData')