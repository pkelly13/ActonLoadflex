#Aggregate total nutrients into daily loads - For Mike

rm(list=ls())
setwd('~')
source('.Rprofile')

#load workspace 
setwd('~/Desktop')
load('totalNutrients.RData')

#Four Mile Creek
#TDN
meta<-metadata(constituent='TDN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

fm.tdn.data<-fm.tdn.comp[,c(1,3,4)]
colnames(fm.tdn.data)<-c('date','fit','se.pred')
fm.tdn.daily<-aggregateSolute(fm.tdn.data,meta,'flux total','day')
fm.tdn.monthly<-aggregateSolute(fm.tdn.data,meta,'flux total','month')

#TDP
meta<-metadata(constituent='TDP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

fm.tdp.data<-fm.tdp.comp[,c(1,3,4)]
colnames(fm.tdp.data)<-c('date','fit','se.pred')
fm.tdp.daily<-aggregateSolute(fm.tdp.data,meta,'flux total','day')
fm.tdp.monthly<-aggregateSolute(fm.tdp.data,meta,'flux total','month')

#PN
meta<-metadata(constituent='PN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

fm.pn.data<-fm.pn.comp[,c(1,3,4)]
colnames(fm.pn.data)<-c('date','fit','se.pred')
fm.pn.daily<-aggregateSolute(fm.pn.data,meta,'flux total','day')
fm.pn.monthly<-aggregateSolute(fm.pn.data,meta,'flux total','month')

#PP
meta<-metadata(constituent='PP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Four Mile Creek')

fm.pp.data<-fm.pp.comp[,c(1,3,4)]
colnames(fm.pp.data)<-c('date','fit','se.pred')
fm.pp.daily<-aggregateSolute(fm.pp.data,meta,'flux total','day')
fm.pp.monthly<-aggregateSolute(fm.pp.data,meta,'flux total','month')

#Little Four Mile Creek
#TDN
meta<-metadata(constituent='TDN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Little Four Mile Creek')

lf.tdn.data<-lf.tdn.comp[,c(1,3,4)]
colnames(lf.tdn.data)<-c('date','fit','se.pred')
lf.tdn.daily<-aggregateSolute(lf.tdn.data,meta,'flux total','day')
lf.tdn.monthly<-aggregateSolute(lf.tdn.data,meta,'flux total','month')

#TDP
meta<-metadata(constituent='TDP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Little Four Mile Creek')

lf.tdp.data<-lf.tdp.comp[,c(1,3,4)]
colnames(lf.tdp.data)<-c('date','fit','se.pred')
lf.tdp.daily<-aggregateSolute(lf.tdp.data,meta,'flux total','day')
lf.tdp.monthly<-aggregateSolute(lf.tdp.data,meta,'flux total','month')

#PN
meta<-metadata(constituent='PN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Little Four Mile Creek')

lf.pn.data<-lf.pn.comp[,c(1,3,4)]
colnames(lf.pn.data)<-c('date','fit','se.pred')
lf.pn.daily<-aggregateSolute(lf.pn.data,meta,'flux total','day')
lf.pn.monthly<-aggregateSolute(lf.pn.data,meta,'flux total','month')

#PP
meta<-metadata(constituent='PP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Little Four Mile Creek')

lf.pp.data<-lf.pp.comp[,c(1,3,4)]
colnames(lf.pp.data)<-c('date','fit','se.pred')
lf.pp.daily<-aggregateSolute(lf.pp.data,meta,'flux total','day')
lf.pp.monthly<-aggregateSolute(lf.pp.data,meta,'flux total','month')

#Marshalls Branch
#TDN
meta<-metadata(constituent='TDN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Marshalls Branch')

mb.tdn.data<-mb.tdn.comp[,c(1,3,4)]
colnames(mb.tdn.data)<-c('date','fit','se.pred')
mb.tdn.daily<-aggregateSolute(mb.tdn.data,meta,'flux total','day')
mb.tdn.monthly<-aggregateSolute(mb.tdn.data,meta,'flux total','month')

#TDP
meta<-metadata(constituent='TDP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Marshalls Branch')

mb.tdp.data<-mb.tdp.comp[,c(1,3,4)]
colnames(mb.tdp.data)<-c('date','fit','se.pred')
mb.tdp.daily<-aggregateSolute(mb.tdp.data,meta,'flux total','day')
mb.tdp.monthly<-aggregateSolute(mb.tdp.data,meta,'flux total','month')

#PN
meta<-metadata(constituent='PN.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Marshalls Branch')

mb.pn.data<-mb.pn.comp[,c(1,3,4)]
colnames(mb.pn.data)<-c('date','fit','se.pred')
mb.pn.daily<-aggregateSolute(mb.pn.data,meta,'flux total','day')
mb.pn.monthly<-aggregateSolute(mb.pn.data,meta,'flux total','month')

#PP
meta<-metadata(constituent='PP.ug.L',flow='Q.m3.sec',dates='dateTime',conc.units='ug L^-1',flow.units='cms',load.units='kg',load.rate.units='kg d^-1',station='Marshalls Branch')

mb.pp.data<-mb.pp.comp[,c(1,3,4)]
colnames(mb.pp.data)<-c('date','fit','se.pred')
mb.pp.daily<-aggregateSolute(mb.pp.data,meta,'flux total','day')
mb.pp.monthly<-aggregateSolute(mb.pp.data,meta,'flux total','month')

#Combine by nutrient
tdn<-merge(fm.tdn.daily,lf.tdn.daily,by='Day')
tdn<-merge(tdn,mb.tdn.daily,by='Day')
tdn<-tdn[,-c(4,9,14)]
colnames(tdn)<-c('Day','FM.TDN.load.ug.L','SE','CI_lower','CI_upper','LF.TDN.load.ug.L','SE','CI_lower','CI_upper','MB.TDN.load.ug.L','SE','CI_lower','CI_upper')

tdp<-merge(fm.tdp.daily,lf.tdp.daily,by='Day')
tdp<-merge(tdp,mb.tdp.daily,by='Day')
tdp<-tdp[,-c(4,9,14)]
colnames(tdp)<-c('Day','FM.TDP.load.ug.L','SE','CI_lower','CI_upper','LF.TDP.load.ug.L','SE','CI_lower','CI_upper','MB.TDP.load.ug.L','SE','CI_lower','CI_upper')

pn<-merge(fm.pn.daily,lf.pn.daily,by='Day')
pn<-merge(pn,mb.pn.daily,by='Day')
pn<-pn[,-c(4,9,14)]
colnames(pn)<-c('Day','FM.PN.load.ug.L','SE','CI_lower','CI_upper','LF.PN.load.ug.L','SE','CI_lower','CI_upper','MB.PN.load.ug.L','SE','CI_lower','CI_upper')

pp<-merge(fm.pp.daily,lf.pp.daily,by='Day')
pp<-merge(pp,mb.pp.daily,by='Day')
pp<-pp[,-c(4,9,14)]
colnames(pp)<-c('Day','FM.PP.load.ug.L','SE','CI_lower','CI_upper','LF.PP.load.ug.L','SE','CI_lower','CI_upper','MB.PP.load.ug.L','SE','CI_lower','CI_upper')


#Combine monthly data by nutrient
month.tdn<-merge(fm.tdn.monthly,lf.tdn.monthly,by='Month')
month.tdn<-merge(month.tdn,mb.tdn.monthly,by='Month')
month.tdn<-month.tdn[,-c(4,9,14)]
colnames(month.tdn)<-c('Month','FM.TDN.load.ug.L','SE','CI_lower','CI_upper','LF.TDN.load.ug.L','SE','CI_lower','CI_upper','MB.TDN.load.ug.L','SE','CI_lower','CI_upper')

month.tdp<-merge(fm.tdp.monthly,lf.tdp.monthly,by='Month')
month.tdp<-merge(month.tdp,mb.tdp.monthly,by='Month')
month.tdp<-month.tdp[,-c(4,9,14)]
colnames(month.tdp)<-c('Month','FM.TDP.load.ug.L','SE','CI_lower','CI_upper','LF.TDP.load.ug.L','SE','CI_lower','CI_upper','MB.TDP.load.ug.L','SE','CI_lower','CI_upper')

month.pn<-merge(fm.pn.monthly,lf.pn.monthly,by='Month')
month.pn<-merge(month.pn,mb.pn.monthly,by='Month')
month.pn<-month.pn[,-c(4,9,14)]
colnames(month.pn)<-c('Month','FM.PN.load.ug.L','SE','CI_lower','CI_upper','LF.PN.load.ug.L','SE','CI_lower','CI_upper','MB.PN.load.ug.L','SE','CI_lower','CI_upper')

month.pp<-merge(fm.pp.monthly,lf.pp.monthly,by='Month')
month.pp<-merge(month.pp,mb.pp.monthly,by='Month')
month.pp<-month.pp[,-c(4,9,14)]
colnames(month.pp)<-c('Month','FM.PP.load.ug.L','SE','CI_lower','CI_upper','LF.PP.load.ug.L','SE','CI_lower','CI_upper','MB.PP.load.ug.L','SE','CI_lower','CI_upper')

#Write to excel files
setwd('~/Desktop')
write.xlsx(tdn,'TotalNutrients.xlsx',sheetName='TDN',row.names=F)
write.xlsx(tdp,'TotalNutrients.xlsx',sheetName='TDP',append=T,row.names=F)
write.xlsx(pn,'TotalNutrients.xlsx',sheetName='PN',append=T,row.names=F)
write.xlsx(pp,'TotalNutrients.xlsx',sheetName='PP',append=T,row.names=F)

#Write monthly data to csv files
setwd('~/Documents/Miami U/Acton data/monthly load estimates')
write.csv(month.tdn,'TDN.csv')
write.csv(month.tdp,'TDP.csv')
write.csv(month.pn,'PN.csv')
write.csv(month.pp,'PP.csv')

