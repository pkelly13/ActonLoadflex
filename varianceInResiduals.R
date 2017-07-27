#Get summary statistics from new "best model" runs of loads

rm(list=ls())
setwd('~')
source('.Rprofile')

library(loadflex)
library(rloadest)

#get data
setwd('~/Documents/Miami U/loadflex/best model_data')
#start with Four Mile Creek

load('TRIMMED.FM.loads.RData')
load('allData_bestModelRegression.RData')

#get residuals
#regression model
fm.nh4.resid.lr<-getResiduals(fm.nh4_lrAll,'flux','absolute')
var(fm.nh4.resid.lr$Resid)

fm.no3.resid.lr<-getResiduals(fm.no3_lrAll,'flux','absolute')
var(fm.no3.resid.lr$Resid)

fm.srp.resid.lr<-getResiduals(fm.srp_lrAll,'flux','absolute')
var(fm.srp.resid.lr$Resid)

fm.ss.resid.lr<-getResiduals(fm.ss_lrAll,'flux','absolute')
var(fm.ss.resid.lr$Resid)

lfm.nh4.resid.lr<-getResiduals(lfm.nh4_lrAll,'flux','absolute')
var(lfm.nh4.resid.lr$Resid)

lfm.no3.resid.lr<-getResiduals(lfm.no3_lrAll,'flux','absolute')
var(lfm.no3.resid.lr$Resid)

lfm.srp.resid.lr<-getResiduals(lfm.srp_lrAll,'flux','absolute')
var(lfm.srp.resid.lr$Resid)

lfm.ss.resid.lr<-getResiduals(lfm.ss_lrAll,'flux','absolute')
var(lfm.ss.resid.lr$Resid)

mb.nh4.resid.lr<-getResiduals(mb.nh4_lrAll,'flux','absolute')
var(mb.nh4.resid.lr$Resid)

mb.no3.resid.lr<-getResiduals(mb.no3_lrAll,'flux','absolute')
var(mb.no3.resid.lr$Resid)

mb.srp.resid.lr<-getResiduals(mb.srp_lrAll,'flux','absolute')
var(mb.srp.resid.lr$Resid)

mb.ss.resid.lr<-getResiduals(mb.ss_lrAll,'flux','absolute')
var(mb.ss.resid.lr$Resid)


#Now for composite models
fm.nh4.resid<-rbind(summaryStats.lc.resid[[1]],summaryStats.lc.resid[[2]],summaryStats.lc.resid[[3]],summaryStats.lc.resid[[4]],summaryStats.lc.resid[[5]],summaryStats.lc.resid[[6]],summaryStats.lc.resid[[7]])

var(fm.nh4.resid$Resid)

fm.no3.resid<-rbind(summaryStats.no3.lc.resid[[1]],summaryStats.no3.lc.resid[[2]],summaryStats.no3.lc.resid[[3]],summaryStats.no3.lc.resid[[4]],summaryStats.no3.lc.resid[[5]],summaryStats.no3.lc.resid[[6]],summaryStats.no3.lc.resid[[7]])

var(fm.no3.resid$Resid)

fm.srp.resid<-rbind(summaryStats.srp.lc.resid[[1]],summaryStats.srp.lc.resid[[2]],summaryStats.srp.lc.resid[[3]],summaryStats.srp.lc.resid[[4]],summaryStats.srp.lc.resid[[5]],summaryStats.srp.lc.resid[[6]],summaryStats.srp.lc.resid[[7]])

var(fm.srp.resid$Resid)

fm.ss.resid<-rbind(summaryStats.ss.lc.resid[[1]],summaryStats.ss.lc.resid[[2]],summaryStats.ss.lc.resid[[3]],summaryStats.ss.lc.resid[[4]],summaryStats.ss.lc.resid[[5]],summaryStats.ss.lc.resid[[6]],summaryStats.ss.lc.resid[[7]])

var(fm.ss.resid$Resid)

#Little Four Mile Creek
rm(list=ls())
setwd('~')
source('.Rprofile')

library(loadflex)
library(rloadest)

#get data
setwd('~/Documents/Miami U/loadflex/best model_data')
#start with Four Mile Creek

load('TRIMMED.MB.loads.RData')

lfm.nh4.resid<-rbind(summaryStats.lc.resid[[1]],summaryStats.lc.resid[[2]],summaryStats.lc.resid[[3]],summaryStats.lc.resid[[4]],summaryStats.lc.resid[[5]],summaryStats.lc.resid[[6]],summaryStats.lc.resid[[7]])

var(lfm.nh4.resid$Resid)

lfm.no3.resid<-rbind(summaryStats.no3.lc.resid[[1]],summaryStats.no3.lc.resid[[2]],summaryStats.no3.lc.resid[[3]],summaryStats.no3.lc.resid[[4]],summaryStats.no3.lc.resid[[5]],summaryStats.no3.lc.resid[[6]],summaryStats.no3.lc.resid[[7]])

var(lfm.no3.resid$Resid)

lfm.srp.resid<-rbind(summaryStats.srp.lc.resid[[1]],summaryStats.srp.lc.resid[[2]],summaryStats.srp.lc.resid[[3]],summaryStats.srp.lc.resid[[4]],summaryStats.srp.lc.resid[[5]],summaryStats.srp.lc.resid[[6]],summaryStats.srp.lc.resid[[7]])

var(lfm.srp.resid$Resid)

lfm.ss.resid<-rbind(summaryStats.ss.lc.resid[[1]],summaryStats.ss.lc.resid[[2]],summaryStats.ss.lc.resid[[3]],summaryStats.ss.lc.resid[[4]],summaryStats.ss.lc.resid[[5]],summaryStats.ss.lc.resid[[6]],summaryStats.ss.lc.resid[[7]])

var(lfm.ss.resid$Resid)