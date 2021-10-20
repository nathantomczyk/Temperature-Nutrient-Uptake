#load libraries

library(ggplot2)
library(tidyverse)
library(lubridate)
library(reshape2)
library("cowplot")


######## Read in data
setwd("C:\\Users\\nt78066\\OneDrive - University of Georgia\\Documents\\Warming_uptake_stuff\\Budget_uptake_whole_stream\\")

d<-read.csv("whole_stream_nutrients_master.csv")
q<-read.csv("C:\\Users\\nt78066\\OneDrive - University of Georgia\\Documents\\Warming_uptake_stuff\\whole_stream_uptake\\crews_paireddischarge_all_2july2021_long.csv")
q.diff<-read.csv("discharge_top_and_bottom.csv")
widths<-read.csv("paired.wettedwidth.Master.csv")

# format dates and times
d$sample.date<-ymd(as.character(d$sample.date))
q$datetime<-ymd_hms(q$datetime)
widths$sample.date<-ymd(widths$sample.date)

#### cleaning and aggregating nutrient data
d$location.nutrient<-as.factor(paste(d$stream,d$stream.location,d$analyte,sep="-"))


id<-levels(d$location.nutrient)

d2<-data.frame()

for (i in 1:length(id)){
  s<-d[d$location.nutrient==id[i],]
  s$ignore.value<-"FALSE"
  mean.conc<-mean(s$original.concentration.ugL,na.rm=TRUE)
  sd.conc<-sd(s$original.concentration.ugL,na.rm=TRUE)
  
  s<-s[is.na(s$original.concentration.ugL)==FALSE,]
  
  s[s$original.concentration.ugL>(mean.conc+2*sd.conc),"ignore.value"]<-TRUE
  s[s$original.concentration.ugL<(mean.conc-2*sd.conc),"ignore.value"]<-TRUE
  
  
  d2<-rbind(d2,s)
  
}

d.means<-aggregate(original.concentration.ugL~analyte+stream+stream.location+sample.date,data=d2[d2$ignore.value==FALSE,],mean)
########## aggregating widths

widths2<-aggregate(wet.width.cm~sample.date+stream,data=widths,mean)

##########

## Nutreints were typically sampled from the reference stream around 9:00
q_nine<-q[hour(q$datetime)==9,]
q_nine$date.only<-date(q_nine$datetime)

q_nine_ref<-aggregate(discharge.ls~date.only,data=q_nine[q_nine$stream=="ws55",],mean)



dates<-levels(as.factor(d$sample.date))

##### calculating the change in nutrients in the ref stream.

concentration.changes<-data.frame(stream=rep(NA,200),din.change=rep(NA,200),srp.change=rep(NA,200),date=rep(NA,200),discharge=rep(NA,200))

for (i in 1:length(dates)){
  srp.change<-NA
  din.change<-NA
  
  try(ns<-d.means[d.means$sample.date==dates[i] & d.means$stream=="ws55",]) # subset nutreint data
  
  discharge<-q_nine_ref[q_nine_ref$date==dates[i],]
  
  srp.change<-ns[ns$analyte=="SRP" & ns$stream.location=="flume","original.concentration.ugL"]-ns[ns$analyte=="SRP" & ns$stream.location=="seep","original.concentration.ugL"]
  din.change<-ns[ns$analyte=="NH4" & ns$stream.location=="flume","original.concentration.ugL"]+
    ns[ns$analyte=="NO3" & ns$stream.location=="flume","original.concentration.ugL"]-
    ns[ns$analyte=="NO3" & ns$stream.location=="seep","original.concentration.ugL"]-
    ns[ns$analyte=="NH4" & ns$stream.location=="seep","original.concentration.ugL"]
  
  concentration.changes$stream[i]<-"ws55"
  try(concentration.changes$din.change[i]<-din.change)
  try(concentration.changes$srp.change[i]<-srp.change)
  concentration.changes$date[i]<-dates[i]
  try(concentration.changes$discharge[i]<-discharge$discharge.ls)

  
}


ggplot(concentration.changes,aes(x=ymd(date),y=din.change))+geom_point()+geom_smooth(se=FALSE)+geom_hline(yintercept=0,lty=2)+
  theme_classic()+geom_line()

ggplot(concentration.changes,aes(x=ymd(date),y=srp.change))+geom_point()+geom_smooth(se=FALSE)+geom_hline(yintercept=0,lty=2)+
  theme_classic()+geom_line()


top.bottom.discharge<-lm(log(meas.q.ls.top.of.reach)~log(meas.q.ls),data=q.diff[q.diff$date!="2018-04-25",]) ## modeling paired measurements of discharge at the top and bottom of the reach
tbd.summary<-summary(top.bottom.discharge)

fit.data<-data.frame(flume.q=seq(0,11,by=0.01))
fit.data$top.of.reach<-exp(coefficients(top.bottom.discharge)[1])*fit.data$flume.q**coefficients(top.bottom.discharge)[2]


ggplot(q.diff,aes(x=meas.q.ls,y=meas.q.ls.top.of.reach))+geom_point()+theme_classic()+theme(text = element_text(size=20))+
  xlim(0,11)+ylim(0,11)+geom_line(data=fit.data,aes(x=flume.q,y=top.of.reach),size=2)+xlab("Downstream discharge L/s")+
  ylab("Upstream discharge L/s")
  
  
                                                                                              
                                                                                              
## Times only

q_noon<-q[hour(q$datetime)==12,]
q_noon$date.only<-date(q_noon$datetime)

q_noon_warm<-aggregate(discharge.ls~date.only,data=q_noon[q_noon$stream=="towr",],mean)
q_noon_warm$upstream.discharge.ls<-exp(coefficients(top.bottom.discharge)[1])*q_noon_warm$discharge.ls**coefficients(top.bottom.discharge)[2]

###### Estimating the difference in concentraiton at the warm stream

for (i in 1:length(dates)){
  srp.change<-NA
  din.change<-NA
  
  try(ns<-d.means[d.means$sample.date==dates[i] & d.means$stream=="towr",]) # subset nutreint data
  
  discharge<-q_noon_warm[q_noon_warm$date.only==dates[i],]
  
  portion.upstream.flow<-discharge$upstream.discharge.ls/discharge$discharge.ls
  try(if (portion.upstream.flow>1){portion.upstream=1})
  
  added.srp<-ns[ns$analyte=="SRP" & ns$stream.location=="upstream","original.concentration.ugL"]*portion.upstream.flow+
    ns[ns$analyte=="SRP" & ns$stream.location=="seep","original.concentration.ugL"]*(1-portion.upstream.flow)
  
  added.din<-(ns[ns$analyte=="NH4" & ns$stream.location=="upstream","original.concentration.ugL"]+
                ns[ns$analyte=="NO3" & ns$stream.location=="upstream","original.concentration.ugL"])*portion.upstream.flow+
    (ns[ns$analyte=="NH4" & ns$stream.location=="seep","original.concentration.ugL"]+
       ns[ns$analyte=="NO3" & ns$stream.location=="seep","original.concentration.ugL"])*(1-portion.upstream.flow)
  
  srp.change<-ns[ns$analyte=="SRP" & ns$stream.location=="flume","original.concentration.ugL"]-added.srp
  din.change<-ns[ns$analyte=="NH4" & ns$stream.location=="flume","original.concentration.ugL"]+
    ns[ns$analyte=="NO3" & ns$stream.location=="flume","original.concentration.ugL"]-
    added.din
  
  concentration.changes$stream[i+length(dates)]<-"towr"
  try(concentration.changes$din.change[i+length(dates)]<-din.change)
  try(concentration.changes$srp.change[i+length(dates)]<-srp.change)
  concentration.changes$date[i+length(dates)]<-dates[i]
  try(concentration.changes$discharge[i+length(dates)]<-discharge$discharge.ls)
  
  
}


ggplot(concentration.changes,aes(x=ymd(date),y=din.change,color=stream))+geom_point()+geom_smooth(se=FALSE)+geom_hline(yintercept=0,lty=2)+
  theme_classic()+geom_line()

ggplot(concentration.changes,aes(x=ymd(date),y=srp.change,color=stream))+geom_point()+geom_smooth(se=FALSE)+geom_hline(yintercept=0,lty=2)+
  theme_classic()+geom_line()

######################################################### widths
widths2$date.only<-widths2$sample.date
ref.widths<-merge(widths2[widths2$stream=="WS55",],q_nine_ref,by=c("date.only"))

ggplot(ref.widths,aes(x=discharge.ls,y=wet.width.cm))+geom_point()+theme_classic()

ref.model<-lm(wet.width.cm~discharge.ls,data=ref.widths)
summary(ref.model)

warm.widths<-merge(widths2[widths2$stream=="TOWR",],q_noon_warm,by=c("date.only"))

warm.model<-lm(wet.width.cm~discharge.ls,data=warm.widths)
summary(warm.model)
ggplot(warm.widths,aes(x=discharge.ls,y=wet.width.cm))+geom_point()+theme_classic()

concentration.changes$width<-NA
concentration.changes[is.na(concentration.changes$discharge)!=TRUE & concentration.changes$stream=="ws55","width"]<-
  concentration.changes[is.na(concentration.changes$discharge)!=TRUE & concentration.changes$stream=="ws55","discharge"]*coefficients(ref.model)[2]+coefficients(ref.model)[1]

concentration.changes[is.na(concentration.changes$discharge)!=TRUE & concentration.changes$stream=="towr","width"]<-
  concentration.changes[is.na(concentration.changes$discharge)!=TRUE & concentration.changes$stream=="towr","discharge"]*coefficients(warm.model)[2]+coefficients(warm.model)[1]


######################################################### Calculating uptake rates

concentration.changes$uptake.rate.srp<-concentration.changes$srp.change*concentration.changes$discharge/((concentration.changes$width/100)*135)
concentration.changes$uptake.rate.din<-concentration.changes$din.change*concentration.changes$discharge/((concentration.changes$width/100)*135)


din.uptake<-ggplot(concentration.changes[complete.cases(concentration.changes),],aes(x=ymd(date),y=uptake.rate.din,color=stream))+
  geom_point()+
  geom_hline(yintercept=0,lty=2)+
  theme_classic()+geom_line()+
  ylab(expression("Net DIN uptake ("*mu*M~m^-2~s^-1*")"))+
  xlab("Date")+
  theme(text =element_text(size=20))+theme(legend.position = c(0.1, 0.91))+
  scale_color_manual(values=c("firebrick2","steelblue3"),name="")+
  geom_vline(xintercept=ymd("2019-04-01"),lty=2)
  

srp.uptake<-ggplot(concentration.changes[complete.cases(concentration.changes),],aes(x=ymd(date),y=uptake.rate.srp,color=stream))+geom_point()+
  geom_point()+
  geom_hline(yintercept=0,lty=2)+
  theme_classic()+geom_line()+
  ylab(expression("Net SRP uptake ("*mu*M~m^-2~s^-1*")"))+
  xlab("Date")+
  theme(text =element_text(size=20))+theme(legend.position = c(0.1, 0.91))+
  scale_color_manual(values=c("firebrick2","steelblue3"),name="")+
  geom_vline(xintercept=ymd("2019-04-01"),lty=2)


#########################################################

warm<-concentration.changes[concentration.changes$stream=="towr",]
ref<-concentration.changes[concentration.changes$stream=="ws55",]

test<-merge(warm,ref,by="date")

test$difference.in.DIN.change<-test$uptake.rate.din.x-test$uptake.rate.din.y
test$difference.in.SRP.change<-test$uptake.rate.srp.x-test$uptake.rate.srp.y


test<-test[complete.cases(test),]
test$treatment<-"reference"

test[ymd(test$date)>ymd("2019-05-06"),"treatment"]<-"warming"


summary(lm(difference.in.SRP.change~treatment,data=test))

summary(lm(difference.in.DIN.change~treatment,data=test))

library(Rmisc)

srp.means<-summarySE(test,measurevar = "difference.in.SRP.change",groupvars = "treatment")
din.means<-summarySE(test,measurevar = "difference.in.DIN.change",groupvars = "treatment")


srp.plot<-ggplot(srp.means,aes(x=treatment,y=-difference.in.SRP.change,fill=treatment))+geom_bar(stat="identity")+
  theme_classic()+
  theme(text =element_text(size=20))+theme(legend.position = c(0.28, 0.91))+
  scale_fill_manual(values=c("firebrick2","steelblue3"),name="")+
  ylab("Difference between streams")+xlab("")+
  geom_errorbar(aes(ymin=-difference.in.SRP.change-se,ymax=-difference.in.SRP.change+se),width=0.00001)
srp.plot

din.plot<-ggplot(din.means,aes(x=treatment,y=-difference.in.DIN.change,fill=treatment))+geom_bar(stat="identity")+
  theme_classic()+
  theme(text =element_text(size=20))+theme(legend.position = c(0.28, 0.91))+
  scale_fill_manual(values=c("firebrick2","steelblue3"),name="")+
  ylab("Difference between streams")+xlab("")+
  geom_errorbar(aes(ymin=-difference.in.DIN.change-se,ymax=-difference.in.DIN.change+se),width=0.00001)
din.plot

tiff(filename="Whole_stream_net_uptake_29sept2021.tiff",units="in",res=800,width=12,height=8,compression="lzw")
plot_grid(srp.uptake,srp.plot,din.uptake,din.plot,rel_widths = c(2, 1),labels="AUTO",cols=2,label_x=0.9,label_y=0.93,label_size = 20)
dev.off()