############ Notes
# I have the code worked up to estimate uptake in the streams and the R2 look pretty decent for all estimates
# Need to check out the discharge estiamtes and figure out what is up with the holes in the data
# Also want to check with phillip about how is was/ wasn't doing stage correction of the discharge
# Once these issues are worked out I *THINK* these estimates may be good to go. 
############



library(dplyr)
library(tidyverse)
library(lubridate)
library(cowplot)

q<-read.csv("./Whole_stream_gross_uptake/crews_paireddischarge_all_23Aug2021_long.csv")

release.info<-read.csv("./Whole_stream_gross_uptake/whole.stream.uptake.times.temps.release.info.csv")

chem<-read.csv("./Whole_stream_gross_uptake/uptake_temperature_nutrients.csv")

chem[chem$sample.id==10 & chem$date=="2021-01-26","po4.p.ug.L"]<-NA

widths<-read.csv("./Whole_stream_gross_uptake/stream_widths.csv")

cond<-read.csv("whole_stream_conductivity.csv")

chem$date<-ymd(chem$date)
widths$date<-ymd(widths$date)
cond$date<-ymd(cond$date)
q$datetime<-ymd_hms(q$datetime)

q$stream<-as.factor(q$stream)
q$stream<-recode_factor(q$stream,'ws55'="ref",'towr'="warm")

cond$cor.cond<-(cond$plt.spc-cond$background.spc)

cond2<-data.frame(date=cond$date,stream=cond$stream,corrected.cond=cond$cor.cond,meter=cond$meter)

rm(min)
rm(mean)

chem2<-aggregate(po4.p.ug.L~date+stream+meter+bg.plt,chem,mean)

plt.only<-chem2[chem2$bg.plt=="plt",]

bg.only<-chem[chem$bg.plt=="bkd",]



names(bg.only)[names(bg.only) == "po4.p.ug.L"] <- "background.PO4"
bg.only$id<-paste(bg.only$date,bg.only$stream,sep="-")
ggplot(bg.only[bg.only$stream=="ref",],aes(x=meter,y=background.PO4,color=id))+geom_point()+geom_smooth(se=FALSE)






chem.bg<-merge(plt.only,bg.only,by=c("stream","meter","date"))

chem3<-merge(chem.bg,cond2,by=c("date","meter","stream"))


##### percent dilution estimates #######

chem3$stream.date<-as.factor(paste(chem3$stream,chem3$date,sep="-"))

id<-levels(chem3$stream.date)

chem4<-data.frame()


for (i in 1:length(id)){
  s<-chem3[chem3$stream.date==id[i],]
  
  eighty<-s[s$meter==80,]
  
  s$dilution<-s$corrected.cond/eighty$corrected.cond
  

  
  chem4<-rbind(chem4,s)
  
}

############## Calculating uptake lengths

chem4$dist.from.eighty<-100-chem4$meter

#date.stream<-unique(data.frame(date=chem4$date,stream=chem4$stream))

uptake.lengths<-data.frame()

#chem4[chem4$date==ymd("2021-01-26") & chem4$stream=="ref" & chem4$dist.from.eighty==60,"bg.corrected.srp"]<-NA
#chem4[chem4$date==ymd("2021-01-26") & chem4$stream=="ref" & chem4$dist.from.eighty==60,"normalized.srp"]<-NA


for (i in 1:length(id)){
  s<-chem4[chem4$stream.date==id[i],]
  
  #bg<-chem4[chem4$stream==s$stream[1],]
  
  t<-test[test$stream==s$stream[1] & test$date==s$date[1],]
  
  s$bg.corrected.srp<-s$po4.p.ug.L -mean(s$background.PO4,na.rm=TRUE)
  
  s$dilution.corrected.srp<-s$bg.corrected.srp/s$dilution
  
  s$srp.remain.proportion<-s$dilution.corrected.srp/t$upstream.conc[1]
  

  #model<-lm(log(dilution.corrected.srp)~dist.from.eighty,s)
  
  model<-lm(log(srp.remain.proportion)~dist.from.eighty+0,s)
  model.sum<-summary(model)
  
  
  uptake.plot<- ggplot(s,aes(x=dist.from.eighty,y=log(srp.remain.proportion)))+geom_point()+theme_classic()+
    geom_abline(intrecept=0,slope=model.sum$coefficients[1])+xlim(0,100)+#ylim(0,log(.0001))+
    xlab("distance")+ylab("Normalized srp concentraiton")+
    annotate(x=80,y=log(0.8),"text",label=paste("R2=",round(model.sum$r.squared,2),sep=" "),size=6)+
    annotate(x=80,y=log(0.6),"text",label=paste(id[i]),size=6)
    #scale_y_log10()
  
  uptake.plot
  
  
  tiff(filename=paste(id[i],".tiff",sep=""),units="in",res=200,width=8,height=8,compression="lzw")
  print(uptake.plot)
  dev.off()
    
  output<-data.frame(stream=s$stream[1],date=s$date[1],k=model.sum$coefficients[1],sw=(-1/model.sum$coefficients[1]),r2=model.sum$r.squared)
  uptake.lengths<-rbind(uptake.lengths,output)
}



  
ggplot(uptake.lengths,aes(x=date,y=sw,color=stream))+geom_point()+geom_line()+theme_classic()  


############### Need to change values for the first three dates because the logger was malfunctioning, going
############### to just redo them all. 

model.q.data<- q %>%filter(datetime>=ymd_hms("2020-01-01 10:12:30") & datetime<=ymd_hms("2020-06-01 10:12:30"))

ggplot(model.q.data,aes(x=datetime,y=discharge.ls,color=stream))+geom_point()


model.q.data$datetime.h<-ymd_h(substring(as.character(model.q.data$datetime),1,13))

model.q.data.hourly<-aggregate(discharge.ls~stream+datetime.h,data=model.q.data,mean)

warm<-model.q.data.hourly[model.q.data.hourly$stream=="warm",]
ref<-model.q.data.hourly[model.q.data.hourly$stream=="ref",]

discharge.model.data<-merge(warm,ref,by="datetime.h")


ggplot(discharge.model.data,aes(y=discharge.ls.x,x=discharge.ls.y))+geom_point()+geom_smooth(method="lm",se=FALSE)+
  theme_classic()+xlab("Reference stream discharge L/s")+ylab("Warm stream discharge L/s")+
  theme(text = element_text(size=20))

discharge.model<-lm(discharge.ls.x~discharge.ls.y,data=discharge.model.data)
summary(discharge.model)

############### projecting discharge in the warm stream

needed.dates<- q %>% filter(datetime>=times$start.time2[1] & datetime<=times$end.time2[18] & stream=="ref")

warm.discharge<-data.frame(X=NA,datetime=needed.dates$datetime,year=2021,discharge.ls=coef(discharge.model)[1]+coef(discharge.model)[2]*needed.dates$discharge.ls,water_temp.c=NA,stream="warm")

q.fixed<-rbind(needed.dates,warm.discharge)

### pulling out right estimates of q

uptake.lengths$q<-NA
uptake.lengths$q.sd<-NA
uptake.lengths$width<-NA
uptake.lengths$temperature<-NA

times<-read.csv("whole.stream.uptake.times.temps.csv")

times$start.time2<-ymd_h(paste(times$date,substring(as.character(times$start.time),1,2),sep=" "))

times$end.time2<-ymd_h(paste(times$date,substring(as.character(times$end.time),1,2),sep=" "))

times$date2<-ymd(times$date)



for (i in 1:nrow(uptake.lengths)){
  times2<-times[times$stream==uptake.lengths$stream[i] & times$date2==uptake.lengths$date[i],]
  
  q2<-q.fixed %>% filter(datetime>=times2$start.time2 & datetime<=times2$end.time2 & stream==times2$stream)

  
  
    
  

  
  
  mean.q<-mean(q2$discharge.ls,na.rm=TRUE)
  sd.q<-sd(q2$discharge.ls,na.rm=TRUE)
  
  width2<-widths[widths$stream==uptake.lengths$stream[i] & widths$date==uptake.lengths$date[i],]
  mean.width<-mean(width2$width.cm,na.rm=TRUE)
  
  
  uptake.lengths$q[i]<-mean.q
  uptake.lengths$width[i]<-mean.width
  uptake.lengths$temperature[i]<-times2$temperature

}




#####

background.conc<-aggregate(po4.p.ug.L~stream+date,chem[chem$bg.plt=="bkd",],mean)


uptake.lengths<-merge(uptake.lengths,background.conc,by=c("stream","date"))

uptake.lengths$uptake.rate<-(uptake.lengths$q*1)/(uptake.lengths$width/100*uptake.lengths$sw)*60 ### estiamtes uptake rate in ug m2 min at a conc of 1 ug/L

ggplot(uptake.lengths,aes(x=date,y=uptake.rate,color=stream))+geom_point()+geom_line()+theme_classic()  

ggplot(uptake.lengths,aes(x=temperature,y=uptake.rate,color=stream))+geom_point()+geom_smooth(method="lm",se=FALSE)+theme_classic()  

Boltz<-8.6173303*10**-5 # boltzmann constant

uptake.lengths$boltz.temp<-1/((273+uptake.lengths$temperature)*Boltz)

uptake.lengths$log.rate<-log(uptake.lengths$uptake.rate)




model<-lm(log.rate~boltz.temp+stream,data=uptake.lengths[uptake.lengths$date!="2021-01-28",])
anova(model)
summary(model)

model2<-lm(log.rate~boltz.temp+stream,data=uptake.lengths)
anova(model2)
summary(model2)

model3<-lm(log.rate~boltz.temp*stream,data=uptake.lengths)
anova(model3)
summary(model3)


modelwarm<-lm(log.rate~boltz.temp,data=uptake.lengths[uptake.lengths$stream=="warm",])
summary(modelwarm)

modelref<-lm(log.rate~boltz.temp,data=uptake.lengths[uptake.lengths$stream=="ref",])
summary(modelref)


whole_stream_botlz<-ggplot(uptake.lengths,aes(x=boltz.temp,y=log.rate,color=stream))+theme_classic()+geom_point()+geom_smooth(method="lm",se=FALSE)+
  xlab("Inverse Temperature 1/kT")+
  ylab(expression("ln(SRP uptake rate ("*mu*g~"SRP hr"^-1~"m"^-2*"))"))+
  theme(text = element_text(size=20),legend.position = c(0.8, 0.2))+
  scale_color_manual(values=c("steelblue3","firebrick2"),name="")
whole_stream_botlz

qplot<-q %>% filter(datetime>=times$start.time2[1] & datetime<=times$end.time2[18])

temp<-read.csv("C:\\Users\\nt78066\\OneDrive - University of Georgia\\Documents\\Warming_uptake_stuff\\whole_stream_uptake\\time.series.temps.plot.csv")
temp$date.time<-ymd_hm(paste(temp$date,temp$time,sep=" "))

temp2<-aggregate(temp~date.time+stream,temp,mean)

vlines<-data.frame(lines=ymd_hms(paste(levels(as.factor(uptake.lengths$date)),"12:00:00",sep=" ")))

flow.temp<-ggplot(q.fixed,aes(x=datetime,y=discharge.ls,color=stream))+geom_line(size=1)+theme_classic()+
  geom_vline(data=vlines,aes(xintercept = as.numeric(lines)), linetype = 4, colour = "black")+
  geom_line(data=temp2,aes(x=date.time,y=temp,color=stream),lty=2,size=1)+
  theme(text = element_text(size=20),legend.position = c(0.28, 1.5))+
  scale_color_manual(values=c("steelblue3","steelblue3","firebrick2","firebrick2"),name="",labels=NULL)+
  xlab("Date")+ylab(expression("Discharge ("*L*s^-1*") or temperature ("^o*"C)"))

flow.temp

tiff(filename="whole_stream_gross_uptake.tiff",units="in",res=800,width=8,height=12,compression="lzw")
plot_grid(flow.temp,whole_stream_botlz,labels="AUTO",cols=1,label_x=0.9,label_y=0.93,label_size = 20)
dev.off()




######## estimating discharge based on release info to check calcs

bottom.cond<-cond[cond$meter==0,]
release.info$date<-ymd(release.info$date)

q.calcs<-merge(release.info,bottom.cond,by=c("date","stream"))

q.calcs$bottom.or.cl<-q.calcs$cor.cond*0.4967 ### Calibration of logger with Cl concentration
q.calcs$total.mg.cl<-q.calcs$nacl.mg*0.59*1000
q.calcs$conc.cl.mg.ml<-q.calcs$total.mg.cl/(q.calcs$water.volume.l*1000)
q.calcs$mg.cl.sec<-q.calcs$conc.cl.mg.ml*(q.calcs$release.rate.ml.min/60)
q.calcs$discharge.l.s<-q.calcs$mg.cl.sec/q.calcs$bottom.or.cl
q.calcs$ug.p.ml<-(0.228*q.calcs$mass.po4.g*10**6/q.calcs$water.volume.l)/1000
q.calcs$srp.ug.sec<-q.calcs$ug.p.ml*q.calcs$release.rate.ml.min/60
q.calcs$upstream.ug.l.p<-q.calcs$srp.ug.sec/q.calcs$discharge.l.s

d<-merge(q.calcs,uptake.lengths,by=c("date","stream"))

d$uptake.rate<-(d$discharge.l.s*1)/(d$width/100*d$sw)*60 ### estiamtes uptake rate in ug m2 min at a conc of 1 ug/L
d$log.rate<-log(d$uptake.rate)


ggplot(d,aes(x=q,y=discharge.l.s,color=stream))+geom_abline(intercept = 0,slope=1)+geom_point()+ylim(0,5)+xlim(0,5)+
  geom_smooth(method="lm",se=FALSE)




model3<-lm(log.rate~boltz.temp*stream,data=d)
anova(model3)
summary(model3)

model4<-lm(log.rate~boltz.temp+stream,data=d)
anova(model4)
summary(model4)

model5<-lm(log.rate~boltz.temp,data=d)
anova(model5)
summary(model5)
#### better model

fit.data<-data.frame(kt=seq(from=40.78272,to=41.52606,by=0.001))
fit.data$log.rate<-fit.data$kt*coef(model5)[2]+coef(model5)[1]

whole_stream_botlz<-ggplot(d,aes(x=boltz.temp,y=log.rate,color=stream))+theme_classic()+geom_point()+
  geom_line(data=fit.data,aes(x=kt,y=log.rate),color="black",size=3)+
  xlab("Inverse Temperature 1/kT")+
  ylab(expression("ln(SRP uptake rate ("*mu*g~"SRP hr"^-1~"m"^-2*"))"))+
  theme(text = element_text(size=20),legend.position = c(0.8, 0.9))+
  scale_color_manual(values=c("steelblue3","firebrick2"),name="")
whole_stream_botlz


AICc(model3,model4,model5)


tiff(filename="whole_stream_gross_uptake.tiff",units="in",res=800,width=8,height=12,compression="lzw")
plot_grid(flow.temp,whole_stream_botlz,labels="AUTO",cols=1,label_x=0.9,label_y=0.93,label_size = 20)
dev.off()

AIC.weights<-(AIC=AICc(model3,model4,model5))
AIC.weights$delta<-AIC.weights$AICc-min(AIC.weights$AICc)
AIC.weights$rel.likli<-exp(-0.5*AIC.weights$delta)
AIC.weights$relative.weights<-AIC.weights$rel.likli/sum(AIC.weights$rel.likli)
