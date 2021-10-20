############ Notes
# I have the code worked up to estimate uptake in the streams and the R2 look pretty decent for all estimates
# Need to check out the discharge estiamtes and figure out what is up with the holes in the data
# Also want to check with phillip about how is was/ wasn't doing stage correction of the discharge
# Once these issues are worked out I *THINK* these estimates may be good to go. 
############

# load packages

library(dplyr)
library(tidyverse)
library(lubridate)
library(cowplot)
library(ggplot2)
library(MuMIn)

#load data

q<-read.csv("./Whole_stream_gross_uptake/crews_paireddischarge_all_23Aug2021_long.csv")

release.info<-read.csv("./Whole_stream_gross_uptake/whole.stream.uptake.times.temps.release.info.csv")

chem<-read.csv("./Whole_stream_gross_uptake/uptake_temperature_nutrients.csv")

chem[chem$sample.id==10 & chem$date=="2021-01-26","po4.p.ug.L"]<-NA

widths<-read.csv("./Whole_stream_gross_uptake/stream_widths.csv")

cond<-read.csv("./Whole_stream_gross_uptake/whole_stream_conductivity.csv")

#format dates
chem$date<-ymd(chem$date)
widths$date<-ymd(widths$date)
cond$date<-ymd(cond$date)
q$datetime<-ymd_hms(q$datetime)
release.info$date<-ymd(release.info$date)

# code factors (stream names)
q$stream<-as.factor(q$stream)
q$stream<-recode_factor(q$stream,'ws55'="ref",'towr'="warm")

# correctiong background conductivity
cond$cor.cond<-(cond$plt.spc-cond$background.spc)

cond2<-data.frame(date=cond$date,stream=cond$stream,corrected.cond=cond$cor.cond,meter=cond$meter)

# calculating discharge from salt rlease

bottom.cond<-cond2[cond2$meter==0,]

q.calcs<-merge(release.info,bottom.cond,by=c("date","stream"))

q.calcs$bottom.or.cl<-q.calcs$corrected.cond*0.4967 ### Calibration of logger with Cl concentration
q.calcs$total.mg.cl<-q.calcs$nacl.mg*0.59*1000
q.calcs$conc.cl.mg.ml<-q.calcs$total.mg.cl/(q.calcs$water.volume.l*1000)
q.calcs$mg.cl.sec<-q.calcs$conc.cl.mg.ml*(q.calcs$release.rate.ml.min/60)
q.calcs$discharge.l.s<-q.calcs$mg.cl.sec/q.calcs$bottom.or.cl
q.calcs$ug.p.ml<-(0.228*q.calcs$mass.po4.g*10**6/q.calcs$water.volume.l)/1000
q.calcs$srp.ug.sec<-q.calcs$ug.p.ml*q.calcs$release.rate.ml.min/60
q.calcs$upstream.ug.l.p<-q.calcs$srp.ug.sec/q.calcs$discharge.l.s

q.calcs.short<-data.frame(date=q.calcs$date,stream=q.calcs$stream,
                          top.of.reach.srp=q.calcs$upstream.ug.l.p,discharge.l.s=q.calcs$discharge.l.s,
                          temperature=q.calcs$temperature)

### aggregating water chemistry data
chem2<-aggregate(po4.p.ug.L~date+stream+meter+bg.plt,chem,mean)

plt.only<-chem2[chem2$bg.plt=="plt",]

bg.only<-chem[chem$bg.plt=="bkd",]

names(bg.only)[names(bg.only) == "po4.p.ug.L"] <- "background.PO4"
bg.only$id<-paste(bg.only$date,bg.only$stream,sep="-")

ggplot(bg.only[bg.only$stream=="ref",],aes(x=meter,y=background.PO4,color=id))+geom_point()+geom_smooth(se=FALSE)

chem.bg<-merge(plt.only,bg.only,by=c("stream","meter","date"))

chem3<-merge(chem.bg,cond2,by=c("date","meter","stream"))
chem4<-merge(chem3,q.calcs.short,by=c("date","stream"))

##### percent dilution estimates #######

chem4$stream.date<-as.factor(paste(chem4$stream,chem4$date,sep="-"))

id<-levels(chem4$stream.date)

chem5<-data.frame()


for (i in 1:length(id)){
  s<-chem4[chem4$stream.date==id[i],]
  
  eighty<-s[s$meter==80,]
  
  s$dilution<-s$corrected.cond/eighty$corrected.cond
  

  
  chem5<-rbind(chem5,s)
  
}

############## Calculating uptake lengths

chem5$dist.from.release<-100-chem5$meter

#date.stream<-unique(data.frame(date=chem4$date,stream=chem4$stream))

uptake.lengths<-data.frame()

#chem4[chem4$date==ymd("2021-01-26") & chem4$stream=="ref" & chem4$dist.from.eighty==60,"bg.corrected.srp"]<-NA
#chem4[chem4$date==ymd("2021-01-26") & chem4$stream=="ref" & chem4$dist.from.eighty==60,"normalized.srp"]<-NA

setwd("./figures/")

for (i in 1:length(id)){
  s<-chem5[chem5$stream.date==id[i],]
  
  #bg<-chem4[chem4$stream==s$stream[1],]

  s$bg.corrected.srp<-s$po4.p.ug.L -mean(s$background.PO4,na.rm=TRUE)
  
  s$dilution.corrected.srp<-s$bg.corrected.srp/s$dilution
  
  s$srp.remain.proportion<-s$dilution.corrected.srp/s$top.of.reach.srp
  

  #model<-lm(log(dilution.corrected.srp)~dist.from.eighty,s)
  
  model<-lm(log(srp.remain.proportion)~dist.from.release+0,s)
  model.sum<-summary(model)
  
  
  uptake.plot<- ggplot(s,aes(x=dist.from.release,y=log(srp.remain.proportion)))+geom_point()+theme_classic()+
    geom_abline(intrecept=0,slope=model.sum$coefficients[1])+xlim(0,100)+#ylim(0,log(.0001))+
    xlab("distance")+ylab("Normalized srp concentraiton")+
    annotate(x=80,y=log(0.8),"text",label=paste("R2=",round(model.sum$r.squared,2),sep=" "),size=6)+
    annotate(x=80,y=log(0.6),"text",label=paste(id[i]),size=6)
    #scale_y_log10()
  
  uptake.plot
  
  
  tiff(filename=paste(id[i],".tiff",sep=""),units="in",res=200,width=8,height=8,compression="lzw")
  print(uptake.plot)
  dev.off()
    
  output<-data.frame(stream=s$stream[1],date=s$date[1],k=model.sum$coefficients[1],
                     sw=(-1/model.sum$coefficients[1]),r2=model.sum$r.squared,discharge.l.s=s$discharge.l.s[1],
                     temperature=s$temperature[1])
  uptake.lengths<-rbind(uptake.lengths,output)
}



  
ggplot(uptake.lengths,aes(x=date,y=sw,color=stream))+geom_point()+geom_line()+theme_classic()  





#####

background.conc<-aggregate(po4.p.ug.L~stream+date,chem[chem$bg.plt=="bkd",],mean)

widths2<-aggregate(width.cm~stream+date,widths,mean)

uptake.lengths<-merge(uptake.lengths,background.conc,by=c("stream","date"))

uptake.lengths<-merge(uptake.lengths,widths2,by=c("stream","date"))

uptake.lengths$uptake.rate<-(uptake.lengths$discharge.l.s*1)/(uptake.lengths$width.cm/100*uptake.lengths$sw)*60 ### estiamtes uptake rate in ug m2 min at a conc of 1 ug/L

ggplot(uptake.lengths,aes(x=date,y=uptake.rate,color=stream))+geom_point()+geom_line()+theme_classic()  

ggplot(uptake.lengths,aes(x=temperature,y=uptake.rate,color=stream))+geom_point()+geom_smooth(method="lm",se=FALSE)+theme_classic()  

Boltz<-8.6173303*10**-5 # boltzmann constant

uptake.lengths$boltz.temp<-1/((273+uptake.lengths$temperature)*Boltz)

uptake.lengths$log.rate<-log(uptake.lengths$uptake.rate)




model3<-lm(log.rate~boltz.temp*stream,data=uptake.lengths)
anova(model3)
summary(model3)

model4<-lm(log.rate~boltz.temp+stream,data=uptake.lengths)
anova(model4)
summary(model4)

model5<-lm(log.rate~boltz.temp,data=uptake.lengths)
anova(model5)
summary(model5)

AIC.weights<-(AIC=AICc(model3,model4,model5))
AIC.weights$delta<-AIC.weights$AICc-min(AIC.weights$AICc)
AIC.weights$rel.likli<-exp(-0.5*AIC.weights$delta)
AIC.weights$relative.weights<-AIC.weights$rel.likli/sum(AIC.weights$rel.likli)




fit.data<-data.frame(kt=seq(from=40.78272,to=41.52606,by=0.001))
fit.data$log.rate<-fit.data$kt*coef(model5)[2]+coef(model5)[1]

whole_stream_botlz<-ggplot(uptake.lengths,aes(x=boltz.temp,y=log.rate,color=stream))+theme_classic()+
  geom_line(data=fit.data,aes(x=kt,y=log.rate),color="black",size=2)+
  xlab("Inverse Temperature 1/kT")+
  ylab(expression("ln(SRP uptake rate ("*mu*g~"SRP hr"^-1~"m"^-2*"))"))+
  theme(text = element_text(size=20),legend.position = c(0.8, 0.9))+
  geom_point(size=2)+
  scale_color_manual(values=c("steelblue3","firebrick2"),name="")
whole_stream_botlz

####### Plotting discharge and temperature during the experiment


############### Depth logger in the warm stream had some issues during the deployment period - fixing data by
############### regressing with the warm stream from a time period where things worked better

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

temp<-read.csv("./Whole_stream_gross_uptake/time.series.temps.plot.csv")
temp$date.time<-ymd_hm(paste(temp$date,temp$time,sep=" "))


needed.dates<- q %>% filter(datetime>=ymd("2021-01-25") & datetime<=ymd("2021-03-10") & stream=="ref")

warm.discharge<-data.frame(X=NA,datetime=needed.dates$datetime,year=2021,discharge.ls=coef(discharge.model)[1]+coef(discharge.model)[2]*needed.dates$discharge.ls,water_temp.c=NA,stream="warm")

q.fixed<-rbind(needed.dates,warm.discharge)

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
plot_grid(flow.temp,whole_stream_botlz,labels="AUTO",ncol=1,label_x=0.9,label_y=0.93,label_size = 20)
dev.off()

