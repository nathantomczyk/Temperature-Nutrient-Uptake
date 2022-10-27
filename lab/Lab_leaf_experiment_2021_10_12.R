##################################################################################
### This code summarizes the results of a lab experiment where we measured uptake and respiraiton rates
### on conditioned leaf litter at different temperatures. We estimate the Ea of both uptake and
### respiration, and compare the molar ratios at which these two were used across temperatures. 
##################################################################################

############ Load libraries 

library(ggplot2)
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(MuMIn)
library(Rmisc)
library(cowplot)




############ loading data
mass<-read.csv("./lab/Lab_leaf_masses_spring2021.csv")

resp<-read.csv("./lab/Lab_respiration_spring2021.csv")

volumes<-read.csv("./lab/vial_volumes_lab_experiment.csv")

##############################################################################################
######################## Resipration #########################################################
##############################################################################################


Boltz<-8.6173303*10**-5 ## boltzmann constant

## converting time format to decimal
resp$time.elapsed2<-as.numeric(substring(resp$time.elapsed,1,1))+as.numeric(substring(resp$time.elapsed,3,4))/60
# merging exact volume of vials used
resp2<-merge(resp,volumes,by="Sample.number")

## mass of dissolved oxygen consumed during the measurement
resp2$DO.consumed.mg<-resp2$start.DO*(resp2$vial.volume.ml/1000)-resp2$end.DO*(resp2$vial.volume.ml/1000)
resp2$DO.consumed.mg.hr<-resp2$DO.consumed.mg/resp2$time.elapsed2

#averages for the vials that had no leaves added(i.e., the blanks)
blanks<-resp2[resp2$bottle=="blank ",]
blanks2<-data.frame(batch=blanks$batch,temperature=blanks$temperature,blank.DO.consumed=blanks$DO.consumed.mg.hr)
blanks3<-aggregate(blank.DO.consumed~batch+temperature,data=blanks2,mean)


samples<-resp2[resp2$bottle!="blank ",]

# merging values for the blanks into the measurements of samples. 
samples2<-merge(samples,blanks3,by=c("temperature","batch"))

# correcting for loss in the blanks
samples2$DO.consumed.mg.corrected<-samples2$DO.consumed.mg.hr-samples2$blank.DO.consumed

#merging in the mass of leaves used in each measurement
resp3<-merge(samples2,mass[mass$Measurement=="Respiration",],by=c("batch","Sample.number"))

#normalizign rates of mass loss for the mass of leaves used
resp3$DO.consumed.mg.g.dm<-resp3$DO.consumed.mg.corrected/resp3$Dry.Mass..g.

######## calculating boltzmann temperature and log transforming rates of respiration
resp3$measured.boltz.temp<-1/(Boltz*(273+resp3$measured_temperature))
resp3$log.rate<-log(resp3$DO.consumed.mg.g.dm) ### this drops one observation at 4C where DO went up marginally - likely a bubble in the vial

resp3$bottle.id<-paste(resp3$batch,resp3$bottle,resp3$temperature,sep="-")

summarized.resp<-aggregate(DO.consumed.mg.g.dm~measured.boltz.temp+bottle.id+batch+bottle+measured_temperature,data=resp3,mean)
summarized.resp$log.rate<-log(summarized.resp$DO.consumed.mg.g.dm)


summarized.resp$centered.boltz.temp<-1/(Boltz*(273+12))-summarized.resp$measured.boltz.temp

####### Estimating Ea of respiraiton
model.all.data<-lm(log.rate~centered.boltz.temp+as.factor(batch),data=summarized.resp)
anova(model.all.data)
model.all.data.sum<-summary(model.all.data)
r.squaredGLMM(model.all.data)
model.all.data.sum


plot(model.all.data) # residuals look alright

label.position<-0.02
label.position.r<-0.06

label.pos.resp<-(range(summarized.resp$DO.consumed.mg.g.dm)[2]-range(summarized.resp$DO.consumed.mg.g.dm)[1])*label.position+range(summarized.resp$DO.consumed.mg.g.dm)[1]
label.pos.r<-(range(summarized.resp$DO.consumed.mg.g.dm)[2]-range(summarized.resp$DO.consumed.mg.g.dm)[1])*label.position.r+range(summarized.resp$DO.consumed.mg.g.dm)[1]


respiration.plot<-ggplot(summarized.resp,aes(x=measured_temperature,y=DO.consumed.mg.g.dm))+
  geom_smooth(method="lm",se=FALSE,size=2,color="grey60")+
  geom_point(size=3)+
  theme_classic()+theme(text = element_text(size=20))+
  #xlab(expression("Inverse Temperature (1/kT - 1/kT"[0]*")"))+
  #scale_x_continuous(sec.axis=~1/(Boltz*(273+12))-1/((.+273)*Boltz),name="new")+
  scale_x_continuous(sec.axis=sec_axis(~1/(Boltz*(273+12))-1/((.+273)*Boltz),
                                       name=expression("Standardized temperature "*"(1/["*italic("kT")[12]-italic(kT) *"])")))+
  #xlab(expression("Temperature "^o*"C"))+
  xlab("")+
  ylab(expression("Respiration rate (mg O"[2]~"hr"^-1~"g"^-1*")"))+
  scale_y_log10()+
  annotate("text",18,label.pos.resp,label=expression(italic(Ea)*" = 1.02 "*italic("eV")),size=8)+
  annotate("text",18,label.pos.r,label=expression(italic(R)^2*" = 0.92"),size=8)

respiration.plot



##############################################################################################
######################## Uptake ##############################################################
##############################################################################################

# load data
uptake.times<-read.csv("./lab/Lab_uptake_spring2021.csv")
srp<-read.csv("./lab/lab_uptake_srp_data.csv")


######### one outlier

srp[srp$Sample.number==29 & srp$batch==3,"srp.ug.l"]<-NA  # This observation is really close to the blanks and is much different than the other samples in this treatment

# incubation times
uptake.times$time.elapsed<-hm(uptake.times$End.Time)-hm(uptake.times$Start.Time)

uptake.data<-merge(srp,uptake.times,by=c("batch","Sample.number"))

# Loss of SRP in blanks

uptake.data.blanks<-uptake.data[uptake.data$Bottle.Letter.x=="blank",]
uptake.data.blanks2<-data.frame(date=uptake.data.blanks$date,Temperature.x=uptake.data.blanks$Temperature.x,blank.srp=uptake.data.blanks$srp.ug.l)

mean.blanks<-aggregate(blank.srp~Temperature.x*date,uptake.data.blanks2,mean)

# merge the blanks back in
uptake.data.2<-merge(uptake.data[uptake.data$Bottle.Letter.x!="blank",],mean.blanks,by=c("date","Temperature.x"))

uptake.data.2$srp.loss<-uptake.data.2$blank.srp*0.04-uptake.data.2$srp.ug.l*0.04 ## division converts concentraiton to mass (40mL of water in chambers)

## converting time format to decimal
uptake.data.2$time.elapsed3<-uptake.data.2$time.elapsed@hour + uptake.data.2$time.elapsed@minute / 60


uptake.data.3<-merge(uptake.data.2,mass[mass$Measurement=="P uptake",],by=c("batch","Sample.number")) # merging lea fragment mass in

uptake.data.3$srp.uptake.ug.srp.g.hr<-uptake.data.3$srp.loss/uptake.data.3$Dry.Mass..g./uptake.data.3$time.elapsed3

## Transforming rates and temperatures
uptake.data.3$log.rate<-log(uptake.data.3$srp.uptake.ug.srp.g.hr)
uptake.data.3$measured.boltz.temp<-1/(Boltz*(273+uptake.data.3$measured_temperature))


## summarizing data for plot

uptake.data.3$bottle.id<-paste(uptake.data.3$batch,uptake.data.3$Bottle.Letter,uptake.data.3$Temperature.x,sep="-")


uptake.summary<-aggregate(srp.uptake.ug.srp.g.hr~measured.boltz.temp+batch+bottle.id+Bottle.Letter+measured_temperature,data=uptake.data.3,mean)
uptake.summary$log.rate<-log(uptake.summary$srp.uptake.ug.srp.g.hr)

uptake.summary$centered.boltz.temp<-1/(Boltz*(273+12))-uptake.summary$measured.boltz.temp

## Estimting Ea of SRP uptake
uptake.model<-lm(log.rate~centered.boltz.temp+batch,data=uptake.summary)
uptake.model.sum<-summary(uptake.model)
uptake.model.sum
anova(uptake.model)
plot(uptake.model) # residuals look alright 

label.pos.uptake<-(range(uptake.summary$srp.uptake.ug.srp.g.hr)[2]-range(uptake.summary$srp.uptake.ug.srp.g.hr)[1])*label.position+range(uptake.summary$srp.uptake.ug.srp.g.hr)[1]
label.pos.uptake.r<-(range(uptake.summary$srp.uptake.ug.srp.g.hr)[2]-range(uptake.summary$srp.uptake.ug.srp.g.hr)[1])*label.position.r+range(uptake.summary$srp.uptake.ug.srp.g.hr)[1]


uptake_plot<-ggplot(uptake.summary,aes(x=measured_temperature,y=srp.uptake.ug.srp.g.hr))+
  geom_smooth(method="lm",se=FALSE,size=2,color="grey60")+
  geom_point(size=3)+theme_classic()+theme(text=element_text(size=20))+
  #xlab(expression("Inverse Temperature (1/kT - 1/kT"[0]*")"))+
  #xlab(expression("Temperature "^o*"C"))+
  ylab(expression("SRP uptake rate ("*mu*g~"hr"^-1~"g"^-1*")"))+
  xlab("")+
  scale_x_continuous(sec.axis=sec_axis(~1/(Boltz*(273+12))-1/((.+273)*Boltz),labels=NULL))+
  #scale_x_continuous(sec.axis=sec_axis(~1/(Boltz*(273+12))-1/((.+273)*Boltz),
  #                                     name=expression("Standardized temperature "*"(1/["*italic("kT")[12]-italic(kT) *"])")))+
  
  scale_y_log10()+
  #xlim(-1.13,1.23)+
  #scale_y_continuous(breaks=c(1,0))+
  annotate("text",18,label.pos.uptake,label=expression(italic(Ea)*" = 0.48 "*italic("eV")),size=8)+
  annotate("text",18,label.pos.uptake.r,label=expression(italic(R)^2*" = 0.81"),size=8)
uptake_plot


  


################################################ Are the slopes actually different???

uptake_data_short<-data.frame(boltzman.temp=uptake.summary$measured.boltz.temp,rate=uptake.summary$log.rate,
                              measurement="uptake",batch=uptake.summary$batch,id=uptake.summary$bottle.id)
resp_data_short<-data.frame(boltzman.temp=summarized.resp$measured.boltz.temp,rate=summarized.resp$log.rate
                            ,measurement="respiration",batch=summarized.resp$batch,id=summarized.resp$bottle.id)

slope.test.data<-rbind(uptake_data_short,resp_data_short)
slope.test.data$measurement<-as.factor(slope.test.data$measurement)

## are the slopes different

difference.model<-lm(rate~boltzman.temp*measurement+as.factor(batch),data=slope.test.data)
anova(difference.model)
difference.model.sum<-summary(difference.model)
difference.model.sum
r.squaredGLMM(difference.model)



###################### looking at temperature effect on molar ratios
uptake_data_short<-data.frame(batch=uptake.summary$batch,
                              rate=uptake.summary$srp.uptake.ug.srp.g.hr,measurement="uptake",
                              letter=uptake.summary$Bottle.Letter,temperature=uptake.summary$measured_temperature,
                              id=uptake.summary$bottle.id)
resp_data_short<-data.frame(boltzman.temp=summarized.resp$measured.boltz.temp,batch=summarized.resp$batch,
                            rate=summarized.resp$DO.consumed.mg.g.dm,measurement="respiration",
                            letter=summarized.resp$bottle,temperature=summarized.resp$measured_temperature,
                            id=summarized.resp$bottle.id)



resp_data_short$resp.rate.molar<-resp_data_short$rate/32*0.85 ## converting to mmols of C hr-1 g-1 0.85 is respirtory quotient from Berggen et al 2012

uptake_data_short$uptake.rate.molar<-uptake_data_short$rate/1000/31 ## converting to mmols of P hr-1 g-1


stoich.data<-merge(resp_data_short,uptake_data_short,by=c("id"))

stoich.data$use.ratio<-stoich.data$resp.rate.molar/stoich.data$uptake.rate.molar

#temps<-aggregate(actual.temp~letter+temp.treatment+batch,data=uptake_data_short,mean)

#stoich.data.2<-merge(stoich.data,temps,by=c("letter","temp.treatment","batch"))


stoich.model<-lm(log(use.ratio)~boltzman.temp+batch.x,data=stoich.data)
anova(stoich.model)
stoich.mod.sum<-summary(stoich.model)
stoich.mod.sum

stoich.data$log.ratio<-log(stoich.data$use.ratio)

stoich.data$centered.botlz.temp<-stoich.data$boltzman.temp-1/((273+12)*Boltz)

stoich.label.pos<-(range(stoich.data$use.ratio)[2]-range(stoich.data$use.ratio)[1])*label.position+range(stoich.data$use.ratio)[1]
stoich.label.pos.r<-(range(stoich.data$use.ratio)[2]-range(stoich.data$use.ratio)[1])*label.position.r+range(stoich.data$use.ratio)[1]


stoich_plot<-ggplot(stoich.data,aes(x=temperature.x,y=use.ratio))+geom_point(size=3)+
  theme_classic()+
  geom_smooth(method="lm",se=FALSE,size=2,color="grey60")+
  theme(text=element_text(size=20))+
  xlab(expression("Temperature ("^o*"C)"))+
  ylab("C:P respiration to uptake")+
  scale_x_continuous(sec.axis=sec_axis(~1/(Boltz*(273+12))-1/((.+273)*Boltz),labels=NULL))+
  #scale_x_continuous(sec.axis=sec_axis(~1/(Boltz*(273+12))-1/((.+273)*Boltz),
  #                                     name=expression("Standardized temperature "*"(1/["*italic("kT")[12]-italic(kT) *"])")))+
  scale_y_log10()+
  annotate("text",18,stoich.label.pos,label=expression(italic(Ea)*" = 0.54 "*italic("eV")),size=8)+
  annotate("text",18,stoich.label.pos.r,label=expression(italic(R)^2*" = 0.56"),size=8)
stoich_plot



  xlab(expression("Temperature "^o*"C"))+ylab(expression("Respiration rate (mg O"[2]~"hr"^-1~"g"^-1*")"))+
  


stoich.data$untransformed.temp<-stoich.data$boltzman.temp**-1/Boltz-273
stoich.model<-lm(use.ratio~temperature.x+batch.x,data=stoich.data)
summary(stoich.model)
anova(stoich.model)
########################
tiff(filename="./figures/resp_data_poster.tiff",units="in",res=800,width=12,height=8,compression="lzw")
respiration.plot
dev.off()

tiff(filename="./figures/p_uptake_data.tiff",units="in",res=800,width=12,height=8,compression="lzw")
uptake_plot
dev.off()

tiff(filename="./figures/stoich_data.tiff",units="in",res=800,width=12,height=8,compression="lzw")
stoich_plot
dev.off()


################ Plotting data

tiff(filename="./figures/lab_experiment_data_feb2022.tiff",units="in",res=800,width=8,height=16,compression="lzw")
plot_grid(respiration.plot,uptake_plot,stoich_plot,labels=c("(a)","(b)","(c)"),ncol=1,label_x=0.12,label_y=0.94,label_size = 20,align="v")
dev.off()

#################
