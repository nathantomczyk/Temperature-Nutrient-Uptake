##################################################################################
### This code summarizes the results of a lab experiment where we measured uptake and respiraiton rates
### on conditioned leaf litter at different temperatures. We estimate the Ea of both uptake and
### respiraiton, and compare the molar ratios at which these two were used across temperatures. 


############ Load libraries 

library(ggplot2)
library(tidyverse)
library(lubridate)
library(lme4)
library(lmerTest)
library(MuMIn)
library(Rmisc)

## Set working directory
setwd("C:\\Users\\nt78066\\OneDrive - University of Georgia\\Documents\\Warming_uptake_stuff\\Lab_data\\")

############ loading data
mass<-read.csv("Lab_leaf_masses_spring2021.csv")

resp<-read.csv("Lab_respiration_spring2021.csv")

volumes<-read.csv("vial_volumes_lab_experiment.csv")

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

####### Estimating Ea of repiraiton
model.all.data<-lmer(log.rate~measured.boltz.temp+(1|bottle.id),data=resp3)
anova(model.all.data)
model.all.data.sum<-summary(model.all.data)
r.squaredGLMM(model.all.data)
model.all.data.sum


plot(model.all.data) # residuals look alright

## line of best fit
fit.data1<-data.frame(kt=seq(from=39.80249,to=42.04423,by=0.0001))
fit.data1$y<-fit.data1$kt*model.all.data.sum$coefficients[2,1]+model.all.data.sum$coefficients[1,1]


resp.summary<-summarySE(resp3,measurevar = "log.rate",groupvars = c("measured.boltz.temp","batch","bottle.id"),na.rm=TRUE)

respiration.plot<-ggplot(resp.summary,aes(x=measured.boltz.temp,y=log.rate))+
  geom_line(data=fit.data1,aes(x=kt,y=y),size=2,color="blue")+
  geom_point(size=3)+
  theme_classic()+theme(text = element_text(size=20))+geom_errorbar(aes(ymin=log.rate-se,ymax=log.rate+se),width=0.01)+
  xlab("Inverse Temperature 1/kT")+
  ylab(expression("ln(Respiration rate ("*"mg"~O[2]~"hr"^-1~"g"^-1*"))"))
  #geom_line(data=fit.data1,aes(x=kt,y=y),size=2,color="black")+

respiration.plot

##############################################################################################
######################## Uptake ##############################################################
##############################################################################################

# load data
uptake.times<-read.csv("Lab_uptake_spring2021.csv")
srp<-read.csv("lab_uptake_srp_data.csv")


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
uptake.summary<-summarySE(uptake.data.3,measurevar = "log.rate",groupvars = c("measured.boltz.temp","batch","Bottle.Letter.x"),na.rm=TRUE)

uptake.data.3$bottle.id<-paste(uptake.data.3$Temperature.x,uptake.data.3$Bottle.Letter,uptake.data.3$batch,sep="-")

## Estimting Ea of SRP uptake
uptake.model<-lmer(log.rate~measured.boltz.temp+(1|bottle.id),data=uptake.data.3,REML = TRUE)
uptake.model.sum<-summary(uptake.model)
r.squaredGLMM(uptake.model)
uptake.model.sum
plot(uptake.model) # residuals look alright 

fit.data4<-data.frame(kt=seq(from=39.52777,to=41.92,by=0.0001))
fit.data4$y<-fit.data4$kt*uptake.model.sum$coefficients[2,1]+uptake.model.sum$coefficients[1,1]

#position_dodge(width=0.1)

uptake_plot<-ggplot(uptake.summary[uptake.summary$batch!=1,],aes(x=measured.boltz.temp,y=log.rate))+
  geom_line(data=fit.data4,aes(x=kt,y=y),size=2,color="blue")+
  geom_point(size=2,position=position_dodge(width = .5))+theme_classic()+theme(text=element_text(size=20))+
  geom_errorbar(aes(ymin=log.rate-se,ymax=log.rate+se),position=position_dodge(width = .5),width=0.0001)+
  xlab("Inverse Temperature 1/kT")+
  ylab(expression("ln(SRP uptake rate ("*mu*g~"hr"^-1~"g"^-1*"))"))
  
uptake_plot

################################################ Are the slopes actually different???

uptake_data_short<-data.frame(boltzman.temp=uptake.data.3$measured.boltz.temp,batch=uptake.data.3$batch,
                              temp.treatment=uptake.data.3$Temperature.x,letter=uptake.data.3$Bottle.Letter,
                              rate=uptake.data.3$log.rate,measurement="uptake")
resp_data_short<-data.frame(boltzman.temp=resp3$measured.boltz.temp,batch=resp3$batch,temp.treatment=resp3$temperature,
                            letter=resp3$Bottle.Letter,rate=resp3$log.rate,measurement="respiration")

slope.test.data<-rbind(uptake_data_short,resp_data_short)
slope.test.data$measurement<-as.factor(slope.test.data$measurement)

## are the slopes different

difference.model<-lm(rate~boltzman.temp*measurement,data=slope.test.data)
anova(difference.model)
difference.model.sum<-summary(difference.model)
difference.model.sum
r.squaredGLMM(difference.model)



###################### looking at temperature effect on molar ratios
uptake_data_short<-data.frame(boltzman.temp=uptake.data.3$measured.boltz.temp,batch=uptake.data.3$batch,
                              temp.treatment=uptake.data.3$Temperature.x,letter=uptake.data.3$Bottle.Letter,
                              rate=uptake.data.3$srp.uptake.ug.srp.g.hr,measurement="uptake",actual.temp=uptake.data.3$measured_temperature)
resp_data_short<-data.frame(boltzman.temp=resp3$measured.boltz.temp,batch=resp3$batch,temp.treatment=resp3$temperature,
                            letter=resp3$Bottle.Letter,rate=resp3$DO.consumed.mg.g.dm,measurement="respiration")



resp_data_short$resp.rate.molar<-resp_data_short$rate/32/1000 ## converting to mols of C hr-1 g-1

uptake_data_short$uptake.rate.molar<-uptake_data_short$rate/1000/1000/31 ## converting to mols of P hr-1 g-1

uptake.short.summary<-summarySE(uptake_data_short,measurevar = "uptake.rate.molar",groupvars=c("batch","temp.treatment","letter"))
resp.short.summary<-summarySE(resp_data_short,measurevar = "resp.rate.molar",groupvars=c("batch","temp.treatment","letter"))

stoich.data<-merge(uptake.short.summary,resp.short.summary,by=c("batch","temp.treatment","letter"))

stoich.data$use.ratio<-stoich.data$resp.rate.molar/stoich.data$uptake.rate.molar

temps<-aggregate(actual.temp~letter+temp.treatment+batch,data=uptake_data_short,mean)

stoich.data.2<-merge(stoich.data,temps,by=c("letter","temp.treatment","batch"))


stoich.model<-lmer(log(use.ratio)~actual.temp+(1|batch),data=stoich.data.2)
anova(stoich.model)
stoich.mod.sum<-summary(stoich.model)
r.squaredGLMM(stoich.model)

stoich.model2<-lmer(use.ratio~actual.temp+(1|batch),data=stoich.data.2)
anova(stoich.model2)
stoich.mod.sum2<-summary(stoich.model2)
r.squaredGLMM(stoich.model2)


fit.1<-data.frame(temp=seq(from=4,to=20,by=0.1))
fit.1$y<-exp(fit.1$temp*stoich.mod.sum$coefficients[2,1]+stoich.mod.sum$coefficients[1,1])



stoich_plot<-ggplot(stoich.data.2[stoich.data.2$batch!=1,],aes(x=actual.temp,y=use.ratio))+geom_point(size=2)+
  scale_y_log10()+
  theme_classic()+
  geom_line(data=fit.1,aes(x=temp,y=y),color="blue",size=2)+
  theme(text=element_text(size=20))+
  xlab(expression("Temperature "^o*"C"))+
  ylab("C:P respiration to uptake ratio")
stoich_plot
########################

################ Plotting data

tiff(filename="lab_experiment_data.tiff",units="in",res=800,width=8,height=16,compression="lzw")
plot_grid(respiration.plot,uptake_plot,stoich_plot,labels="AUTO",cols=1,label_x=0.9,label_y=0.93,label_size = 20)
dev.off()

#################
