library("ggplot2")
library("Rmisc")
library("cowplot")

Boltz<-8.6173303*10**-5 ## boltzmann constant

PO4uptake_data<-read.csv("./lab2/Concentration_uptake_data.csv", header=TRUE)

#
PO4uptake_data$PO4.ug.l.dilution<-PO4uptake_data$PO4.ug.l*PO4uptake_data$dilution # correcting for dilution at the analysis stage

#seperating and specifying final and initial samples
blanks<-PO4uptake_data[which(PO4uptake_data$sample=="blank"),]

#are pre and post blanks similar?
blank.test<-lm(PO4.ug.l.dilution~concentraiton.treatment*group,data=blanks)
anova(blank.test)
# pre and post blanks are not different so they will be averaged for analysis
blanks2<-aggregate(PO4.ug.l.dilution~temp.treatment+concentraiton.treatment,data=blanks,mean)
blanks3<-data.frame(temp.treatment=blanks2$temp.treatment,concentraiton.treatment=blanks2$concentraiton.treatment,blank.concentraiton.srp=blanks2$PO4.ug.l.dilution)

#
samples<-PO4uptake_data[which(PO4uptake_data$sample=="final"),]

#making new dataframe with only the critical data from the blanks

#recombining dataframes into one, defining conc change and then loss rate per leaf mass
d<-merge(samples,blanks3, by=c("temp.treatment", "concentraiton.treatment"))
d$srp.change<-d$blank.concentraiton.srp*0.2-d$PO4.ug.l.dilution*0.2 ## 0.2 Liters of water in each incubation
d$srp.loss.rate.ug.g.hr<-d$srp.change/d$time.elapsed/d$leaf_dry_mass_g

## renaming target concentrations based on actual measurements (for plotting purposes)

new.concentrations<-round(aggregate(blank.concentraiton.srp~concentraiton.treatment,data=d,mean),0)

names(new.concentrations)[names(new.concentrations) == "blank.concentraiton.srp"] <- "actual.initial.conc.ug.l"


d$log.rate<-log(d$srp.loss.rate.ug.g.hr)  ## two samples had negative rates of SRP loss
d$temp.kt<-1/(Boltz*(273+d$temp.measured.c))

d2<-merge(d,new.concentrations,by="concentraiton.treatment")

d2$centered.boltz.temp<-d2$temp.kt-1/((273+12)*Boltz)

custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

opposite.plot<-ggplot(d2,aes(x=temp.measured.c,y=srp.loss.rate.ug.g.hr,color=as.factor(actual.initial.conc.ug.l)))+theme_classic()+geom_point(size=3)+
  geom_smooth(method="lm",se=FALSE,size=2)+
  theme(text = element_text(size=20),legend.position = c(0.24, 0.18))+
  scale_colour_brewer(palette="Dark2",direction=1,name = expression(PO[4]~mu*g~L^-1))+
  scale_x_continuous(sec.axis=sec_axis(~1/(Boltz*(273+12))-1/((.+273)*Boltz),
                                       name=expression("Standardized temperature "*"(1/["*italic("kT")[12]-italic(kT) *"])")))+
  xlab(expression("Temperature ("^o*"C)"))+
  scale_y_log10()+
  guides(color=guide_legend(ncol=2))+
  ylab(expression("SRP uptake rate ("*mu*g~"hr"^-1~"g"^-1*")"))
opposite.plot

d2$inverse.rate<-1/d2$srp.loss.rate.ug.g.hr
d2$inverse.conc<-1/d2$actual.initial.conc.ug.l
temps<-levels(as.factor(d2$temp.treatment))

d2<-d2[d2$abs!=0.289,] ### one point from the 4C temp needed to be removed to let the model fit

mm_fits<-data.frame(concentraiton=rep(1:225,times=6),temperature=rep(temps,each=225),model.fit=NA)

mm_paramaters<-data.frame(km=NA,km.se=NA,vmax=NA,vmax.se=NA,temperature=temps,km.p=NA,vmax.p=NA,rsq=NA)

for (i in 1:length(temps)){
  one_temp<-d2[d2$temp.treatment==temps[i],]

  model.nls <- nls(srp.loss.rate.ug.g.hr ~ Vm * actual.initial.conc.ug.l / (K + actual.initial.conc.ug.l),one_temp, start=list(Vm=30, K=25))
  a<-summary(model.nls)
  
  mm_fits[mm_fits$temperature==temps[i],"model.fit"]<-mm_fits[mm_fits$temperature==temps[i],"concentraiton"]*coef(model.nls)[1]/(mm_fits[mm_fits$temperature==temps[i],"concentraiton"]+coef(model.nls)[2])
  
  RSS<-sum(residuals(model.nls)^2)
  TSS<-sum((one_temp$srp.loss.rate.ug.g.hr-mean(one_temp$srp.loss.rate.ug.g.hr))^2)
  
  mm_paramaters$km[i]<-coef(model.nls)[2]
  mm_paramaters$vmax[i]<-coef(model.nls)[1]
  mm_paramaters$vmax.se[i]<-a$coefficients[1,2]
  mm_paramaters$km.se[i]<-a$coefficients[2,2]
  mm_paramaters$km.p[i]<-a$coefficients[2,4]
  mm_paramaters$vmax.p[i]<-a$coefficients[1,4]
  mm_paramaters$rsq[i]<-1-RSS/TSS

  
}

write.csv(mm_paramaters,file="./lab2/mm_fits.csv")

nls(srp.loss.rate.ug.g.hr ~ Vm * actual.initial.conc.ug.l/(K+actual.initial.conc.ug.l), data = one_temp, 
    start = list(K = 50, Vm = 10))



temperature.plot<-ggplot(d2,aes(x=actual.initial.conc.ug.l,y=srp.loss.rate.ug.g.hr,color=as.factor(temp.treatment)))+
    theme_classic()+geom_point(size=2)+
  geom_line(data=mm_fits,aes(x=concentraiton,y=model.fit,color=temperature),size=2)+
  #geom_smooth(se=FALSE)+
  theme(text = element_text(size=20),legend.position = c(0.2, 0.75))+
  scale_colour_brewer(palette="Spectral",direction=-1,name = "Temperature")+
  xlab(expression("Initial "*PO[4]~mu*g~L^-1))+
  ylab(expression("SRP uptake rate ("*mu*g~"hr"^-1~"g"^-1*")"))
temperature.plot


mm_paramaters$boltz.temp<-1/((as.numeric(mm_paramaters$temperature)+273)*Boltz)
mm_paramaters$centered.boltz.temp<-mm_paramaters$boltz.temp-1/((273+12)*Boltz)

km_temp_sens<-ggplot(mm_paramaters,aes(x=as.numeric(temperature),y=km))+
  scale_y_log10()+
  theme_classic()+geom_point(size=3)+
  #geom_smooth(method="lm",se=FALSE)+
  theme(text = element_text(size=20),legend.position = c(0.2, 0.8))+
  #xlab(expression("Temperature ("^o*"C)"))+
  xlab("")+
  scale_x_continuous(sec.axis=sec_axis(~1/(Boltz*(273+12))-1/((.+273)*Boltz),
                                       name=expression("Standardized temperature "*"(1/["*italic("kT")[12]-italic(kT) *"])")))+
  ylab(expression("Half saturation constant ("*italic("K"[m])*")"))+
  annotate("text",18,130,label=expression(italic(Ea)*" = NA "),size=8)+
  annotate("text",18,110,label=expression(italic(R^2)*" = 0.00"),size=8)
km_temp_sens

summary(lm(log(km)~boltz.temp,data=mm_paramaters))

umax_temp_sens<-ggplot(mm_paramaters,aes(x=as.numeric(temperature),y=vmax))+
  scale_y_log10()+
  theme_classic()+
  geom_smooth(method="lm",se=FALSE,size=2,color="grey60")+geom_point(size=3)+
  theme(text = element_text(size=20),legend.position = c(0.2, 0.8))+
  xlab(expression("Temperature ("^o*"C)"))+
  scale_x_continuous(sec.axis=sec_axis(~1/(Boltz*(273+12))-1/((.+273)*Boltz),
                                       name=NULL,labels=NULL))+
  ylab(expression("Maximum uptake rate ("*italic("U"[max])*")"))+
  annotate("text",18,12,label=expression(italic(Ea)*" = 0.54 "*italic(e)*"V"),size=8)+
  annotate("text",18,15,label=expression(italic(R^2)*" = 0.69"),size=8)
umax_temp_sens

km.model<-lm(log(km)~boltz.temp,data=mm_paramaters)
anova(km.model)
summary(km.model)

umax.model<-lm(log(vmax)~boltz.temp,data=mm_paramaters)
anova(umax.model)
summary(umax.model)


tiff(filename="./figures/mm_parms_temp_sens.tiff",units="in",res=800,width=8,height=16,compression="lzw")
plot_grid(temperature.plot,km_temp_sens,umax_temp_sens,labels=c("(a)","(b)","(c)"),ncol=1,label_x=0.9,label_y=0.95,label_size = 20,align="v")
dev.off()

jpeg(filename="./figures/cb_check.tiff",units="in",res=50,width=8,height=8)
temperature.plot
dev.off()



tiff(filename="./figures/mm_temps.tiff",units="in",res=800,width=10,height=6,compression="lzw")
temperature.plot
dev.off()

tiff(filename="./figures/conc_uptake.tiff",units="in",res=800,width=10,height=6,compression="lzw")
opposite.plot
dev.off()

d2$centered.temperature<-d2$temp.kt-1/((273+12)*Boltz)

slopes<-data.frame()

for (i in 1:nrow(new.concentrations)){
  one.concentration.subset<-d2[d2$actual.initial.conc.ug.l==new.concentrations$actual.initial.conc.ug.l[i],]  
  model<-lm(log.rate~centered.temperature,data=one.concentration.subset)
  model.sum<-summary(model)
  
  output<-data.frame(concentration=new.concentrations$actual.initial.conc.ug.l[i],activation.energy=-model.sum$coefficients[2,1],
                     activation.energy.se=model.sum$coefficients[2,2],
                     rsq=model.sum$r.squared,p.value=model.sum$coefficients[2,4],
                     intercept=model.sum$coefficients[1,1],intercept.se=model.sum$coefficients[1,2])
  
  slopes<-rbind(slopes,output)
  
  plot<-ggplot(one.concentration.subset,aes(x=temp.kt,y=log.rate))+geom_point()+geom_smooth(method="lm")
  print(plot)
  
}

write.csv(slopes,file="./lab2/model_fits_slopes.csv")

library(MuMIn)
mm.ea.model<-nls(activation.energy ~ Vm * concentration/(K+concentration), data = slopes, 
    start = list(K = 50, Vm = 10))
conc_ea_model<-lm(activation.energy~concentration,data=slopes)
null<-lm(activation.energy~1,data=slopes)
AICc(mm.ea.model,conc_ea_model,null)

best.fit<-data.frame(conc=seq(11,217,by=1))
best.fit$fit<-coef(mm.ea.model)[2]*best.fit$conc/(coef(mm.ea.model)[1]+best.fit$conc)

RSS<-sum(residuals(mm.ea.model)^2)
TSS<-sum((slopes$activation.energy-mean(slopes$activation.energy))^2)
rsq<-1-RSS/TSS

summary(mm.ea.model)

r.squaredGLMM(conc_ea_model)

slopes.plot<-ggplot(slopes,aes(x=concentration,y=activation.energy))+theme_classic()+geom_point(size=3)+
  #geom_smooth(method="lm",se=FALSE,color="grey60",size=2)+
  geom_line(data=best.fit,aes(x=conc,y=fit),color="grey60",size=2)+
  theme(text = element_text(size=20))+xlab(expression("Initial "*PO[4]~mu*g~L^-1))+ylab("Activation Energy (eV)")+
  geom_errorbar(aes(ymin=activation.energy-activation.energy.se,ymax=activation.energy+activation.energy.se),width=0.001)+
  annotate("text",180,-0.27,label=expression(italic(R)^2*" = 0.75"),size=8)+
  annotate("text",180,-0.2,label=expression(italic(k[s])*" = 61.4"~mu*g~L^-1),size=8)+
  annotate("text",180,-0.13,label=expression(italic(U[max])*" = 0.56"~italic("eV")),size=8)
slopes.plot

tiff(filename="./figures/concentration_uptake_ea_presentation.tiff",units="in",res=800,width=12,height=8,compression="lzw")
slopes.plot
dev.off()

tiff(filename="./figures/concentration_uptake_ea_presentation1.tiff",units="in",res=800,width=12,height=8,compression="lzw")
opposite.plot
dev.off()



conc_ea_model<-lm(activation.energy~concentration,data=slopes)
anova(conc_ea_model)
summary(conc_ea_model)

tiff(filename="./figures/concentraiton_uptake_ea.tiff",units="in",res=800,width=8,height=14,compression="lzw")
plot_grid(opposite.plot,slopes.plot,labels=c("(a)","(b)"),cols=1,label_x=0.88,label_y=0.97,label_size = 20,align = "v")
dev.off()

tiff(filename="./figures/mm_plot.tiff",units="in",res=800,width=8,height=8,compression="lzw")
temperature.plot
dev.off()

tiff(filename="./figures/ea_by_conc.tiff",units="in",res=800,width=10,height=6,compression="lzw")
slopes.plot
dev.off()

tiff(filename="./figures/temps_and_conc.tiff",units="in",res=800,width=8,height=10,compression="lzw")
plot_grid(opposite.plot,slopes.plot,labels="AUTO",cols=1,label_x=0.9,label_y=0.96,label_size = 20)
dev.off()

model1<-lm(log.rate~temp.kt*actual.initial.conc.ug.l,data=d2)
anova(model1)
summary(model1)


model2<-lm(log.rate~temp.kt*as.factor(actual.initial.conc.ug.l),data=d2)
anova(model2)
summary(model2)

###################################################

############################ MM models of the data

