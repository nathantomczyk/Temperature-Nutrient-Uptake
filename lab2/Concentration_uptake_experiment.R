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



custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

opposite.plot<-ggplot(d2,aes(x=temp.kt,y=log.rate,color=as.factor(actual.initial.conc.ug.l)))+theme_classic()+geom_point(size=2)+
  geom_smooth(method="lm",se=FALSE)+
  theme(text = element_text(size=20),legend.position = c(0.8, 0.18))+
  scale_colour_brewer(palette="Spectral",direction=-1,name = expression(PO[4]~mu*g~L^-1))+
  xlab("Inverse Temperature 1/kT")+
  guides(color=guide_legend(ncol=2))+
  ylab(expression("ln(SRP uptake rate ("*mu*g~"hr"^-1~"g"^-1*"))"))
opposite.plot


temperature.plot<-ggplot(d2,aes(x=actual.initial.conc.ug.l,y=log.rate,color=as.factor(temp.treatment)))+
  theme_classic()+geom_point(size=2)+
  geom_smooth(se=FALSE)+
  theme(text = element_text(size=20),legend.position = c(0.8, 0.23))+
  scale_colour_brewer(palette="Spectral",direction=-1,name = "Temperature")+
  xlab(expression("Initial "*PO[4]~mu*g~L^-1))+
  ylab(expression("ln(SRP uptake rate ("*mu*g~"hr"^-1~"g"^-1*"))"))
temperature.plot


slopes<-data.frame()

for (i in 1:nrow(new.concentrations)){
  one.concentration.subset<-d2[d2$actual.initial.conc.ug.l==new.concentrations$actual.initial.conc.ug.l[i],]  
  model<-lm(log.rate~temp.kt,data=one.concentration.subset)
  model.sum<-summary(model)
  
  output<-data.frame(concentration=new.concentrations$actual.initial.conc.ug.l[i],activation.energy=-model.sum$coefficients[2,1],
                     activation.energy.se=model.sum$coefficients[2,2])
  
  slopes<-rbind(slopes,output)
  
  plot<-ggplot(one.concentration.subset,aes(x=temp.kt,y=log.rate))+geom_point()+geom_smooth(method="lm")
  print(plot)
  
}

slopes.plot<-ggplot(slopes,aes(x=concentration,y=activation.energy))+theme_classic()+geom_point(size=2)+
  theme(text = element_text(size=20))+xlab(expression("Initial "*PO[4]~mu*g~L^-1))+ylab("Activation Energy (eV)")+
  geom_errorbar(aes(ymin=activation.energy-activation.energy.se,ymax=activation.energy+activation.energy.se),width=0.001)
        
slopes.plot

tiff(filename="./figures/concentraiton_uptake_ea.tiff",units="in",res=800,width=8,height=14,compression="lzw")
plot_grid(opposite.plot,slopes.plot,labels="AUTO",cols=1,label_x=0.9,label_y=0.93,label_size = 20)
dev.off()

tiff(filename="./figures/mm_plot.tiff",units="in",res=800,width=8,height=8,compression="lzw")
temperature.plot
dev.off()

tiff(filename="./figures/temps_and_conc.tiff",units="in",res=800,width=8,height=16,compression="lzw")
plot_grid(opposite.plot,slopes.plot,temperature.plot,labels="AUTO",cols=1,label_x=0.9,label_y=0.96,label_size = 20)
dev.off()

model1<-lm(log.rate~temp.kt*actual.initial.conc.ug.l,data=d2)
anova(model1)
summary(model1)


model2<-lm(log.rate~temp.kt*as.factor(actual.initial.conc.ug.l),data=d2)
summary(model2)

