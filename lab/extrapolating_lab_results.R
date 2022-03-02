

####### funciton rate of for carbon use

c_rate<-function(temp,Ea_02,reference_rate_O2){




rq<-1.2 ## respiratory quotent 

Boltz<-8.6173303*10**-5 ## boltzmann constant

centered.temp<-1/((temp+273)*Boltz)-1/(Boltz*(273+12))

oxygen_use_rate<-reference_rate_O2*exp(Ea_O2*centered.temp)/0.45*24 ## oxygen use in mg per g of leaf 0.45 converts
## rates measured in mg o2 per g leaf to mg 02 per g leaf carbon -- assuming leaves are about 45% carbon
## 24 converts hours to daily rates

## units are not mg 02 per g leaf C per day

carbon_use_rate<-oxygen_use_rate/32*12 /1000 *rq ### converts to grams of carbon lost per gram of leaf carbon per day
return(carbon_use_rate)
}

########### function for P uptake temperature sensitivity


Ea_p<-(-0.47735)
Ea_p<-(-0.2)

reference_rate_p<-exp(0.66205)

p_rate<-function(temp,Ea_p,reference_rate_p){


  
Boltz<-8.6173303*10**-5

centered.temp<-1/((temp+273)*Boltz)-1/(Boltz*(273+12))

p_uptake_rate<- reference_rate_p*exp(Ea_p*centered.temp)*24/1000*0.45 ## converts to milligrams, and the 0.45 converts dry mass to carbon
return(p_uptake_rate) ## mg of p taken up per gram of leaf per day
}


simulate_uptake_and_decomp<-function(temp1,temp2,Ea_p,Ea_O2,reference_rate_p,reference_rate_O2,nutrient_effect){

initial_c<-250 ### initial mass of leaves in the stream
  


sim_length<-10000

time_series<-data.frame(time=1:sim_length,c.remaining=rep(NA,time=sim_length),p.uptake=rep(NA,time=sim_length),
                        c.remaining.hot=rep(NA,time=sim_length),p.uptake.hot=rep(NA,time=sim_length))


time_series$c.remaining[1]<-initial_c
time_series$c.remaining.hot[1]<-initial_c
time_series$p.uptake[1]<-0
time_series$p.uptake.hot[1]<-0





for (i in 2:sim_length){
  time_series$c.remaining[i]<-time_series$c.remaining[i-1]-time_series$c.remaining[i-1]*c_rate(temp1,Ea_O2,reference_rate_O2)*nutrient_effect
  
  time_series$p.uptake[i]<-time_series$p.uptake[i-1]+time_series$c.remaining[i-1]*p_rate(temp1,Ea_p,reference_rate_p)
  
  
  time_series$c.remaining.hot[i]<-time_series$c.remaining.hot[i-1]-time_series$c.remaining.hot[i-1]*c_rate(temp2,Ea_02,reference_rate_O2)*nutrient_effect
  
  time_series$p.uptake.hot[i]<-time_series$p.uptake.hot[i-1]+time_series$c.remaining.hot[i-1]*p_rate(temp2,Ea_p,reference_rate_p)
  


  
}

cold_sim_low_mass<-time_series[time_series$c.remaining<2.5,]
warm_sim_low_mass<-time_series[time_series$c.remaining.hot<2.5,]

#response_ratio<-min(warm_sim_low_mass$p.uptake.hot)/min(cold_sim_low_mass$p.uptake)

return(c(min(cold_sim_low_mass$p.uptake),min(warm_sim_low_mass$p.uptake.hot)))
}


set.seed(207)
nsims<-100

null_model<-simulate_uptake_and_decomp(12,16,-0.65,-0.65,exp(0.66205),exp(-3.23757),1)


output.row.names<-c("low.nutreint.cold","low.concentration.hot","high.concentraiton.cold",
                    "high.concentraiton.hot","high.concentraiton.cold.c.effect",
                    "high.concentraiton.hot.c.effect")

nutrient.status<-c("low","low","high","high","high","high")
temperature.status<-c("cold","hot","cold","hot","cold","hot")
c.effect<-c("no","no","no","no","yes","yes")

simulated_data<-data.frame()


for (i in 1:nsims){
  
  
  ###### Uptake happening at low nutrient concentration
  uptake_ea<-rnorm(1,mean=0.19492,sd=0.075*12**0.5)
  respiration_ea<-rnorm(1,mean=1.02,sd=0.06*25**0.5)
  
  low.concentration<-simulate_uptake_and_decomp(temp1=10,temp2=14,Ea_p=-uptake_ea,Ea_O2=-respiration_ea,
                                   reference_rate_p=exp(0.66205),reference_rate_O2=exp(-3.23757),
                                   nutrient_effect=1)
  

  
  ###### Uptake happening at high nutrient concentration
  uptake_ea<-rnorm(1,mean=0.52200,sd=0.09911*12**0.5)
  respiration_ea<-rnorm(1,mean=1.02,sd=0.06*25**0.5)
  
  high.concentration<-simulate_uptake_and_decomp(temp1=10,temp2=14,Ea_p=-uptake_ea,Ea_O2=-respiration_ea,
                                                reference_rate_p=exp(2.14915),reference_rate_O2=exp(-3.23757),
                                                nutrient_effect=1)
  
  
  high.concentration.c.effect<-simulate_uptake_and_decomp(temp1=10,temp2=14,Ea_p=-uptake_ea,Ea_O2=-respiration_ea,
                                                 reference_rate_p=exp(2.14915),reference_rate_O2=exp(-3.23757),
                                                 nutrient_effect=1.32)
  
  
  
  
  output<-data.frame(iteration=i,simulation=output.row.names,temperature=temperature.status,nutrients=nutrient.status,
                     c.effect=c.effect,cummulative.uptake=c(low.concentration,high.concentration,high.concentration.c.effect))
  simulated_data<-rbind(simulated_data,output)
  
}



simulated_data$order<-rep((2:7),times=nsims)

null_model<-simulate_uptake_and_decomp(10,14,-0.65,-0.65,exp(0.66205),exp(-3.23757),1)
null.plot.data<-data.frame(iteration=1:10,simulation="Null",temperature=c("cold","hot"),
                           nutrients="Null",c.effect="no",cummulative.uptake=null_model,
                           order=1)

simulated_data<-rbind(simulated_data,null.plot.data)

plot_data<-aggregate(cummulative.uptake~temperature+nutrients+c.effect,simulated_data,quantile,probs=c(0.025,0.5,0.975))
plot_data$groups<-as.factor(paste(plot_data$nutrients,plot_data$c.effect,sep="-"))

plot_data$groups<-factor(plot_data$groups,levels=c("Null-no","low-no","high-no","high-yes"))



scenerios<-ggplot(plot_data,aes(x=groups,y=cummulative.uptake[,2]/1000,fill=temperature))+
  geom_bar(position="dodge", stat="identity")+
  theme_classic()+ylab(expression("Cummulative SRP uptake (g m"^-2*")"))+xlab("")+
  theme(text =element_text(size=20),legend.position = c(0.18, 0.89))+
  scale_fill_manual(values=c("steelblue","firebrick"))+
  geom_errorbar(position=position_dodge(width=0.9),aes(ymin=cummulative.uptake[,1]/1000,ymax=cummulative.uptake[,3]/1000),width=0.0001)+
  scale_y_log10()+
  scale_x_discrete(labels=c("Null-no" = " Low P 
Ea=0.65", "low-no" = "Low P
measured Ea", "high-no" = "High P
measured Ea","high-yes"="High P
measured Ea
+ P effect on respiration"))

scenerios


############## Example plot

temp1<-10
temp2<-14
Ea_O2<-(-1.02)
reference_rate_O2<-exp(-3.23757)
Ea_p<-(-0.52200)
reference_rate_p<-exp(2.14915)




time_series<-data.frame(time=1:sim_length,c.remaining=rep(NA,time=sim_length),p.uptake=rep(NA,time=sim_length),
                        c.remaining.hot=rep(NA,time=sim_length),p.uptake.hot=rep(NA,time=sim_length))


time_series$c.remaining[1]<-initial_c
time_series$c.remaining.hot[1]<-initial_c
time_series$p.uptake[1]<-0
time_series$p.uptake.hot[1]<-0





for (i in 2:sim_length){
  time_series$c.remaining[i]<-time_series$c.remaining[i-1]-time_series$c.remaining[i-1]*c_rate(temp1,Ea_O2,reference_rate_O2)
  
  time_series$p.uptake[i]<-time_series$p.uptake[i-1]+time_series$c.remaining[i-1]*p_rate(temp1,Ea_p,reference_rate_p)
  
  
  time_series$c.remaining.hot[i]<-time_series$c.remaining.hot[i-1]-time_series$c.remaining.hot[i-1]*c_rate(temp2,Ea_02,reference_rate_O2)
  
  time_series$p.uptake.hot[i]<-time_series$p.uptake.hot[i-1]+time_series$c.remaining.hot[i-1]*p_rate(temp2,Ea_p,reference_rate_p)
  
  
  
  
}


extrapolation<-ggplot(time_series[time_series$time<3000,],aes(x=time,y=c.remaining.hot))+geom_line(color="firebrick",size=2)+theme_classic()+
  xlab("Day of simulation")+
  geom_line(aes(x=time,y=c.remaining),color="steelblue",size=2)+
  geom_line(aes(x=time,y=p.uptake/100),color="steelblue",linetype="dashed",size=2)+
  geom_line(aes(x=time,y=p.uptake.hot/100),color="firebrick",linetype="dashed",size=2)+
  scale_y_continuous( name = expression("Carbon remaining ("*g~m^-2*")"),
                      sec.axis = sec_axis( trans=~.*0.1, 
                                           name= expression("Cumulative SRP uptake ("*g~m^-2*")")))+
  theme(text=element_text(size=20))

extrapolation

tiff(filename="./figures/extrapolating_results_new.tiff",units="in",res=800,width=9,height=11,compression="lzw")
plot_grid(scenerios,extrapolation,labels="AUTO",ncol=1,label_x=0.9,label_y=0.92,label_size = 20)
dev.off()

