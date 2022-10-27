library(ggplot2)
library(cowplot)

####### function rate of for carbon use

c_rate<-function(temp,Ea_02,reference_rate_O2){




rq<-0.85 ## respiratory quotent 

Boltz<-8.6173303*10**-5 ## boltzmann constant

centered.temp<-1/((temp+273)*Boltz)-1/(Boltz*(273+12))

oxygen_use_rate<-reference_rate_O2*exp(Ea_O2*centered.temp)/0.45*24*14 ## oxygen use in mg per g of leaf 0.45 converts
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

  initial_c<-700 * 0.45 ### initial mass of leaves in the stream
  


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


############## Example plot
sim_length<-3000

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


extrapolation<-ggplot(time_series[time_series$time<600,],aes(x=time,y=c.remaining.hot))+geom_line(color="firebrick",size=2)+theme_classic()+
  xlab("Day of simulation")+
  geom_line(aes(x=time,y=c.remaining),color="steelblue",size=2)+
  geom_line(aes(x=time,y=p.uptake/10),color="steelblue",linetype="dashed",size=2)+
  geom_line(aes(x=time,y=p.uptake.hot/10),color="firebrick",linetype="dashed",size=2)+
  scale_y_continuous( name = expression("Carbon remaining ("*g~m^-2*")"),
                      sec.axis = sec_axis( trans=~.*0.01, 
                                           name= expression("Cumulative SRP uptake ("*g~m^-2*")")))+
  theme(text=element_text(size=20))

extrapolation

test<-time_series[time_series$p.uptake>time_series$p.uptake.hot,]

tiff(filename="./figures/extrapolating_results_faster.tiff",units="in",res=800,width=9,height=6,compression="lzw")
extrapolation
dev.off()

