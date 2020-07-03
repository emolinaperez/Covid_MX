covid_UMLE<-function(x,verbose=FALSE){

x<-as.numeric(x)

times <- seq(0 , #inicial time #days
             max(data_real$time), #end time #days
             0.2 ); #time step #days

intg.method<-c("rk4")


InitialConditions <- c(
                       population.susceptible.to.COVID =  as.numeric(unique(data_real$Pop))/1e6,
                       population.infected.with.COVID = (max(c(1/1e6,min(data_real$Confirmed.Cases)/1e6)))*round(as.numeric(x[4])),
                       perceived.population.infected.with.COVID =(max(c(1/1e6,min(data_real$Confirmed.Cases)/1e6))),
                       population.recovered.from.COVID = 0,
                       deaths=0,
                       Isolated.population=0,
                       perceived.Infection.Rate=(max(c(1/1e6,min(data_real$Hist.Infection.Rate.hp)/1e6))),
                       perceived.Death.Rate=0,
                       Confirmed.Cases=(max(c(1/1e6,min(data_real$Confirmed.Cases)/1e6)))*round(as.numeric(x[4])),
                       perceived.Confirmed.Cases=(max(c(1/1e6,min(data_real$Confirmed.Cases)/1e6))),
                       perceived.deaths=0
                     )

parameters<-c( Infectivity = 0.2*max(c(0,round(x[1],4))), # [1] dimmensionless
               mortality.rate.base  =.10*max(c(0,round(x[2],4))), #*as.numeric(x['mortality.rate']), # [1] dimmensionless
               Contact.Frequency = max(c(1.0,2.0*max(c(0,round(x[8],4))))), # people/person/day
               Total.Population = as.numeric(unique(data_real$Pop))/1e6, # million people,
               average.delay.time =  15*max(c(0,round(x[3],4))), #days
               #Social.distancing.trigger.s1 =  max(data_real$time)*min(c(1.0,max(c(0,round(x[7],4))))),
               Contact.traicing.trigger.s1 = 1000,#max(data_real$time)*min(c(1.0,max(c(0,round(x[10],4))))),
               #Social.distancing.policy = min(c(1.0,max(c(0,round(x[4],4))))),
               Contact.traicing.policy = 0,#min(c(1.0,max(c(0,round(x[11],4))))),
               Health.system.carrying.capacity = as.numeric(unique(data_real$Hbedsx1000))*(as.numeric(unique(data_real$Pop))/1e6/1000),
               Average.Duration.Of.Infectivity = 10*max(c(0,round(x[6],4))), # days #this parameter cannot change too much !!!!!!!!!!!!!!
               hospitalization.rate = 0.30*max(c(0,round(x[7],4))), #[1] dimmensionless
               overburden.impact=0.03*max(c(0,round(x[5],4)))
              )

out <- ode(y = InitialConditions,
           times = times,
           func = covid.epidemic,
           parms = parameters,
           method =intg.method )

out<-data.frame(out)
#subset for calibration comparison
out<-subset(out,time%in%c(0:max(times)))


#compare real data, versus simulated data

calib.times<-subset(data_real,Confirmed.Cases>W)$time

if (length(calib.times>0)) {
  calib.times<-calib.times } else {
  calib.times<-subset(data_real,Confirmed.Cases>200)$time
  }


#compare real data, versus simulated data (considering only rates)
#cases
#  cases_real<-data_real$Hist.Infection.Rate.hp[data_real$time%in%calib.times]
#  cases_simulated<-out$perceived.Infection.Rate[out$time%in%calib.times]*1e6
#deaths
#  deaths_real<-data_real$Hist.Death.Rate.hp[data_real$time%in%calib.times]
  deaths_simulated<-out$perceived.Death.Rate[out$time%in%calib.times]*1e6

#compare cumulative values
#
#cases
  cases_real<-data_real$Confirmed.Cases[data_real$time%in%calib.times]
  cases_simulated<-out$perceived.Confirmed.Cases[out$time%in%calib.times]*1e6
#deaths
  deaths_real<-data_real$C.Deaths[data_real$time%in%calib.times]
  deaths_simulated<-out$perceived.deaths[out$time%in%calib.times]*1e6


#set weights
#for cases
#  w_c<-exp(W*c(1:length(cases_real)))
#  w_c<-w_c/sum(w_c)
#for deaths
#  w_d<-exp(W*c(1:length(deaths_real)))
#  w_d<-w_d/sum(w_d)

#create weighted time series
#  cases_real<-cases_real*w_c
#  cases_simulated<-cases_simulated*w_c
#  deaths_real<-deaths_real*w_d
#  deaths_simulated<-deaths_simulated*w_d

#Use Theil's decomposition for minimizing model's error
#cases
 cases.diff<-cases_real-cases_simulated
 cases.diff<-ifelse(cases.diff>100,100,cases.diff)
 cases.diff<-ifelse(cases.diff<(100*-1),-100,cases.diff)
 MSE.1<-mean(cases.diff^2)
 sd.diff.1<-min(c(100,sd(cases_real)-sd(cases_simulated)))
 U_S.1<-((sd.diff.1)^2)/MSE.1
 m.diff.1<-min(c(100,mean(cases_real)-mean(cases_simulated)))
 U_M.1<-((m.diff.1)^2)/MSE.1

#deaths
 death.diff<-deaths_real-deaths_simulated
 death.diff<-ifelse(death.diff>100,100,death.diff)
 death.diff<-ifelse(death.diff<(100*-1),-100,death.diff)
 MSE.2<-mean(death.diff^2)
 sd.diff.2<-min(c(100,sd(deaths_real)-sd(deaths_simulated)))
 U_S.2<-((sd.diff.2)^2)/MSE.2
 m.diff.2<-min(c(100,mean(deaths_real)-mean(deaths_simulated)))
 U_M.2<-((m.diff.2)^2)/MSE.2

#Objective function
  # U<-0.5*U_M.1+0.5*U_M.2
   U<-0.25*U_M.1+0.25*U_M.2+0.25*U_S.1+0.25*U_S.2

 rm(out)
 rm(deaths_real)
 rm(deaths_simulated)
 return(U)
 }
