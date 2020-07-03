compare.calib<-function(data.all,region,W,x) {
#region<-params$Region[1]
#W<-params$W[1]
data_real<-subset(data_all,Edo_code==region)
data_real<-data_real[order(data_real$Year,data_real$Month,data_real$Day),]
#data_real<-subset(data_real,Confirmed.Cases>0)
data_real$time<-c(0:(nrow(data_real)-1))
data_real$Hbedsx1000<-data_real$Hbedsx10000/10

library(mFilter)
#smooth variables to be calibrated
data_real$Hist.Confirmed.Cases<-data_real$Confirmed.Cases
data_real$Hist.C.Deaths<-data_real$C.Deaths
data_real$Confirmed.Cases<- hpfilter(data_real$Confirmed.Cases,freq=100)$trend
data_real$C.Deaths<- hpfilter(data_real$C.Deaths,freq=100)$trend
data_real$Hist.Infection.Rate.hp<-hpfilter(data_real$Hist.Infection.Rate,freq=100)$trend
data_real$Hist.Death.Rate.hp<- hpfilter(data_real$Hist.Death.Rate,freq=100)$trend
#data_real$Hist.Infection.Rate.hp<- c(NA,hpfilter(data_real$Hist.Infection.Rate[2:length(data_real$Hist.Infection.Rate)],freq=100)$trend)
#data_real$Hist.Death.Rate.hp<- c(NA,hpfilter(data_real$Hist.Death.Rate[2:length(data_real$Hist.Death.Rate)],freq=100)$trend)
data_real$Mov_Index.hp<- hpfilter(data_real$Mov_Index,freq=100)$trend

#load mov_Data
 Mov.Data <<- approxfun( x =  data_real$time, #times
                    y = data_real$Mov_Index.hp, # Actual data points
                    method = "linear",
                    rule = 2)

#paramaters vector
#  x<-params[1,c('Infectivity.e','mortality.rate.base.e','average.delay.time.e','population.infected.with.COVID.e','overburden.impact.e','Average.Duration.Of.Infectivity.e','hospitalization.rate.e','Contact.Frequency.e')]
#  x<-as.numeric(x)
#  covid_UMLE(x)
#
times <- seq(0 , #inicial time #days
             max(data_real$time),#+30, #end time #days
             0.2 ); #time step #days

intg.method<-c("rk4")

#

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


library(deSolve)
out <- ode(y = InitialConditions,
           times = times,
           func = covid.epidemic,
           parms = parameters,
           method =intg.method )

out<-data.frame(out)
#subset for calibration comparison
out<-subset(out,time%in%c(0:max(times)))


#+++++++++++++++++++++++++++++++++++++
calib.times<-subset(data_real,Confirmed.Cases>0)$time

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




#+++++++++++++++++++++++++++++++++++++++



#rbind with rest of data
#append dates

dates<-as.Date(paste(2020,
c(
 paste("03",as.character(c(1:31)),sep="-"),
 paste("04",as.character(c(1:30)),sep="-"),
 paste("05",as.character(c(1:31)),sep="-"),
 paste("06",as.character(c(1:30)),sep="-"),
 paste("07",as.character(c(1:31)),sep="-"),
 paste("08",as.character(c(1:31)),sep="-"),
 paste("09",as.character(c(1:30)),sep="-")
),sep="-"))


#cases_simulated<-out$perceived.Infection.Rate
#deaths
#deaths_real<-data_real$Hist.Death.Rate.hp/1e6
#deaths_simulated<-out$Death.Rate


#
data_simulated<-data.frame(
                            Date_new=subset(dates,dates>=min(data_real$Date_new))[1:length(out$time)],
                            Confirmed.Cases=out$perceived.Confirmed.Cases*1e6,
                            C.Deaths=out$perceived.deaths*1e6,
                            Infection.Rate=out$perceived.Infection.Rate*1e6,
                            Death.Rate=out$perceived.Death.Rate*1e6,
                            R0=out$R0,
                            type='simulation'
                            )

#
data_simulated_real<-data.frame(
                            Date_new=subset(dates,dates>=min(data_real$Date_new))[1:length(out$time)],
                            Confirmed.Cases=out$Confirmed.Cases*1e6,
                            C.Deaths=out$deaths*1e6,
                            Infection.Rate=out$perceived.Infection.Rate*1e6,
                            Death.Rate=out$perceived.Death.Rate*1e6,
                            R0=out$R0,
                            type='simulation real'
                            )

#
 data_hist<-data_real[,c("Date_new","Hist.Confirmed.Cases","Hist.C.Deaths","Hist.Infection.Rate","Hist.Death.Rate")]
 #data_hist<-data_real[,c("Date_new","Confirmed.Cases","C.Deaths","Hist.Infection.Rate.hp","Hist.Death.Rate.hp")]
 colnames(data_hist)<-c("Date_new","Confirmed.Cases","C.Deaths","Infection.Rate","Death.Rate")
 data_hist$R0<-0
 data_hist$type<-'real'
 data<-rbind(data_simulated,data_hist,data_simulated_real)
 data$Region<-region
 data$W<-W
 return(data)
}

data_calib<-apply(params,1,function(x){compare.calib(
                                                     data_all,
                                                     as.numeric(x['Region']),
                                                     as.numeric(x['W']),
                                                      c(
                                                        as.numeric(x['Infectivity.e']),
                                                        as.numeric(x['mortality.rate.base.e']),
                                                        as.numeric(x['average.delay.time.e']),
                                                        #as.numeric(x['Social.distancing.policy.e']),
                                                        as.numeric(x['population.infected.with.COVID.e']),
                                                        as.numeric(x['overburden.impact.e']),
                                                        #as.numeric(x['Social.distancing.trigger.s1.e']),
                                                        as.numeric(x['Average.Duration.Of.Infectivity.e']),
                                                        as.numeric(x['hospitalization.rate.e']),
                                                        as.numeric(x['Contact.Frequency.e'])
                                                       )
                                                      )
                                        }
                    )
