#=================================================================================================================
# PART 1: Set up directory structure
#=================================================================================================================
#Set root directory
#in pc
 root<-"C:~\\CovidModel_MX\\"
 root<-"C:\\Users\\L03054557\\OneDrive\\Edmundo-ITESM\\1.Articulos\\14. Covid Analysis\\CovidModel_MX\\"
#in cloud
 root<-"D:\\2. Papers\\2. CovidAnalysis\\CovidModel_MX\\"

#Source model
 dir.model<-paste0(root,"\\SupportingFunctions\\")
 model.version<-"Covid19MX_2020_07_03.r"
 source(paste(dir.model,model.version,sep=""))

#Source MLE function
 dir.mle<-paste0(root,"\\SupportingFunctions\\")
 mle.version<-"Covid_MLE_2020_07_03.r"
 source(paste(dir.mle,mle.version,sep=""))

#Data repositories
  dir.Indata<-paste0(root,"\\SupportingData\\")
  calibration.date<-"2020_07_04"
  dir.harness<-paste0(root,"params_mx_all_",calibration.date,"\\")
  dir.Outdata<-paste0(root,"\\OutData\\")

#Data files
 io.table<-"IO_table_0630.csv"
 mov.table<-"junio27_indicemov.csv"
 pop.table<-"edades_final.csv"

#===========================================================================================================
# PART TWO: Load historical data and estimate MLE parameters
#===========================================================================================================
#Load data
 data_all<-read.csv(paste0(dir.Indata,io.table,sep=""))
 data_all[,c("Day","Month","Year")]<-do.call("rbind",lapply(strsplit(as.character(data_all$Date),"/"),function(x) { c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3]))   }))
 data_all$Date_new<-paste(data_all$Year,data_all$Month,data_all$Day,sep="-")
 data_all$Date_new<-as.Date(data_all$Date_new)
 data_all<-data_all[order(data_all$Year,data_all$Month,data_all$Day),]

#read age cohors and adjust population numbers
 pop_cohorts<-read.csv(paste0(dir.Indata,pop.table,sep=""))

#merge
 dim(data_all)
 data_all<-merge(data_all,pop_cohorts[,c("Edo_code","Share_above20")],by="Edo_code")
 dim(data_all)
 data_all<-data_all[order(data_all$Year,data_all$Month,data_all$Day),]

#update popualtion
 data_all$Pop<-data_all$Pop*data_all$Share_above20

#now we begin the iteration
 region<-unique(data_all$Edo_code)

#Calculate actual infention and death rates
 rates<- apply(
                data.frame(region=region),1,function(x){ pivot<-subset(data_all,Edo_code==as.numeric(x['region']))[,c("Edo_code","Confirmed.Cases","C.Deaths","Date_new")];
                                                         colnames(pivot)<-c("Edo_code","Hist.Infection.Rate","Hist.Death.Rate","Date_new");
                                                         pivot$Hist.Infection.Rate<-c(NA,diff(pivot$Hist.Infection.Rate));
                                                         pivot$Hist.Death.Rate<-c(NA,diff(pivot$Hist.Death.Rate));
                                                         pivot
                                                         }
              )

  rates<-do.call("rbind",rates)

#merge back with historial data
  dim(data_all)
  data_all<-Reduce(function(...) { merge(..., all=TRUE) }, list(data_all, rates))
  dim(data_all)


#load movility data
  mov_data<-read.csv(paste0(dir.Indata,mov.table,sep=""))
  mov_data$movilidad<-rowMeans(mov_data[,c(
                                           "retail_and_recreation_percent_change_from_baseline",
                                           "grocery_and_pharmacy_percent_change_from_baseline",
                                           "transit_stations_percent_change_from_baseline",
                                           "workplaces_percent_change_from_baseline"
                                           )
                                        ]
                                )/100
  mov_data$Mov_Index<-mov_data$movilidad+1.0
  mov_data[,c("Day","Month","Year")]<-do.call("rbind",lapply(strsplit(as.character(mov_data$date),"/"),function(x) { c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3]))   }))
  mov_data$Date_new<-paste(mov_data$Year,mov_data$Month,mov_data$Day,sep="-")
  mov_data$Date_new<-as.Date(mov_data$Date_new)
  mov_data$Edo_code<-mov_data$ENTIDAD
  mov_data$entidadFed<-NULL
  mov_data$date<-NULL
  mov_data$movilidad<-NULL

#merge movility data and historic data
 dim(data_all)
  data_all<-Reduce(function(...) { merge(...) }, list(data_all, mov_data))
 dim(data_all)

#remove NAs,
  data_all<-subset(data_all,is.na(Hist.Infection.Rate)==FALSE)
  dim(data_all)

#re-order data
  data_all<-data_all[order(data_all$Year,data_all$Month,data_all$Day),]

#set time series thersholds for calibration
  Ws<-c(0.0,50,100,200,500,1000,2000)

#for (j in 1:length(Ws))
#{
 j<-1

 W<-Ws[j]
for (i in 1:length(region))
{
 #i<-9
 data_real<-subset(data_all,Edo_code==region[i])
 data_real<-data_real[order(data_real$Year,data_real$Month,data_real$Day),]
 data_real$time<-c(0:(nrow(data_real)-1))
 data_real$Hbedsx1000<-data_real$Hbedsx10000/10

library(mFilter)
#smooth variables to be calibrated
 data_real$Confirmed.Cases<- hpfilter(data_real$Confirmed.Cases,freq=100)$trend
 data_real$C.Deaths<- hpfilter(data_real$C.Deaths,freq=100)$trend
 data_real$Hist.Infection.Rate.hp<-hpfilter(data_real$Hist.Infection.Rate,freq=100)$trend
 data_real$Hist.Death.Rate.hp<- hpfilter(data_real$Hist.Death.Rate,freq=100)$trend
 #data_real$Hist.Infection.Rate.hp<- c(NA,hpfilter(data_real$Hist.Infection.Rate[2:length(data_real$Hist.Infection.Rate)],freq=100)$trend)
 #data_real$Hist.Death.Rate.hp<- c(NA,hpfilter(data_real$Hist.Death.Rate[2:length(data_real$Hist.Death.Rate)],freq=100)$trend)
 data_real$Mov_Index.hp<- hpfilter(data_real$Mov_Index,freq=100)$trend

#test1<-data_real[,c("Date_new","Hist.Infection.Rate","Hist.Death.Rate","Mov_Index")]
#colnames(test1)<-c("Date_new","Hist.Infection.Rate","Hist.Death.Rate","Mov_Index")
#test1$type<-"real"
#test2<-data_real[,c("Date_new","Hist.Infection.Rate.hp","Hist.Death.Rate.hp","Mov_Index.hp")]
#colnames(test2)<-c("Date_new","Hist.Infection.Rate","Hist.Death.Rate","Mov_Index")
#test2$type<-"hp"
#test<-rbind(test1,test2)
#library(ggplot2)
#ggplot(test,aes(x=Date_new,y=Hist.Infection.Rate,colour=type))+geom_line()
#ggplot(test,aes(x=Date_new,y=Hist.Death.Rate,colour=type))+geom_line()
#ggplot(test,aes(x=Date_new,y=Mov_Index,colour=type))+geom_line()

#load mov_Data,
 Mov.Data <<- approxfun( x =  data_real$time, #times
                     y = data_real$Mov_Index.hp, # Actual data points
                     method = "linear",
                     rule = 2)

library(deSolve)
library(snow)
#Specify numbers of cores available for calibration
  nCore<- 25
#Define cluster
 cl <- makeSOCKcluster(names = rep('localhost',nCore))
 global.elements<-list("ode","data_real","covid.epidemic","W","Mov.Data")
 clusterExport(cl,global.elements,envir=environment())

#Load libary genoud
library(rgenoud)
#set seed for genetic optimization
set.seed(55555)
#Execute the optimization
out<-genoud(covid_UMLE,max=FALSE,
     nvars=8,
     starting.values = c(
                         1.0,  #infectivity
                         1.0,  #mortality.rate.base
                         1.0,  #average.delay.time
                         10,   #population.infected.with.COVID
                         1.20, #overburden.impact
                         0.5,   #Average.Duration.Of.Infectivity
                         1.0,    #hospitalization.rate
                         1.0#,    #Contact.Frequency
                        # 0.10,   #Contact.traicing.trigger.s1
                        # 0.05   #Contact.traicing.policy
                         ),
     pop.size=10000,
     Domains=matrix(c(#inferior limits
                         0.01, #infectivity
                         0.01, #mortality.rate.base
                         0.01, #average.delay.time
                         1.0,  #population.infected.with.COVID
                         0.01, #overburden.impact
                         0.01,  #Average.Duration.Of.Infectivity
                         0.01,   #hospitalization.rate
                         0.01,   #Contact.Frequency
                         #0.01,   #Contact.traicing.trigger.s1
                         #0.00,   #Contact.traicing.policy
                      #superior limits
                          10, #infectivity
                          10, #mortality.rate.base
                          10, #average.delay.time
                          500,#population.infected.with.COVID
                          10, #overburden.impact
                          10, #Average.Duration.Of.Infectivity
                          10, #hospitalization.rate
                          10#, #Contact.Frequency
                          #1.0, #Contact.traicing.trigger.s1
                          #0.50 #Contact.traicing.policy
                           ),
                          ncol=2),
     cluster=cl,
     print.level=1,
     boundary.enforcement=2,
     solution.tolerance=0.005) #original 0.001
   stopCluster(cl)
 out
#save result
  x<-c(
       round(out$par[1],4),
       round(out$par[2],4),
       round(out$par[3],4),
       round(out$par[4]),
       round(out$par[5],4),
       round(out$par[6],4),
       round(out$par[7],4),
       round(out$par[8],4)#,
       #round(out$par[9],4)#,
       #round(out$par[10],4),
       #round(out$par[11],4)
       )
  calib.params<-data.frame(t(x))
  colnames(calib.params)<-c(
                             'Infectivity.e',
                              'mortality.rate.base.e',
                              'average.delay.time.e',
                              'population.infected.with.COVID.e',
                              'overburden.impact.e',
                              'Average.Duration.Of.Infectivity.e',
                              'hospitalization.rate.e',
                              #'Contact.traicing.trigger.s1.e',
                              #'Contact.traicing.policy.e',
                              'Contact.Frequency.e'
                            )
  calib.params$Region<-region[i]
  calib.params$W<-W
  calib.params$value<-out$value
 write.csv(calib.params,paste(dir.harness,"params_",as.character(region[i]),"_",as.character(W),".csv",sep=""),row.names=FALSE)
 }

#}

#
#===========================================================================================================
# PART THREE: Load estimated parameters and compare calibration against historical record
#===========================================================================================================
#
#load parameters
filenames <- list.files(dir.harness, pattern="*.csv", full.names=FALSE)
params<-lapply(filenames, function (x) {read.csv(paste(dir.harness,x,sep=""))})
params<-do.call(rbind,params)

#find the run with the lowest error value
 params<-params[order(params$value),]


# Load data
data_all<-read.csv(paste0(dir.Indata,"IO_table_0612.csv",sep=""))
data_all[,c("Day","Month","Year")]<-do.call("rbind",lapply(strsplit(as.character(data_all$Date),"/"),function(x) { c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3]))   }))
data_all$Date_new<-paste(data_all$Year,data_all$Month,data_all$Day,sep="-")
data_all$Date_new<-as.Date(data_all$Date_new)
data_all<-data_all[order(data_all$Year,data_all$Month,data_all$Day),]

#read age cohors and adjust population numbers
pop_cohorts<-read.csv(paste0(dir.Indata,"edades_final.csv",sep=""))

#merge
dim(data_all)
data_all<-merge(data_all,pop_cohorts[,c("Edo_code","Share_above20")],by="Edo_code")
dim(data_all)
data_all<-data_all[order(data_all$Year,data_all$Month,data_all$Day),]

#update popualtion
data_all$Pop<-data_all$Pop*data_all$Share_above20

#now we begin the iteration
region<-unique(data_all$Edo_code)

#Calculate actual infention and death rates
rates<- apply(
               data.frame(region=region),1,function(x){ pivot<-subset(data_all,Edo_code==as.numeric(x['region']))[,c("Edo_code","Confirmed.Cases","C.Deaths","Date_new")];
                                                        colnames(pivot)<-c("Edo_code","Hist.Infection.Rate","Hist.Death.Rate","Date_new");
                                                        pivot$Hist.Infection.Rate<-c(NA,diff(pivot$Hist.Infection.Rate));
                                                        pivot$Hist.Death.Rate<-c(NA,diff(pivot$Hist.Death.Rate));
                                                        pivot
                                                        }
             )

 rates<-do.call("rbind",rates)

#merge back with historial data
 dim(data_all)
 data_all<-Reduce(function(...) { merge(..., all=TRUE) }, list(data_all, rates))
 dim(data_all)


#load movility data
  mov_data<-read.csv(paste0(dir.Indata,mov.table,sep=""))
  mov_data$movilidad<-rowMeans(mov_data[,c(
                                           "retail_and_recreation_percent_change_from_baseline",
                                           "grocery_and_pharmacy_percent_change_from_baseline",
                                           "transit_stations_percent_change_from_baseline",
                                           "workplaces_percent_change_from_baseline"
                                           )
                                        ]
                                )/100
  mov_data$Mov_Index<-mov_data$movilidad+1.0
  mov_data[,c("Day","Month","Year")]<-do.call("rbind",lapply(strsplit(as.character(mov_data$date),"/"),function(x) { c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3]))   }))
  mov_data$Date_new<-paste(mov_data$Year,mov_data$Month,mov_data$Day,sep="-")
  mov_data$Date_new<-as.Date(mov_data$Date_new)
  mov_data$Edo_code<-mov_data$ENTIDAD
  mov_data$entidadFed<-NULL
  mov_data$date<-NULL
  mov_data$movilidad<-NULL


#merge movility data and historic data
dim(data_all)
 data_all<-Reduce(function(...) { merge(...) }, list(data_all, mov_data))
dim(data_all)

#
#remove NAs,
  data_all<-subset(data_all,is.na(Hist.Infection.Rate)==FALSE)
  dim(data_all)

#re order data
  data_all<-data_all[order(data_all$Year,data_all$Month,data_all$Day),]

#Source compare calibration function
  dir.compare<-paste0(root,"\\SupportingFunctions\\")
  compare.version<-"CompareCalib_2020_07_03.r"
  source(paste(dir.compare,compare.version,sep=""))

#run calibration 
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


#run compare calib
 data_calib<-do.call('rbind',data_calib)

#after inspection, select most appropriate calibrations and print both params table and calibration table
#keep succesful calibration runs
 params$index<-paste(params$Region,params$W,sep="_")
 keeps<-c(
            "1_200",
            "2_1000",
            "3_200",
            "4_500",
            "5_500",
            "6_50",
            "7_500",
            "8_500",
            "9_50",
            "10_50",
            "11_500",
            "12_50",
            "13_1000",
            "14_1000",
            "15_0",
            "16_1000",
            "17_1000",
            "18_50",
            "19_100",
            "20_500",
            "21_1000",
            "22_500",
            "23_500",
            "24_200",
            "25_0",
            "26_200",
            "27_500",
            "28_1000",
            "29_1000",
            "30_200",
            "31_500",
            "32_1000"
          )

params<-subset(params,index%in%keeps)

params.names<-c(
                "Infectivity.e",
                "mortality.rate.base.e",
                "average.delay.time.e",
                "population.infected.with.COVID.e",
                "overburden.impact.e",
                "Average.Duration.Of.Infectivity.e",
                "hospitalization.rate.e",
                "Contact.Frequency.e"
                )
params<-params[,c(params.names,"Region","W")]

write.csv(params,paste(dir.Outdata,"params_2020_06_20.csv",sep=""),row.names=FALSE)
write.csv(data_calib,paste(dir.Outdata,"data_calib.csv",sep=""),row.names=FALSE)

#===========================================================================================================
# PART FOUR: Produce baseline conditions runs
#===========================================================================================================

#read params
 params<-read.csv(paste0(dir.Outdata,"params_2020_06_20.csv",sep=""))

 data_all<-read.csv(paste0(dir.data,"IO_table_0612.csv",sep=""))
 data_all[,c("Day","Month","Year")]<-do.call("rbind",lapply(strsplit(as.character(data_all$Date),"/"),function(x) { c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3]))   }))
 data_all$Date_new<-paste(data_all$Year,data_all$Month,data_all$Day,sep="-")
 data_all$Date_new<-as.Date(data_all$Date_new)
 data_all<-data_all[order(data_all$Year,data_all$Month,data_all$Day),]

 #read age cohors and adjust population numbers
 pop_cohorts<-read.csv(paste0(dir.data,"edades_final.csv",sep=""))

 #merge
 dim(data_all)
 data_all<-merge(data_all,pop_cohorts[,c("Edo_code","Share_above20")],by="Edo_code")
 dim(data_all)
 data_all<-data_all[order(data_all$Year,data_all$Month,data_all$Day),]

 #update popualtion
 data_all$Pop<-data_all$Pop*data_all$Share_above20

 #now we begin the iteration
 region<-unique(data_all$Edo_code)

 #Calculate actual infection and death rates
 rates<- apply(
                data.frame(region=region),1,function(x){ pivot<-subset(data_all,Edo_code==as.numeric(x['region']))[,c("Edo_code","Confirmed.Cases","C.Deaths","Date_new")];
                                                         colnames(pivot)<-c("Edo_code","Hist.Infection.Rate","Hist.Death.Rate","Date_new");
                                                         pivot$Hist.Infection.Rate<-c(NA,diff(pivot$Hist.Infection.Rate));
                                                         pivot$Hist.Death.Rate<-c(NA,diff(pivot$Hist.Death.Rate));
                                                         pivot
                                                         }
              )

  rates<-do.call("rbind",rates)

 #merge back with historial data
  dim(data_all)
  data_all<-Reduce(function(...) { merge(..., all=TRUE) }, list(data_all, rates))
  dim(data_all)


 #load movility data
  mov_data<-read.csv(paste0(dir.data,"junio7_indicemov.csv",sep=""))
  mov_data$Mov_Index<-mov_data$movilidad+1.0
  mov_data[,c("Day","Month","Year")]<-do.call("rbind",lapply(strsplit(as.character(mov_data$date),"/"),function(x) { c(as.numeric(x[1]),as.numeric(x[2]),as.numeric(x[3]))   }))
  mov_data$Date_new<-paste(mov_data$Year,mov_data$Month,mov_data$Day,sep="-")
  mov_data$Date_new<-as.Date(mov_data$Date_new)
  mov_data$Edo_code<-mov_data$entidadFed
  mov_data$entidadFed<-NULL
  mov_data$date<-NULL
  mov_data$movilidad<-NULL

 #merge movility data and historic data
 dim(data_all)
  data_all<-Reduce(function(...) { merge(...) }, list(data_all, mov_data))
 dim(data_all)

 #
 #remove NAs,
   data_all<-subset(data_all,is.na(Hist.Infection.Rate)==FALSE)
   dim(data_all)

 #re order data
   data_all<-data_all[order(data_all$Year,data_all$Month,data_all$Day),]

#
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
#
data_calib<-do.call('rbind',data_calib)

#add edo names
edo.names<-read.csv(paste0(dir.data,"IO_table_0531.csv",sep=""))
edo.names<-unique(edo.names[,c("Edo_code","edo")])
#merge
 data_calib$Edo_code<-data_calib$Region
dim(data_calib)
 data_calib<-merge(data_calib,edo.names,by="Edo_code")
dim(data_calib)
write.csv(data_calib,paste(dir.Outdata,"data_calib.csv",sep=""),row.names=FALSE)
