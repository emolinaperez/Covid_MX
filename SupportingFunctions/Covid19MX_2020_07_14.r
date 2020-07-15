covid.epidemic <- function(t, state, parameters) {
with(as.list(c(state,parameters)), {
  #Endogenous auxiliary variables
    #Social.distancing<-ifelse(t>Social.distancing.trigger.s1, Social.distancing.policy,1.0)
    Social.distancing<-Mov.Data(t)
    Contact.traicing<-ifelse(t>Contact.traicing.trigger.s1, Contact.traicing.policy,0)
    Probability.of.Contact.with.Infected.Person<-population.infected.with.COVID/Total.Population #[1]
    Susceptible.Contacts<-population.susceptible.to.COVID*Contact.Frequency*Social.distancing #[people/day]=[people]*[1/day]*[1]
    Contacts.between.Infected.and.Uninfected.People<-Susceptible.Contacts*Probability.of.Contact.with.Infected.Person #[people/day]=[people/day]*[1]
  #Flow variables
    hospitalization.rate<-min(c(1,hospitalization.rate))
    Infection.Rate<-Contacts.between.Infected.and.Uninfected.People*Infectivity #[people/day]=[people/day]*[1]
    population.infected.with.COVID.that.requires.hospitalization<-population.infected.with.COVID*hospitalization.rate #[people]=[people]*[1]
    health.system.load<-population.infected.with.COVID.that.requires.hospitalization/Health.system.carrying.capacity #[1]=[people]/[people]
    effect.of.overburden.on.mortality.rate<-1.0+atan((1/overburden.impact)*health.system.load) #[1]
    mortality.rate<-mortality.rate.base*effect.of.overburden.on.mortality.rate #[1]=[1]*[1]
    mortality.rate<-min(c(1,mortality.rate)) #[1]
    isolations.that.require.hospitalization<-Isolated.population*hospitalization.rate #[people]=[people]*[1]
    Recovered.isolations.without.hospitalization<-Isolated.population*(1-hospitalization.rate) #[people]=[people]*[1]
    Recovered.isolations.after.hospitalization<-isolations.that.require.hospitalization*(1-mortality.rate) #[people]=[people]*[1]
    Recovered.population.without.hospitalization<-population.infected.with.COVID*(1-hospitalization.rate) #[people]=[people]*[1]
    Recovered.population.after.hospitalization<-population.infected.with.COVID.that.requires.hospitalization*(1-mortality.rate) #[people]=[people]*[1]

    #Recovery.Rate.infected.wo.hospitalization<-Recovered.population.without.hospitalization/15
    #Recovery.Rate.infected.after.hospitalization<-Recovered.population.after.hospitalization/Average.Duration.Of.Infectivity
    #Recovery.Rate.infected<-  Recovery.Rate.infected.wo.hospitalization+Recovery.Rate.infected.after.hospitalization
    Recovery.Rate.infected<-(Recovered.population.without.hospitalization+Recovered.population.after.hospitalization)/Average.Duration.Of.Infectivity #[people]=[people]/[day]


    Recovery.Rate.isolated<-(Recovered.isolations.without.hospitalization+Recovered.isolations.after.hospitalization)/Average.Duration.Of.Infectivity #[people]=[people]/[day]
    Fatalities.infected.population<-population.infected.with.COVID.that.requires.hospitalization*mortality.rate #[people]=[people]*[1]
    Fatalities.isolated.population<-isolations.that.require.hospitalization*mortality.rate #[people]=[people]*[1]

    Death.Rate<-(Fatalities.infected.population+Fatalities.isolated.population)/Average.Duration.Of.Infectivity #[people/day]=[people]/[day]

    R0<-Contact.Frequency*Social.distancing*Infectivity*Average.Duration.Of.Infectivity*((population.susceptible.to.COVID)/Total.Population)
    Isolation.Rate<-min(c(perceived.population.infected.with.COVID*Contact.Frequency*Contact.traicing,population.susceptible.to.COVID)) # [people/day]=[people]*[1/day]*[1]
  #Information delays
    dperceived.Infection.Rate<-(Infection.Rate-perceived.Infection.Rate)/average.delay.time
    dperceived.population.infected.with.COVID<-(population.infected.with.COVID-perceived.population.infected.with.COVID)/average.delay.time
    dperceived.Death.Rate<-(Death.Rate-perceived.Death.Rate)/average.delay.timeD
    dperceived.Confirmed.Cases<-(Confirmed.Cases-perceived.Confirmed.Cases)/average.delay.time
    dperceived.deaths<-(deaths-perceived.deaths)/average.delay.timeD
  #State variables
    ddeaths<-Death.Rate
    dpopulation.susceptible.to.COVID<-(-1)*min(c(Infection.Rate+Isolation.Rate,population.susceptible.to.COVID)) # [people/day]
    dpopulation.infected.with.COVID<-max(c(Infection.Rate-Recovery.Rate.infected,-1*population.infected.with.COVID)) # [people/day]
    dpopulation.recovered.from.COVID<-Recovery.Rate.infected+Recovery.Rate.isolated # [people/day]
    dIsolated.population<-Isolation.Rate # [people/day]
    dConfirmed.Cases<-Infection.Rate
  list(c(dpopulation.susceptible.to.COVID,
         dpopulation.infected.with.COVID,
         dperceived.population.infected.with.COVID,
         dpopulation.recovered.from.COVID,
         ddeaths,
         dIsolated.population,
         dperceived.Infection.Rate,
         dperceived.Death.Rate,
         dConfirmed.Cases,
         dperceived.Confirmed.Cases,
         dperceived.deaths
       ),
        Recovery.Rate.infected=Recovery.Rate.infected,
        Recovered.population.without.hospitalization=Recovered.population.without.hospitalization,
        Recovered.population.after.hospitalization=Recovered.population.after.hospitalization,
        Isolation.Rate=Isolation.Rate,
        Death.Rate=Death.Rate,
        R0=R0
        )
})
}
