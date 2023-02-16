### Comparison of normal and binomial distributions for accuracy in forced choice tasks
#     Copyright (C) 2022  Leonardo Jost
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

#load libraries
library(lme4)
library(optimx)
library(ggplot2)

###functions

#generate dataset
#n - number of participants
#numberOfTrials - number of measurements per category
#intercept - mean intercept
generateData=function(n,numberOfTrials,intercept=0) {
  testdata=data.frame(ids=as.factor(rep(1:n,each=2)), #2 measurements per ID for within subjects factor
                      factor1=(rep(c(-0.5,0.5),each=n)), #between subjects
                      factor2=(rep(c(-0.5,0.5),n))) #within subjects
  #random intercept for id (normally distributed)
  testdata$randomIntercept=rep(rnorm(n),each=2)
  #means for log odds (true skill level anywhere between probability 0 and 1)
  effect1=1
  effect2=1
  effectInteraction=1
  testdata$logOdds=intercept+
    testdata$randomIntercept+
    effect1*testdata$factor1+
    effect2*testdata$factor2+
    effectInteraction*testdata$factor1*testdata$factor2
  #convert to probability for binomial distribution
  testdata$prob=exp(testdata$logOdds)/(1+exp(testdata$logOdds))
  #generate dependent variable (with random error)
  testdata$correctResponses=rbinom(2*n,numberOfTrials,prob=testdata$prob)
  #add chance level guessing to incorrectly answered trials
  testdata$correctResponses=testdata$correctResponses+rbinom(2*n,numberOfTrials-testdata$correctResponses,0.5)
  #correct for guessing but set negative values to 0
  testdata$correctResponsesGuessing=ifelse(testdata$correctResponses>numberOfTrials/2,testdata$correctResponses-numberOfTrials/2,0)
  #set weights and normalize data to 0..1
  testdata$weights=numberOfTrials
  testdata$correctResponses=testdata$correctResponses/testdata$weights
  testdata$weightsGuessing=numberOfTrials/2
  testdata$correctResponsesGuessing=testdata$correctResponsesGuessing/testdata$weightsGuessing
  return(testdata)
}


#get p values and effect sizes of fixed effects
getDataGlmer=function(testdata){
  #catch errors if all values 0 or 1 for high and low probabilities
  result=tryCatch(
    {
      #model containing all effects
      glmerModel=glmer(correctResponses~factor1*factor2+(1|ids),
                       family=binomial(),data=testdata,weights=weights,
                       control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
      #comparison models each not containing one fixed effect of interest
      glmerModel1=update(glmerModel,formula = ~. -factor1)
      glmerModel2=update(glmerModel,formula = ~. -factor2)
      glmerModel3=update(glmerModel,formula = ~. -factor1:factor2)
      return(c(anova(glmerModel,glmerModel1)$"Pr(>Chisq)"[2],
               anova(glmerModel,glmerModel2)$"Pr(>Chisq)"[2],
               anova(glmerModel,glmerModel3)$"Pr(>Chisq)"[2],
               fixef(glmerModel)))
    },
    error=function(cond) {
      message("Error with glmerModel:")
      message(cond)
      message("\n")
      return(rep(NA,7))
    }
  )
  return(result)
}
#same for guessing adjustment
getDataGlmerGuessing=function(testdata){
  result=tryCatch(
    {
      #testdata=testdata[which(testdata$correctResponsesGuessing>=0),] #exclude negative scores due to correcting for guessing
      glmerModel=glmer(correctResponsesGuessing~factor1*factor2+(1|ids),
                       family=binomial(),data=testdata,weights=weightsGuessing,
                       control = glmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
      glmerModel1=update(glmerModel,formula = ~. -factor1)
      glmerModel2=update(glmerModel,formula = ~. -factor2)
      glmerModel3=update(glmerModel,formula = ~. -factor1:factor2)
      return(c(anova(glmerModel,glmerModel1)$"Pr(>Chisq)"[2],
               anova(glmerModel,glmerModel2)$"Pr(>Chisq)"[2],
               anova(glmerModel,glmerModel3)$"Pr(>Chisq)"[2],
               fixef(glmerModel)))
    },
    error=function(cond) {
      message("Error with glmerModel with adjustment:")
      message(cond)
      message("\n")
      return(rep(NA,7))
    }
  )
  return(result)
}
#same for normal distribution
getDataLmer=function(testdata){
  result=tryCatch(
    {
      lmerModel=lmer(correctResponses~factor1*factor2+(1|ids),
                   data=testdata,REML=FALSE,
                   control = lmerControl(optimizer = "optimx",optCtrl = list(method = "bobyqa")))
      lmerModel1=update(lmerModel,formula = ~. -factor1)
      lmerModel2=update(lmerModel,formula = ~. -factor2)
      lmerModel3=update(lmerModel,formula = ~. -factor1:factor2)
      return(c(anova(lmerModel,lmerModel1)$"Pr(>Chisq)"[2],
               anova(lmerModel,lmerModel2)$"Pr(>Chisq)"[2],
               anova(lmerModel,lmerModel3)$"Pr(>Chisq)"[2],
               fixef(lmerModel)))
    },
    error=function(cond) {
      message("Error with lmerModel:")
      message(cond)
      message("\n")
      return(rep(NA,7))
    }
  )
  return(result)
}
#iterate over random data
#ns - vector of number of participants
#numberOfTrialsVector - vector of numbers of trials 
#reps - number of simulations
randSim=function(ns,numberOfTrialsVector=c(20),intercepts=c(0),reps=1000){
  numTrials=length(numberOfTrialsVector)
  #generate data frame for saving data
  dataOfSims=data.frame(rep=rep(1:reps,length(ns)*numTrials,each=9), #3(effects)*3(type) values saved per repetition,
                        intercept=rep(intercepts,each=9*reps,length(ns)*length(numberOfTrialsVector)),
                        numberOfTrials=rep(numberOfTrialsVector,each=9*reps*length(intercepts),length(ns)),
                        n=rep(ns,each=9*reps*numTrials*length(intercepts)),
                        effects=rep(c("factor1","factor2","factor1:factor2"),length(ns)*numTrials*length(intercepts)*3*reps),
                        pValue=rep(0,length(ns)*9*numTrials*length(intercepts)*reps),
                        effectSize=rep(0,length(ns)*9*numTrials*length(intercepts)*reps),
                        type=rep(c("normal","binomial","binomialGuessing"),length(ns)*numTrials*length(intercepts)*reps,each=3))
  #loop over combinations of number of Trials and participants and intercepts
  for(intercept in intercepts){
    #output to get some sense of progress
    print(paste("intercept:",intercept))
    for(numberOfTrials in numberOfTrialsVector){
      print(paste("numberOfTrials:",numberOfTrials))
      for(n in ns) {
        #output to get some sense of progress
        print(paste("n:",n))
        #repeat simulations and save data
        for(i in 1:reps){
          testdata=generateData(n,numberOfTrials,intercept)
          glmerData=getDataGlmer(testdata)
          glmerDataGuessing=getDataGlmerGuessing(testdata)
          lmerData=getDataLmer(testdata)
          #save to data frame
          #lmerData
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1")]=
            lmerData[1]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor2")]=
            lmerData[2]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1:factor2")]=
            lmerData[3]
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1")]=
            lmerData[5]      
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor2")]=
            lmerData[6]      
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1:factor2")]=
            lmerData[7]
          #glmerData
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1")]=
            glmerData[1]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor2")]=
            glmerData[2]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1:factor2")]=
            glmerData[3]
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1")]=
            glmerData[5]     
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor2")]=
            glmerData[6]      
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1:factor2")]=
            glmerData[7]
          #glmerDataGuessing
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1")]=
            glmerDataGuessing[1]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor2")]=
            glmerDataGuessing[2]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1:factor2")]=
            glmerDataGuessing[3]
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1")]=
            glmerDataGuessing[5]     
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor2")]=
            glmerDataGuessing[6]      
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects=="factor1:factor2")]=
            glmerDataGuessing[7]
        }      
      }
    }
  }
  return(dataOfSims)
}

###script
#test
 n=50
 numberOfTrials=20
 intercept=toLogOdds(0.1)
 testdata=generateData(n,numberOfTrials,intercept)
# glmerData=getDataGlmer(testdata)
# lmerData=getDataLmer(testdata)
#generate random seed (this should be random enough)
#sample(0:100000,1)
#52766
set.seed(52766)
#show different powers for different numbers of participants
ns=c(100)
numberOfTrialsVector=c(10,20,40)
intercepts=toLogOdds(c(0.01,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.99))
reps=100
dataOfSims=randSim(ns,numberOfTrialsVector,intercepts,reps)
#summarize data
for(effect in unique(dataOfSims$effects)){
  for(intercept in intercepts){
    for(numberOfTrials in numberOfTrialsVector){
      for(n in ns) {
        for(type in unique(dataOfSims$type)) {
          print(type)
          print(effect)
          print(n)
          print(numberOfTrials)
          print(intercept)
          print(mean(dataOfSims$effectSize[which(dataOfSims$type==type & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects==effect)]))
          print(sd(dataOfSims$effectSize[which(dataOfSims$type==type & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects==effect)]))
          print(sum(dataOfSims$pValue[which(dataOfSims$type==type & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$n==n & dataOfSims$effects==effect)]<0.05)/reps)
        }
      }
    }
  }
}



#plot by n
ggplot(significantDataFrame,aes(x=n,y=propSignificant,color=effects)) +
  geom_point() + geom_line() + labs(y="simulated power", x="N") +
  scale_x_continuous(breaks = ns) +
  theme_classic()
ggsave(paste("figs/SimulPowerN2.png",sep=""))
#save data
write.table(significantDataFrame,file="simulPowerN2.csv",sep=";", row.names = FALSE)

#get power for simulated binomial distributions
set.seed(88735)
numberOfEachTrials=c(12,16,24,32,48,64)
significantDataFrame=randSim(ns=c(100),F,numberOfEachTrials,1000,0.73,0)
#plot by number of trials
ggplot(significantDataFrame,aes(x=numberOfEachTrial,y=propSignificant,color=effects)) +
  geom_point() + geom_line() + labs(y="simulated power", x="trials per test") +
  scale_x_continuous(breaks = numberOfEachTrials) +
  theme_classic()
ggsave(paste("figs/SimulPowerTrials2.png",sep=""))
#save data
write.table(significantDataFrame,file="simulPowerTrials2.csv",sep=";", row.names = FALSE)