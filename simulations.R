### Comparison of normal and binomial distributions for accuracy in forced choice tasks
#     Copyright (C) 2023  Leonardo Jost
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
source("functions/helpers.R")

#generate dataset
#n - number of participants
#numberOfTrials - number of measurements per category
#intercept - mean intercept
#betweenEffectSize - mean effect for between subjects factor (factor1)
#withinEffectSize - mean effect for within subjects factor (factor2)
#interactionEffectSize - mean effect for interaction of factors
#randomIntercept - multiplicator for normally distributed random intercept
#randomError - multiplicator for normally distributed additional random error (some error is already introduced by binomial distribution of answers)
generateData=function(n,numberOfTrials,intercept=0,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0) {
  testdata=data.frame(ids=as.factor(rep(1:n,each=2)), #2 measurements per ID for within subjects factor
                      factor1=(rep(c(-0.5,0.5),each=n)), #between subjects
                      factor2=(rep(c(-0.5,0.5),n))) #within subjects
  #random intercept for id (normally distributed)
  testdata$randomIntercept=rep(rnorm(n),each=2)*randomIntercept
  #means for log odds (true skill level anywhere between probability 0 and 1)
  effect1=betweenEffectSize
  effect2=withinEffectSize
  effectInteraction=interactionEffectSize
  testdata$logOdds=intercept+
    testdata$randomIntercept+
    effect1*testdata$factor1+
    effect2*testdata$factor2+
    effectInteraction*testdata$factor1*testdata$factor2+
    rnorm(2*n)*randomError
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
               fixef(glmerModel),
               isSingular(glmerModel)))
    },
    error=function(cond) {
      message("Error with glmerModel:")
      message(cond)
      message("\n")
      return(rep(NA,8))
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
               fixef(glmerModel),
               isSingular(glmerModel)))
    },
    error=function(cond) {
      message("Error with glmerModel with adjustment:")
      message(cond)
      message("\n")
      return(rep(NA,8))
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
               fixef(lmerModel),
               isSingular(lmerModel)))
    },
    error=function(cond) {
      message("Error with lmerModel:")
      message(cond)
      message("\n")
      return(rep(NA,8))
    }
  )
  return(result)
}

#iterate over random data
#Ns - vector of total number of measurements (number of participants * numbers of trials * levels of within factor)
#numberOfTrialsVector - vector of numbers of trials 
#reps - number of simulations
randSim=function(Ns,numberOfTrialsVector=c(20),intercepts=c(0),reps=1000,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0){
  numTrials=length(numberOfTrialsVector)
  #generate data frame for saving data
  dataOfSims=data.frame(rep=rep(1:reps,length(Ns)*numTrials,each=9), #3(effects)*3(type) values saved per repetition,
                        intercept=rep(intercepts,each=9*reps,length(Ns)*length(numberOfTrialsVector)),
                        numberOfTrials=rep(numberOfTrialsVector,each=9*reps*length(intercepts),length(Ns)),
                        N=rep(Ns,each=9*reps*numTrials*length(intercepts)),
                        effects=rep(c("factor1","factor2","factor1:factor2"),length(Ns)*numTrials*length(intercepts)*3*reps),
                        type=rep(c("normal","binomial","binomialGuessing"),length(Ns)*numTrials*length(intercepts)*reps,each=3),
                        pValue=rep(0,length(Ns)*9*numTrials*length(intercepts)*reps),
                        effectSize=rep(0,length(Ns)*9*numTrials*length(intercepts)*reps),
                        singularFit=rep(F,length(Ns)*9*numTrials*length(intercepts)*reps))
  #loop over combinations of number of Trials and participants and intercepts
  for(intercept in intercepts){
    #output to get some sense of progress
    print(paste("intercept:",intercept))
    for(numberOfTrials in numberOfTrialsVector){
      print(paste("numberOfTrials:",numberOfTrials))
      for(N in Ns) {
        #output to get some sense of progress
        print(paste("N:",N))
        #calculate number of participants
        n=N/numberOfTrials/2
        print(paste("n:",n))
        #repeat simulations and save data
        for(i in 1:reps){
          testdata=generateData(n,numberOfTrials,intercept,betweenEffectSize,withinEffectSize,interactionEffectSize,randomIntercept,randomError)
          glmerData=getDataGlmer(testdata)
          glmerDataGuessing=getDataGlmerGuessing(testdata)
          lmerData=getDataLmer(testdata)
          #save to data frame
          #lmerData
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            lmerData[1]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            lmerData[2]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            lmerData[3]
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            lmerData[5]      
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            lmerData[6]      
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            lmerData[7]
          dataOfSims$singularFit[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N)]=
            lmerData[8]
          #glmerData
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            glmerData[1]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            glmerData[2]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            glmerData[3]
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            glmerData[5]     
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            glmerData[6]      
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            glmerData[7]
          dataOfSims$singularFit[which(dataOfSims$rep==i & dataOfSims$type=="binomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N)]=
            glmerData[8]
          #glmerDataGuessing
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            glmerDataGuessing[1]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            glmerDataGuessing[2]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            glmerDataGuessing[3]
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            glmerDataGuessing[5]     
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            glmerDataGuessing[6]      
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            glmerDataGuessing[7]
          dataOfSims$singularFit[which(dataOfSims$rep==i & dataOfSims$type=="binomialGuessing" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N)]=
            glmerDataGuessing[8]
        }      
      }
    }
  }
  return(dataOfSims)
}

###script
# #test
# N=50
# numberOfTrials=20
# intercept=toLogOdds(0.01)
# betweenEffectSize=1
# withinEffectSize=1
# interactionEffectSize=1
# randomIntercept=1
# randomError=0
# testdata=generateData(N,numberOfTrials,intercept,betweenEffectSize,withinEffectSize,interactionEffectSize,randomIntercept,randomError)
# glmerData=getDataGlmer(testdata)
# glmerDataGuessing=getDataGlmerGuessing(testdata)
# lmerData=getDataLmer(testdata)
# 
# Ns=c(1200)
# numberOfTrialsVector=c(20)
# intercepts=toLogOdds(0.9)
# reps=10
# dataOfSims11110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0)

#generate random seed (this should be random enough)
#sample(0:100000,1)
#52766
set.seed(52766)
#different combinations of participants and measurements per participant (Ns need to be multiple of numberOfTrialsVector*2)
Ns=c(1200)
numberOfTrialsVector=c(10,20,40)
intercepts=toLogOdds(c(0.01,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.99))
reps=100
dataOfSims11110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0)
#larger random intercept
dataOfSims11120=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=2,randomError=0)
#negative interactions
dataOfSims11_110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=-1,randomIntercept=1,randomError=0)
#increase N for larger errors
Ns=c(2000)
dataOfSims1111=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=1)


dataOfSims=dataOfSims1110
#plot p-value separated by methods by intercept
ggplot(dataOfSims,aes(x=toAcc(intercept),y=pValue,color=type,shape=as.factor(numberOfTrials))) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
  facet_wrap(~effects) +
  labs(y="p value", x="intercept") +
  theme_classic()
#plot effect size separated by  methods by intercept
ggplot(dataOfSims,aes(x=toAcc(intercept),y=effectSize,color=type,shape=as.factor(numberOfTrials))) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
  facet_wrap(~effects) +
  labs(y="p value", x="intercept") +
  theme_classic()
#power
library(plyr)
dataOfSimsPower=ddply(dataOfSims,
                      .(intercept,type,numberOfTrials,effects,N),
                      summarize,
                      power=sum(pValue<.05)/reps)
ggplot(dataOfSimsPower,aes(x=toAcc(intercept),y=power,color=type,shape=as.factor(numberOfTrials))) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
  facet_wrap(~effects) +
  labs(y="power", x="intercept") +
  theme_classic()
#for specific effect and effect size
effect="factor2"
numberOfTrials=40
dataOfSimsPowerSpecific=dataOfSimsPower[which(dataOfSimsPower$effects==effect & dataOfSimsPower$numberOfTrials==numberOfTrials),]
ggplot(dataOfSimsPowerSpecific,aes(x=toAcc(intercept),y=power,color=type)) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
  labs(y="power", x="intercept") +
  theme_classic()
