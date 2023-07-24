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
          #calculate log odds effect size for normal distribution
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            toLogOdds(lmerData[4]+lmerData[5])-toLogOdds(lmerData[4])    
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            toLogOdds(lmerData[4]+lmerData[6])-toLogOdds(lmerData[4])       
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            toLogOdds(lmerData[4]+lmerData[7])-toLogOdds(lmerData[4]) 
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

#accumulate random datasets to show effects of distribution
#Ns - vector of total number of measurements (number of participants * numbers of trials * levels of within factor)
#numberOfTrialsVector - vector of numbers of trials 
#reps - number of simulations
randDatasets=function(Ns,numberOfTrialsVector=c(20),intercepts=c(0),reps=1000,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0){
  numTrials=length(numberOfTrialsVector)
  #generate data frame for saving data
  dataOfSims=data.frame(rep=rep(1:reps,length(Ns)*numTrials,each=4), #4(effects) values saved per repetition,
                        intercept=rep(intercepts,each=4*reps,length(Ns)*length(numberOfTrialsVector)),
                        numberOfTrials=rep(numberOfTrialsVector,each=4*reps*length(intercepts),length(Ns)),
                        N=rep(Ns,each=4*reps*numTrials*length(intercepts)),
                        effects=rep(c("none","factor1","factor2","factor1:factor2"),length(Ns)*numTrials*length(intercepts)*reps),
                        logOdds=rep(0,length(Ns)*4*numTrials*length(intercepts)*reps),
                        prob=rep(0,length(Ns)*4*numTrials*length(intercepts)*reps),
                        correctResponses=rep(0,length(Ns)*4*numTrials*length(intercepts)*reps),
                        correctResponsesGuessing=rep(F,length(Ns)*4*numTrials*length(intercepts)*reps))
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
          #save to data frame
          #log odds
          dataOfSims$logOdds[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="none")]=
            mean(testdata$logOdds[which(testdata$factor1==-0.5 & testdata$factor2==-0.5)])
          dataOfSims$logOdds[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            mean(testdata$logOdds[which(testdata$factor1==0.5 & testdata$factor2==-0.5)])
          dataOfSims$logOdds[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            mean(testdata$logOdds[which(testdata$factor1==-0.5 & testdata$factor2==0.5)])
          dataOfSims$logOdds[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            mean(testdata$logOdds[which(testdata$factor1==0.5 & testdata$factor2==0.5)])
          #prob
          dataOfSims$prob[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="none")]=
            mean(testdata$prob[which(testdata$factor1==-0.5 & testdata$factor2==-0.5)])
          dataOfSims$prob[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            mean(testdata$prob[which(testdata$factor1==0.5 & testdata$factor2==-0.5)])
          dataOfSims$prob[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            mean(testdata$prob[which(testdata$factor1==-0.5 & testdata$factor2==0.5)])
          dataOfSims$prob[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            mean(testdata$prob[which(testdata$factor1==0.5 & testdata$factor2==0.5)])
          #correctResponses
          dataOfSims$correctResponses[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="none")]=
            mean(testdata$correctResponses[which(testdata$factor1==-0.5 & testdata$factor2==-0.5)])
          dataOfSims$correctResponses[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            mean(testdata$correctResponses[which(testdata$factor1==0.5 & testdata$factor2==-0.5)])
          dataOfSims$correctResponses[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            mean(testdata$correctResponses[which(testdata$factor1==-0.5 & testdata$factor2==0.5)])
          dataOfSims$correctResponses[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            mean(testdata$correctResponses[which(testdata$factor1==0.5 & testdata$factor2==0.5)])
          #correctResponsesGuessing
          dataOfSims$correctResponsesGuessing[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="none")]=
            mean(testdata$correctResponsesGuessing[which(testdata$factor1==-0.5 & testdata$factor2==-0.5)])
          dataOfSims$correctResponsesGuessing[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            mean(testdata$correctResponsesGuessing[which(testdata$factor1==0.5 & testdata$factor2==-0.5)])
          dataOfSims$correctResponsesGuessing[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            mean(testdata$correctResponsesGuessing[which(testdata$factor1==-0.5 & testdata$factor2==0.5)])
          dataOfSims$correctResponsesGuessing[which(dataOfSims$rep==i & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            mean(testdata$correctResponsesGuessing[which(testdata$factor1==0.5 & testdata$factor2==0.5)])
        }      
      }
    }
  }
  return(dataOfSims)
}


#function to generate plots
#maxEffectSize - cutoff for maximum effect size to plot
saveplots=function(dataOfSims,name,maxEffectSize=5){
  #remove singular fits
  dataOfSimsNoSingularFits=dataOfSims[which(dataOfSims$singularFit==0),]
  #plot p-value separated by methods by intercept
  ggplot(dataOfSimsNoSingularFits,aes(x=toAcc(intercept),y=pValue,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    facet_wrap(~effects) +
    labs(y="p value", x="intercept") +
    theme_classic()
  ggsave(paste("figs/",name,"pValue.png",sep=''))
  #plot effect size separated by  methods by intercept
  ggplot(dataOfSimsNoSingularFits,aes(x=toAcc(intercept),y=effectSize,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    facet_wrap(~effects) +
    labs(y="effect size", x="intercept") + coord_cartesian(ylim(-maxEffectSize,maxEffectSize))  + 
    theme_classic()
  ggsave(paste("figs/",name,"effectSize.png",sep=''))
  #power
  library(plyr)
  dataOfSimsSummarizedReps=ddply(dataOfSims,
                                 .(intercept,type,numberOfTrials,effects,N),
                                 summarize,
                                 reps=sum(singularFit==0),
                                 singularFits=sum(singularFit==1),
                                 powerPositive=sum(pValue<.05 & effectSize>0 & !is.na(effectSize) & singularFit==0)/reps,
                                 powerNegative=sum(pValue<.05 & effectSize<0 & !is.na(effectSize) & singularFit==0)/reps,
                                 nsEffect=sum(pValue>=.05 & singularFit==0)/reps)
  #power to detect effect in positive direction
  ggplot(dataOfSimsSummarizedReps,aes(x=toAcc(intercept),y=powerPositive,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    facet_wrap(~effects) +
    labs(y="power", x="intercept") +
    theme_classic()
  ggsave(paste("figs/",name,"powerPositive.png",sep=''))
  #power to detect effect in negative direction
  ggplot(dataOfSimsSummarizedReps,aes(x=toAcc(intercept),y=powerNegative,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    facet_wrap(~effects) +
    labs(y="power", x="intercept") +
    theme_classic()
  ggsave(paste("figs/",name,"powerNegative.png",sep=''))
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
# intercepts=toLogOdds(c(0.1,0.9))
# reps=10
# dataOfSims11110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0)

#generate random seed (this should be random enough)
#sample(0:100000,1)
#52766
set.seed(52766)
#different combinations of participants and measurements per participant (Ns need to be multiple of numberOfTrialsVector*2)
Ns=c(1600) #use 1600 due to Brysbaert & Stevens despite of arbitrary effect sizes here to emphasize this value as minimal number of measurements
numberOfTrialsVector=c(10,20,40)
intercepts=toLogOdds(c(0.01,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.99))
reps=100
dataOfSims11110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0)
#larger random intercept
dataOfSims11120=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=2,randomError=0)
#no interactions
dataOfSims11010=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=0,randomIntercept=1,randomError=0)
#negative interactions
dataOfSims11_110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=-1,randomIntercept=1,randomError=0)
#effects of smaller effect size
dataOfSims55010=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=0.5,withinEffectSize=0.5,interactionEffectSize=0,randomIntercept=1,randomError=0)
#effects of larger effect size
dataOfSims22010=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=2,withinEffectSize=2,interactionEffectSize=0,randomIntercept=1,randomError=0)

#increase N for larger errors
# Ns=c(2000)
# dataOfSims11111=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=1)

#generate some plots
saveplots(dataOfSims11110, "11110")
saveplots(dataOfSims11120, "11120")
saveplots(dataOfSims11010, "11010")
saveplots(dataOfSims11_110, "11-110")
saveplots(dataOfSims55010, "55010")
saveplots(dataOfSims22010, "22010")

#comparison
set.seed(527661)
#different combinations of participants and measurements per participant (Ns need to be multiple of numberOfTrialsVector*2)
Ns=c(1600) #use 1600 due to Brysbaert & Stevens despite of arbitrary effect sizes here to emphasize this value as minimal number of measurements
numberOfTrialsVector=c(10,20,40)
intercepts=toLogOdds(c(0.01,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.99))
reps=100
dataOfSims11110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0)
#larger random intercept
dataOfSims11120=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=2,randomError=0)
#no interactions
dataOfSims11010=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=0,randomIntercept=1,randomError=0)
#negative interactions
dataOfSims11_110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=-1,randomIntercept=1,randomError=0)
#effects of smaller effect size
dataOfSims55010=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=0.5,withinEffectSize=0.5,interactionEffectSize=0,randomIntercept=1,randomError=0)
#effects of larger effect size
dataOfSims22010=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=2,withinEffectSize=2,interactionEffectSize=0,randomIntercept=1,randomError=0)

#increase N for larger errors
# Ns=c(2000)
# dataOfSims11111=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=1)

#generate some plots
saveplots(dataOfSims11110, "111102")
saveplots(dataOfSims11120, "111202")
saveplots(dataOfSims11010, "110102")
saveplots(dataOfSims11_110, "11-1102")
saveplots(dataOfSims55010, "550102")
saveplots(dataOfSims22010, "220102")

#show average of randomly generated datasets
datasets11010=randDatasets(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=0,randomIntercept=1,randomError=0)
ggplot(datasets11010,aes(x=toAcc(intercept),y=logOdds,color=effects,shape=as.factor(numberOfTrials))) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
  labs(y="logOdds", x="intercept") +
  theme_classic()
ggplot(datasets11010,aes(x=toAcc(intercept),y=prob,color=effects,shape=as.factor(numberOfTrials))) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
  labs(y="prob", x="intercept") +
  theme_classic()
ggplot(datasets11010,aes(x=toAcc(intercept),y=correctResponses,color=effects,shape=as.factor(numberOfTrials))) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
  labs(y="correctResponses", x="intercept") +
  theme_classic()
ggplot(datasets11010,aes(x=toAcc(intercept),y=correctResponsesGuessing,color=effects,shape=as.factor(numberOfTrials))) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
  labs(y="correctResponsesGuessing", x="intercept") +
  theme_classic()

#specific plots
#for specific effect and effect size
effect="factor2"
numberOfTrials=40
dataOfSimsPowerSpecific=dataOfSimsSummarizedReps[which(dataOfSimsSummarizedReps$effects==effect & dataOfSimsSummarizedReps$numberOfTrials==numberOfTrials),]
ggplot(dataOfSimsPowerSpecific,aes(x=toAcc(intercept),y=powerPositive,color=type)) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
  labs(y="power", x="intercept") +
  theme_classic()
