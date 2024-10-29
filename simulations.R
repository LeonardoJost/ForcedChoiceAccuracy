### Comparison of normal and binomial distributions for accuracy in forced choice tasks
#     Copyright (C) 2024  Leonardo Jost
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
  testdata$correctResponsesGuessing=ifelse(testdata$correctResponses>numberOfTrials/2,testdata$correctResponses-numberOfTrials/2,0)*2
  #set weights and normalize data to 0..1
  testdata$weights=numberOfTrials
  testdata$correctResponses=testdata$correctResponses/testdata$weights
  #testdata$weightsGuessing=numberOfTrials/2 #same effect as multiplying scores by 2
  testdata$weightsGuessing=numberOfTrials
  testdata$correctResponsesGuessing=testdata$correctResponsesGuessing/testdata$weightsGuessing
  return(testdata)
}

#get p values and effect sizes of fixed effects
getDataGlmer=function(testdata){
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
           isSingular(glmerModel),
           confint(glmerModel,method="Wald")))
  
}
#same for guessing adjustment
getDataGlmerGuessing=function(testdata){
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
           isSingular(glmerModel),
           confint(glmerModel,method="Wald")))
}
#same for normal distribution
getDataLmer=function(testdata){
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
           isSingular(lmerModel),
           confint(lmerModel,method="Wald")))
}

#iterate over random data
#Ns - vector of total number of measurements (number of participants * numbers of trials * levels of within factor)
#numberOfTrialsVector - vector of numbers of trials 
#intercepts - vector of intercepts
#reps - number of simulations
#betweenEffectSize, withinEffectSize, interactionEffectSize, randomIntercept, randomError - values for data generation
randSim=function(Ns,numberOfTrialsVector=c(20),intercepts=c(0),reps=1000,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0){
  #set vectors/length for creation of dataset
  numTrials=length(numberOfTrialsVector)
  effectsVector=c("factor1","factor2","factor1:factor2")
  numberOfEffects=length(effectsVector)
  typesVector=c("normal","binomial","binomialGuessing","postHocBinomial")
  numberOfTypes=length(typesVector)
  #generate data frame for saving data
  dataOfSims=data.frame(rep=rep(1:reps,length(Ns)*numTrials,each=numberOfEffects*numberOfTypes), #3(effects)*3(type) values saved per repetition,
                        intercept=rep(intercepts,each=numberOfEffects*numberOfTypes*reps,length(Ns)*length(numberOfTrialsVector)),
                        numberOfTrials=rep(numberOfTrialsVector,each=numberOfEffects*numberOfTypes*reps*length(intercepts),length(Ns)),
                        N=rep(Ns,each=numberOfEffects*numberOfTypes*reps*numTrials*length(intercepts)),
                        effects=rep(effectsVector,length(Ns)*numTrials*length(intercepts)*numberOfTypes*reps),
                        type=rep(typesVector,length(Ns)*numTrials*length(intercepts)*reps,each=numberOfEffects),
                        pValue=rep(0,length(Ns)*numberOfEffects*numberOfTypes*numTrials*length(intercepts)*reps),
                        effectSize=rep(0,length(Ns)*numberOfEffects*numberOfTypes*numTrials*length(intercepts)*reps),
                        singularFit=rep(F,length(Ns)*numberOfEffects*numberOfTypes*numTrials*length(intercepts)*reps))
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
          #generate dataset
          testdata=generateData(n,numberOfTrials,intercept,betweenEffectSize,withinEffectSize,interactionEffectSize,randomIntercept,randomError)
          #perform analyses
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
            setMinMax(toLogOdds(lmerData[4]+0.5*lmerData[5])-toLogOdds(lmerData[4]-0.5*lmerData[5]))
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            setMinMax(toLogOdds(lmerData[4]+0.5*lmerData[6])-toLogOdds(lmerData[4]-0.5*lmerData[6]))      
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="normal" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            setMinMax(toLogOdds(lmerData[4]+0.5*lmerData[7])-toLogOdds(lmerData[4]-0.5*lmerData[7]))
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
          #post hoc binomial correction
          #calculate values
          interceptFixEf=glmerData[4]
          factor1FixEf=glmerData[5] 
          factor2FixEf=glmerData[6]
          interactionFixEf=glmerData[7]
          #interceptLowerCi=glmerData[10]
          #factor1LowerCi=glmerData[11]
          #factor2LowerCi=glmerData[12]
          interactionLowerCi=glmerData[13]
          #interceptUpperCi=glmerData[15]
          #factor1UpperCi=glmerData[16]
          #factor2UpperCi=glmerData[17]
          interactionUpperCi=glmerData[18]
          #corrected effect size for main effects, 0.5 due to scaling of original effects
          correctedFac1=correctChanceLevelLogOdds(interceptFixEf+0.5*factor1FixEf)-correctChanceLevelLogOdds(interceptFixEf-0.5*factor1FixEf)
          correctedFac2=correctChanceLevelLogOdds(interceptFixEf+0.5*factor2FixEf)-correctChanceLevelLogOdds(interceptFixEf-0.5*factor2FixEf)
          #interaction
          correctedFullInteraction=correctChanceLevelLogOdds(interceptFixEf+0.5*factor1FixEf+0.5*factor2FixEf+0.5*0.5*interactionFixEf)+
            correctChanceLevelLogOdds(interceptFixEf-0.5*factor1FixEf-0.5*factor2FixEf+0.5*0.5*interactionFixEf)-
            correctChanceLevelLogOdds(interceptFixEf+0.5*factor1FixEf-0.5*factor2FixEf-0.5*0.5*interactionFixEf)-
            correctChanceLevelLogOdds(interceptFixEf-0.5*factor1FixEf+0.5*factor2FixEf-0.5*0.5*interactionFixEf)
          #confidence interval
          correctedFullInteractionLowerCI=correctChanceLevelLogOdds(interceptFixEf+0.5*factor1FixEf+0.5*factor2FixEf+0.5*0.5*interactionLowerCi)+
            correctChanceLevelLogOdds(interceptFixEf-0.5*factor1FixEf-0.5*factor2FixEf+0.5*0.5*interactionLowerCi)-
            correctChanceLevelLogOdds(interceptFixEf+0.5*factor1FixEf-0.5*factor2FixEf-0.5*0.5*interactionLowerCi)-
            correctChanceLevelLogOdds(interceptFixEf-0.5*factor1FixEf+0.5*factor2FixEf-0.5*0.5*interactionLowerCi)
          correctedFullInteractionUpperCI=correctChanceLevelLogOdds(interceptFixEf+0.5*factor1FixEf+0.5*factor2FixEf+0.5*0.5*interactionUpperCi)+
            correctChanceLevelLogOdds(interceptFixEf-0.5*factor1FixEf-0.5*factor2FixEf+0.5*0.5*interactionUpperCi)-
            correctChanceLevelLogOdds(interceptFixEf+0.5*factor1FixEf-0.5*factor2FixEf-0.5*0.5*interactionUpperCi)-
            correctChanceLevelLogOdds(interceptFixEf-0.5*factor1FixEf+0.5*factor2FixEf-0.5*0.5*interactionUpperCi)
          #only calculate significance instead of approximated p-value
          #significanceInteraction=ifelse(sign(correctedFullInteractionLowerCI)==sign(correctedFullInteractionUpperCI),0,1)
          #account better for NaNs due to both multiple Inf values (if effect is <0 and upperCi<0 or effect >0 and lowerCI>0)
          #only considering upperCI<0 or lowerCI>0 could also be sufficient, but the effect size needs a sign to detect the direction anyway
          significanceInteraction=ifelse((sign(correctedFullInteraction)<0 && sign(correctedFullInteractionUpperCI<0)) ||
                                           (sign(correctedFullInteraction)>0 && sign(correctedFullInteractionLowerCI>0)),0,1)
          #save values
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="postHocBinomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            glmerData[1]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="postHocBinomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            glmerData[2]
          dataOfSims$pValue[which(dataOfSims$rep==i & dataOfSims$type=="postHocBinomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            significanceInteraction
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="postHocBinomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1")]=
            correctedFac1
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="postHocBinomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor2")]=
            correctedFac2
          dataOfSims$effectSize[which(dataOfSims$rep==i & dataOfSims$type=="postHocBinomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N & dataOfSims$effects=="factor1:factor2")]=
            correctedFullInteraction
          dataOfSims$singularFit[which(dataOfSims$rep==i & dataOfSims$type=="postHocBinomial" & dataOfSims$intercept==intercept & dataOfSims$numberOfTrials==numberOfTrials & dataOfSims$N==N)]=
            glmerData[8]
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
#dataOfSims - dataset to plot
#name - name for saving (string)
#maxEffectSize - cutoff for maximum effect size to plot
saveplots=function(dataOfSims,name,maxEffectSize=2){
  #rename type for better legends
  dataOfSims$type[which(dataOfSims$type=="postHocBinomial")]="Overall Correction"
  dataOfSims$type[which(dataOfSims$type=="binomialGuessing")]="Individual Correction"
  dataOfSims$type[which(dataOfSims$type=="binomial")]="Binomial"
  dataOfSims$type[which(dataOfSims$type=="normal")]="Normal"
  #remove singular fits
  dataOfSimsNoSingularFits=dataOfSims[which(dataOfSims$singularFit==0),]
  #set infinite effect sizes to +-5 for plotting
  dataOfSimsNoSingularFits$effectSize=setMinMax(dataOfSimsNoSingularFits$effectSize)
  #plot p-value separated by methods by intercept
  ggplot(dataOfSimsNoSingularFits,aes(x=toAcc(intercept),y=pValue,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    facet_wrap(~effects) +
    labs(y="p value", x="intercept", color = "type of analysis", shape = "number of trials") +
    theme_bw()+theme(legend.position = "none")
  ggsave(paste("figs/",name,"pValue.png",sep=''))
  #plot effect size separated by  methods by intercept
  ggplot(dataOfSimsNoSingularFits,aes(x=toAcc(intercept),y=effectSize,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    facet_wrap(~effects) +
    labs(y="effect size", x="intercept", color = "type of analysis", shape = "number of trials") + coord_cartesian(ylim=c(-maxEffectSize,maxEffectSize))  + 
    theme_bw()+theme(legend.position = "none")
  ggsave(paste("figs/",name,"effectSize.png",sep=''))
  #power
  library(plyr)
  dataOfSimsSummarizedReps=ddply(dataOfSims,
                                 .(intercept,type,numberOfTrials,effects,N),
                                 summarize,
                                 reps=sum(singularFit==0 & !is.na(pValue)),
                                 singularFits=sum(singularFit==1 | is.na(pValue)),
                                 powerPositive=sum(!is.na(pValue) & pValue<.05 & effectSize>0 & !is.na(effectSize) & singularFit==0)/reps,
                                 powerNegative=sum(!is.na(pValue) & pValue<.05 & effectSize<0 & !is.na(effectSize) & singularFit==0)/reps,
                                 nsEffect=sum(!is.na(pValue) & pValue>=.05 & singularFit==0)/reps)
  #power to detect effect in positive direction
  ggplot(dataOfSimsSummarizedReps,aes(x=toAcc(intercept),y=powerPositive,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    #stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    facet_wrap(~effects) +
    labs(y="power", x="intercept", color = "type of analysis", shape = "number of trials") +
    theme_bw()+theme(legend.position = "none")
  ggsave(paste("figs/",name,"powerPositive.png",sep=''))
  #power to detect effect in negative direction
  ggplot(dataOfSimsSummarizedReps,aes(x=toAcc(intercept),y=powerNegative,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    #stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    facet_wrap(~effects) +
    labs(y="power", x="intercept", color = "type of analysis", shape = "number of trials") +
    theme_bw()+theme(legend.position = "none")
  ggsave(paste("figs/",name,"powerNegative.png",sep=''))
  #number of singular fits
  ggplot(dataOfSimsSummarizedReps,aes(x=toAcc(intercept),y=singularFits,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    #stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    #facet_wrap(~effects) +
    labs(y="Singular Fits", x="Intercept", color = "Type of Analysis", shape = "Number of Trials") +
    theme_bw()+theme(legend.position = "none")
  #+theme(legend.position = c(0.5,0.7), legend.background = element_rect(fill="white",linewidth=0.5,linetype="solid",color="black"))
  ggsave(paste("figs/",name,"singularFits.png",sep=''),width=1000, height=1000,unit="px",dpi=200)
}

#function to generate plots specific to one effect
#dataOfSims - dataset to plot
#name - name for saving (string)
#minEffectSize, maxEffectSize - cutoff for minimum/maximum effect size to plot
#maxPowerPos, maxPowerNeg - cutoff for maximum values of power to plot for detection of positive and negative effects
savePlotsSpecific=function(dataOfSims,name,specificEffect="factor1:factor2",minEffectSize=-1,maxEffectSize=1,maxPowerPos=0.25,maxPowerNeg=0.25){
  #rename type for better legends
  dataOfSims$type[which(dataOfSims$type=="postHocBinomial")]="Overall Correction"
  dataOfSims$type[which(dataOfSims$type=="binomialGuessing")]="Individual Correction"
  dataOfSims$type[which(dataOfSims$type=="binomial")]="Binomial"
  dataOfSims$type[which(dataOfSims$type=="normal")]="Normal"
  #only specific effect
  dataOfSimsSpecificEffect=dataOfSims[which(dataOfSims$effects==specificEffect),]
  #remove singular fits
  dataOfSimsNoSingularFits=dataOfSimsSpecificEffect[which(dataOfSimsSpecificEffect$singularFit==0),]
  #set infinite effect sizes to +-5 for plotting
  dataOfSimsNoSingularFits$effectSize=setMinMax(dataOfSimsNoSingularFits$effectSize)
  #plot effect size separated by  methods by intercept
  ggplot(dataOfSimsNoSingularFits,aes(x=toAcc(intercept),y=effectSize,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    labs(y="Effect Size", x="Intercept", color = "Type of Analysis", shape = "Number of Trials") + coord_cartesian(ylim=c(minEffectSize,maxEffectSize))  + 
    theme_bw()+theme(legend.position = "none")
  #+theme(legend.position = c(0.5,0.8), legend.box="horizontal",legend.background = element_rect(fill="white",linewidth=0.5,linetype="solid",color="black"))
  ggsave(paste("figs/",name,"effectSize.png",sep=''),width=1000, height=1000,unit="px",dpi=200)
  #power
  library(plyr)
  dataOfSimsSummarizedReps=ddply(dataOfSimsSpecificEffect,
                                 .(intercept,type,numberOfTrials,effects,N),
                                 summarize,
                                 reps=sum(singularFit==0 & !is.na(pValue)),
                                 singularFits=sum(singularFit==1 | is.na(pValue)),
                                 powerPositive=sum(!is.na(pValue) & pValue<.05 & effectSize>0 & !is.na(effectSize) & singularFit==0)/reps,
                                 powerNegative=sum(!is.na(pValue) & pValue<.05 & effectSize<0 & !is.na(effectSize) & singularFit==0)/reps,
                                 nsEffect=sum(!is.na(pValue) & pValue>=.05 & singularFit==0)/reps)
  #power to detect effect in positive direction
  ggplot(dataOfSimsSummarizedReps,aes(x=toAcc(intercept),y=powerPositive,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    #stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    geom_hline(yintercept=0.025) + 
    labs(y="Power", x="Intercept", color = "Type of Analysis", shape = "Number of Trials") + coord_cartesian(ylim=c(0,maxPowerPos))  + 
    theme_bw()+theme(legend.position = "none")
  #+theme(legend.position = c(0.5,0.3), legend.background = element_rect(fill="white",linewidth=0.5,linetype="solid",color="black"))
  ggsave(paste("figs/",name,"powerPositive.png",sep=''),width=1000, height=1000,unit="px",dpi=200)
  #power to detect effect in negative direction
  ggplot(dataOfSimsSummarizedReps,aes(x=toAcc(intercept),y=powerNegative,color=type,shape=as.factor(numberOfTrials))) +
    stat_summary(na.rm=TRUE, fun=mean, geom="line") +
    stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
    #stat_summary(fun.data=mean_se,geom="errorbar",position = "dodge",aes(linetype=NULL)) +
    labs(y="Power", x="Intercept", color = "Type of Analysis", shape = "Number of Trials") + coord_cartesian(ylim=c(0,maxPowerNeg))  + 
    geom_hline(yintercept=0.025) + 
    theme_bw()+theme(legend.position = "none")
  #+theme(legend.position = c(0.5,0.7), legend.background = element_rect(fill="white",linewidth=0.5,linetype="solid",color="black"))
  #+theme(legend.position = c(0.8,0.8), legend.background = element_rect(fill="white",linewidth=0.5,linetype="solid",color="black"))
  ggsave(paste("figs/",name,"powerNegative.png",sep=''),width=1000, height=1000,unit="px",dpi=200)
}

###script

#generate random seed (this should be random enough)
#sample(0:100000,1)
#52766
set.seed(52766)
#different combinations of participants and measurements per participant (Ns need to be multiple of numberOfTrialsVector*2)
Ns=c(4160) #4160=40 trials*52 participants*2conditions> 1600 Brysbaert & Stevens (2018), 52 participants for within subjects (Brysbaert, 2019) (interaction would need more) despite of arbitrary effect sizes here to emphasize this value as minimal number of measurements
numberOfTrialsVector=c(10,40)
intercepts=toLogOdds(c(0.01,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,0.99))
reps=1000
#effect sizes of 1 and positive/zero/negative interaction, no random error
#positive interaction
dataOfSims11110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=1,randomIntercept=1,randomError=0)
#no interaction
dataOfSims11010=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=0,randomIntercept=1,randomError=0)
#negative interaction
dataOfSims11_110=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=1,withinEffectSize=1,interactionEffectSize=-1,randomIntercept=1,randomError=0)

#generate some plots
saveplots(dataOfSims11110, "11110")
savePlotsSpecific(dataOfSims11110,"Interactions/11110","factor1:factor2",0,2,1)
savePlotsSpecific(dataOfSims11110,"factor1/11110","factor1",0,2,1)
savePlotsSpecific(dataOfSims11110,"factor2/11110","factor2",0,2,1)
saveplots(dataOfSims11010, "11010")
savePlotsSpecific(dataOfSims11010,"Interactions/11010","factor1:factor2")
savePlotsSpecific(dataOfSims11010,"factor1/11010","factor1",0,2,1)
savePlotsSpecific(dataOfSims11010,"factor2/11010","factor2",0,2,1)
saveplots(dataOfSims11_110, "11-110")
savePlotsSpecific(dataOfSims11_110,"Interactions/11-110","factor1:factor2",-2,0,0.25,1)
savePlotsSpecific(dataOfSims11_110,"factor1/11-110","factor1",0,2,1)
savePlotsSpecific(dataOfSims11_110,"factor2/11-110","factor2",0,2,1)

#effect sizes of 0.5 and positive/zero/negative interaction, no random error
#positive interaction
dataOfSims55550=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=.5,withinEffectSize=.5,interactionEffectSize=.5,randomIntercept=.5,randomError=0)
#no interaction
dataOfSims55050=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=.5,withinEffectSize=.5,interactionEffectSize=0,randomIntercept=.5,randomError=0)
#negative interaction
dataOfSims55_550=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=.5,withinEffectSize=.5,interactionEffectSize=-.5,randomIntercept=.5,randomError=0)

#generate some plots
saveplots(dataOfSims55550, "55550")
savePlotsSpecific(dataOfSims55550,"Interactions/55550","factor1:factor2",0,1,1)
savePlotsSpecific(dataOfSims55550,"factor1/55550","factor1",0,1,1)
savePlotsSpecific(dataOfSims55550,"factor2/55550","factor2",0,1,1)
saveplots(dataOfSims55050, "55050")
savePlotsSpecific(dataOfSims55050,"Interactions/55050","factor1:factor2",-0.5,0.5)
savePlotsSpecific(dataOfSims55050,"factor1/55050","factor1",0,1,1)
savePlotsSpecific(dataOfSims55050,"factor2/55050","factor2",0,1,1)
saveplots(dataOfSims55_550, "55-550")
savePlotsSpecific(dataOfSims55_550,"Interactions/55-550","factor1:factor2",-1,0,0.25,1)
savePlotsSpecific(dataOfSims55_550,"factor1/55-550","factor1",0,1,1)
savePlotsSpecific(dataOfSims55_550,"factor2/55-550","factor2",0,1,1)

#add random error
#no interaction
dataOfSims55055=randSim(Ns,numberOfTrialsVector,intercepts,reps,betweenEffectSize=.5,withinEffectSize=.5,interactionEffectSize=0,randomIntercept=.5,randomError=.5)
saveplots(dataOfSims55055, "55055")
savePlotsSpecific(dataOfSims55055,"Interactions/55055")

#combine images
#singular fits
combineImages(c("figs/11-110singularFits.png",
                "figs/11010singularFits.png",
                "figs/11110singularFits.png",
                "figs/55-550singularFits.png",
                "figs/55050singularFits.png",
                "figs/55550singularFits.png"),2,3,"figs/singularFits.png")
#effect size interaction
combineImages(c("figs/Interactions/11-110effectSize.png",
                "figs/Interactions/11010effectSize.png",
                "figs/Interactions/11110effectSize.png",
                "figs/Interactions/55-550effectSize.png",
                "figs/Interactions/55050effectSize.png",
                "figs/Interactions/55550effectSize.png"),2,3,"figs/effectSizeInteraction.png")
#power negative interactions
combineImages(c("figs/Interactions/11-110powerNegative.png",
                "figs/Interactions/11010powerNegative.png",
                "figs/Interactions/11110powerNegative.png",
                "figs/Interactions/55-550powerNegative.png",
                "figs/Interactions/55050powerNegative.png",
                "figs/Interactions/55550powerNegative.png"),2,3,"figs/powerNegativeInteraction.png")
#power positive interaction
combineImages(c("figs/Interactions/11-110powerPositive.png",
                "figs/Interactions/11010powerPositive.png",
                "figs/Interactions/11110powerPositive.png",
                "figs/Interactions/55-550powerPositive.png",
                "figs/Interactions/55050powerPositive.png",
                "figs/Interactions/55550powerPositive.png"),2,3,"figs/powerPositiveInteraction.png")
#effect size main effects
combineImages(c("figs/factor1/11-110effectSize.png",
                "figs/factor1/11010effectSize.png",
                "figs/factor1/11110effectSize.png",
                "figs/factor1/55-550effectSize.png",
                "figs/factor1/55050effectSize.png",
                "figs/factor1/55550effectSize.png"),2,3,"figs/effectSizeFactor1.png")
combineImages(c("figs/factor2/11-110effectSize.png",
                "figs/factor2/11010effectSize.png",
                "figs/factor2/11110effectSize.png",
                "figs/factor2/55-550effectSize.png",
                "figs/factor2/55050effectSize.png",
                "figs/factor2/55550effectSize.png"),2,3,"figs/effectSizeFactor2.png")
#power main effects
combineImages(c("figs/factor1/11-110powerPositive.png",
                "figs/factor1/11010powerPositive.png",
                "figs/factor1/11110powerPositive.png",
                "figs/factor1/55-550powerPositive.png",
                "figs/factor1/55050powerPositive.png",
                "figs/factor1/55550powerPositive.png"),2,3,"figs/powerPositiveFactor1.png")
combineImages(c("figs/factor2/11-110powerPositive.png",
                "figs/factor2/11010powerPositive.png",
                "figs/factor2/11110powerPositive.png",
                "figs/factor2/55-550powerPositive.png",
                "figs/factor2/55050powerPositive.png",
                "figs/factor2/55550powerPositive.png"),2,3,"figs/powerPositiveFactor2.png")
