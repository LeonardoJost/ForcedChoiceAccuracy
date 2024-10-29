### Effects of outlier exclusion due to accuracy on participant numbers and effect size
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

library(pwr)
library(ggplot2)

#function to calculate effects
#n - number of trials
#p - probability to solve trials
#c - chance level
#x - proportion of non-guessers
#e - exclusion level (typically chance level)
calculateNewEffects=function(n,p,c,g,e=c){
  #calculate the number of persons performing above chance level
  Nsolvers=c(t((1-pbinom(n*e,n,p))%*%t(g))) #these depend on p (vector) -> convert matrix to vector
  Nguessers=rep((1-pbinom(n*e,n,c))*(1-g),length(p))  #these do not depend on p (no vector) -> repeat to same length
  #calculate new values
  Nnew=Nsolvers+Nguessers #new number of persons after exclusion
  gn=Nsolvers/Nnew #new proportion of non-guessers
  #calculate number of participants necessary for power of .8 in original sample
  NOriginal=rep(0,length(g))
  for(i in 1:length(g)){
    NOriginal[i]=pwr.t.test(power=0.8,d=g[i],type="two.sample",alternative="two.sided")$n
  }
  #new power
  powerNew=pwr.t.test(n=NOriginal*Nnew,d=gn,type="two.sample",alternative="two.sided")$power/0.8
  #data frame for saving (3 values for each p/x combination)
  dataframe=data.frame(numberOfTrials=n,
                       p=rep(p,3,each=length(g)),
                       proportionNonGuessers=as.factor(rep(g,3*length(p))),
                       type=rep(c("N","d","power"),each=length(p)*length(g)),
                       value=c(Nnew,gn,powerNew))
  return(dataframe)
}
p=0.5+c(1:49)/100 #probability to solve tasks
c=0.5 #chance level
g=c(0.8,0.95,1) #proportion of non-guessers
dat10=calculateNewEffects(10,p,c,g)
dat20=calculateNewEffects(20,p,c,g)
dat40=calculateNewEffects(40,p,c,g)
dat406=calculateNewEffects(40,p,c,g,0.6)
dat406$numberOfTrials="40-0.6"
datFull=rbind(dat10,dat20,dat40,dat406)


ggplot(datFull,aes(x=p,y=value,color=type,shape=proportionNonGuessers)) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  facet_wrap(~numberOfTrials) + 
  labs(y="Proportion", x="p", shape="Proportion of Non-guessers", color="Type of Measure") +
  #geom_hline(yintercept=0.8) +
  theme_bw()+theme(legend.position = c(0.8,0.65), legend.background = element_rect(fill="white",linewidth=0.5,linetype="solid",color="black"))
ggsave("figs/powerComparison.png",width=2000, height=2000,unit="px",dpi=200)

