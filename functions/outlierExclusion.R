### Effects of outlier exclusion due to accuracy on participant numbers and effect size
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

library(pwr)
library(ggplot2)

#function to calculate effects
#n - number of trials
#p - probability to solve trials
#c - chance level
#x - proportion of non-guessers
calculateNewEffects=function(n,p,c,x){
  Nsolvers=c(t((1-pbinom(n*c,n,p))%*%t(x))) #these depend on p (vector) -> convert matrix to vector
  Nguessers=rep((1-pbinom(n*c,n,c))*(1-x),length(p))  #these do not depend on p (no vector) -> repeat to same length
  #calculate new values
  N=Nsolvers+Nguessers
  xn=Nsolvers/N
  powerNew=pwr.t.test(n=round(1000*N),d=0.1*xn)
  powerOriginal=pwr.t.test(n=1000,d=0.1*x)
  powerProportion=powerNew$power/powerOriginal$power
  #data frame for saving (3 values for each p/x combination)
  dataframe=data.frame(numberOfTrials=n,
                       p=rep(p,3,each=length(x)),
                       proportionNonGuessers=as.factor(rep(x,3*length(p))),
                       type=rep(c("N","x","power"),each=length(p)*length(x)),
                       value=c(N,xn,powerProportion))
  return(dataframe)
}
p=0.5+c(1:49)/100 #probability to solve tasks
c=0.5 #chance level
x=c(0.8,0.95,1) #proportion of non-guessers
dat10=calculateNewEffects(10,p,c,x)
dat20=calculateNewEffects(20,p,c,x)
dat40=calculateNewEffects(40,p,c,x)
datFull=rbind(dat10,dat20,dat40)


ggplot(datFull,aes(x=p,y=value,color=type,shape=proportionNonGuessers)) +
  stat_summary(na.rm=TRUE, fun=mean, geom="line") +
  stat_summary(na.rm=TRUE, fun=mean, geom="point", size=2) +
  facet_wrap(~numberOfTrials) + 
  labs(y="proportion", x="p") +
  theme_classic()
ggsave("figs/powerComparison.png",width=1920, height=1080,unit="px",dpi=200)

