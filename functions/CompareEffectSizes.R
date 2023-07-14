### Comparing additive effects and interactions
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

#parameters
chanceLevel=.5
baseAcc=c(.51,.55,.6,.65,.7,.75,.8,.85,.89) #base value for accuracy score in [chanceLevel,1]
effectSize=.1 #change in accuracy from base accuracy

###calculations (no interaction)
##normal distribution
#simple addition of probabilities
accuracyOneEffect=baseAcc+effectSize
accuracyTwoEffects=accuracyOneEffect+effectSize

##uncorrected binomial distribution
#calculate log odds for corresponding change by one effect
effectSizeLogOdds=toLogOdds(baseAcc+effectSize)-toLogOdds(baseAcc)
#additions in log odds
baseLogOdds=toLogOdds(baseAcc)
logOddsOneEffect=baseLogOdds+effectSizeLogOdds
logOddsTwoEffects=logOddsOneEffect+effectSizeLogOdds
#back to probability/accuracy
accuracyOneEffectBinomial=toAcc(logOddsOneEffect)
accuracyTwoEffectsBinomial=toAcc(logOddsTwoEffects)

##guessing corrected
#calculate corrected scores for effects
baseAccCorrected=correctChanceLevel(baseAcc,chanceLevel)
baseAccCorrectedOneEffect=correctChanceLevel(baseAcc+effectSize,chanceLevel)
#calculate log odds for corresponding change by one effect
effectSizeLogOddsCorrected=toLogOdds(baseAccCorrectedOneEffect)-toLogOdds(baseAccCorrected)
#additions in log odds
baseLogOddsCorrected=toLogOdds(baseAccCorrected)
logOddsOneEffectCorrected=baseLogOddsCorrected+effectSizeLogOddsCorrected
logOddsTwoEffectsCorrected=logOddsOneEffectCorrected+effectSizeLogOddsCorrected
#back to probability/accuracy and add chance level again
accuracyOneEffectCorrected=addChanceLevel(toAcc(logOddsOneEffectCorrected),chanceLevel)
accuracyTwoEffectsCorrected=addChanceLevel(toAcc(logOddsTwoEffectsCorrected),chanceLevel)


#one effect (these should all be the same)
print(accuracyOneEffect)
print(accuracyOneEffectBinomial)
print(accuracyOneEffectCorrected)

#two effects
print(accuracyTwoEffects)
print(accuracyTwoEffectsBinomial)
print(accuracyTwoEffectsCorrected)
