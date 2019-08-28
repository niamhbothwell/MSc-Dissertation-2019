## Provisioning rate with linear mixed effect models (LMEs) ##
## Key of data ##
#year = year of breeding season
#uniquenest = the nest number given a random number for each nest 
#birdID_male = the name of the male individual
#prov_rate_male = provisioning rate of the focal male estimated by the number of times successful fed their chicks divided by the number of hours observed   
#prop_prov_male = the given proportion of the focal males provisioning calculated by the individuals provisioning rate/total provisioning rate of both parents times 100
#birdID_female = the name of the female individual
#prov_rate_female = provisioning rate of the focal female estimated by the number of times successful fed their chicks divided by the number of hours observed.  
#prop_prov_female = the given proportion of the focal females provisioning calculated by the individuals provisioning rate/total provisioning rate of both parents times 100
#Total_prov_rate = the total provisioning rate of both the parents
#male_age = the age of the male during that year 
#average_age_male = the mean of the male calculated by the adding the ages together and dividing them by the number of times each individual was observed across years 
#delta_age_male = calculated by the difference of the focal male average age.  
#female_age = the age of the focal female during that year 
#average_age_female = the mean of the female calculated by the adding the ages together and dividing them by the number of times each individual was observed across years 
#delta_age_female = calculated by the difference of the focal females average age.  
#pairbond = the duration the pair have bred together. e.g. 0 first year of breeding together, 1 the second year of breeding, etc..
#bs = the broad size, the number of hatched. 
#age-difference = the difference between the pairs ages 


# ******************************************************************************************************************************************************************************* 

## 1. Data tidying ##

## clear workspace ##
rm(list = ls())

## set working directory ## 
setwd("~/Documents/Masters EEC/Research project")
library(readxl)
prov<- read_excel("Data_Niamh .xlsx")

prov<- read.table("Data_NB_.csv",h=T,sep=",") # N = 197
summary(prov)
head(prov)
str(prov) 

## install packages ##
install.packages("dplyr")
require("dplyr")
install.packages("lme4")
require(lme4)
install.packages("lmerTest")
require(lmerTest)
install.packages("lmtest")
require(lmtest)
library(lattice)
library(fitdistrplus)
require(fitdistrplus)

# subset NA's # 
provs<- subset(na.omit(prov))  # N = 160

# convert brood size into a factor #
provs$bs_fac<- as.factor(provs$bs)


# ****************************************************************************************************************************************************************************************

## HYPOTHESIS I ##
## The provisioning rate will be more equally shared between partners as the length of their pair-bond increases and that pair-bond duration is more important than age per se ##
# I would expect the effect of pairbond to be significantly negative if the hypothesis holds #

# dependent variable = absolute difference in proportion of provisioning between male and female #
# random effects = year, bird identity #
# fixed effects = pairbond duration, pairbond duration^2, average age, delta age, average age*delta age, absolute age difference, absolute age difference^2, brood size #

# calculate the absolute difference in proportion of provisioning between male and female #

# proportion of provisioning rate for female #
provs$prop_prov_rate_female<- provs$prov_rate_female/provs$Total_prov_rate
summary(provs$prop_prov_rate_female) # varies between 0 and 1 so it is a proportion
hist(provs$prop_prov_rate_female)
descdist(provs$prop_prov_rate_female, discrete = T) # poisson or normal

# proportion of provisioning rate for male #
provs$prop_prov_rate_male<- provs$prov_rate_male/provs$Total_prov_rate
summary(provs$prop_prov_rate_male) # varies between 0 and 1 so it is a proportion
hist(provs$prop_prov_rate_male)
descdist(provs$prop_prov_rate_male, discrete = T) # poisson or normal

# transform to absolute difference to get whole numbers and no negative values #
provs$abs_diff_prop_prov<- abs(provs$prop_prov_rate_male - provs$prop_prov_rate_female)
summary(provs$abs_diff_prop_prov)
hist(provs$abs_diff_prop_prov)
descdist(provs$abs_diff_prop_prov, discrete = T)

## difference in proportion of provisioning between male and female ##
## negative if the young male does less, positive if the older male does more ##
provs$diff_prop_prov = provs$prop_prov_rate_male - provs$prop_prov_rate_female
summary(provs$diff_prop_prov)
hist(provs$diff_prop_prov) #left-skewed 
descdist(provs$diff_prop_prov, discrete = T) # the better is the normal distribution

# age difference between male and female (negative if the male is younger, positive if the male is older) #
names(provs)
provs$age_diff = provs$male_age - provs$female_age
summary(provs$age_diff) 

#absolute age difference between partners #
provs$abs_age_diff<-abs(provs$age_diff)
summary(provs$abs_age_diff)
hist(provs$abs_age_diff)
provs$abs_age_diff<-abs(provs$age_difference)
summary(provs$abs_age_diff) #range: 0-13 years 
hist(provs$abs_age_diff)

## MODEL 1 ##

# absolute difference in proportion of provisioning in females model #
m1<- lmer(abs_diff_prop_prov ~ pairbond
           + I(pairbond^2) #quadratic effect 
           + average_age_female
           + delta_age_female
           + average_age_female:delta_age_female
           + abs_age_diff
           + I(abs_age_diff^2) #quadratic effect  
           + bs_fac
           + (1|birdID_female)
           + (1|year)
           , data = provs)
summary(m1)

## "FINAL" MODEL using backwards elimination ##

m1a<- lmer(abs_diff_prop_prov ~ #pairbond
    #     + I(pairbond^2) #quadratic effect in R 
          + average_age_female
          + delta_age_female
    #     + average_age_female:delta_age_female
    #     + abs_age_diff
    #     + I(abs_age_diff^2) 
          + bs_fac
          + (1|birdID_female)
          + (1|year)
          , data = provs)
summary(m1a)

## result shows only average female age significant ##

## for pairbond duration against absolute difference in the proportion of provisioning between males and females ##

# total of newly-formed paribond and previous pair-bonds between 2015-2019 #

sum(provs$pairbond == "0") #6 
sum(provs$pairbond > "0") #154

# length of pairbond #

sum(provs$pairbond < "5") # 137 pair-bonds were together for less than 5 years 
sum(provs$pairbond >= "5") #23 were in a pair-bond for 5 years and above 

# set the graph #
min(provs$pairbond, na.rm=T) # 0
max(provs$pairbond, na.rm=T) # 18

min(provs$abs_diff_prop_prov, na.rm=T) # 0
max(provs$abs_diff_prop_prov, na.rm=T) # 1

# plot #
par(mar=c(6.5, 8.1, 1.1, 1.1) + 0.5)
plot(provs$pairbond, provs$abs_diff_prop_prov, 
     pch=16, cex=1.0, col="black", xlab="", ylab="", 
     axes = F, tck=-0.015, cex.axis = 1.3, ylim=c(-0.1,1), xlim = c(0,18))


# add axis #
axis(side = 2, cex.axis=1.0, las=1)
axis(side = 1, cex.axis=1.0, las=1)
box()


# add axis labels #
mtext("Absolute difference in provisioning \nbetween partners (proportion)",2,4,cex = 1.0)
mtext("Pairbond length (years)",1,3,cex = 1.0) 

# plot within- and between- average and delta for females and males#

# graph for female average and delta age #

# x= female age #
# y = absolute difference in the proportion of provisioning rate between males and females #

# get the fit for the average age effect for females #
library(effects)
S0<- effect("average_age_female", m1a, se=T, partial.residuals=T) 
plot(S0)
S0<- as.data.frame(S0)
S0

# set graph minimum and maximum #
min(provs$average_age_female, na.rm=T) # 3
max(provs$average_age_female, na.rm=T) # 22.5

min(provs$abs_diff_prop_prov, na.rm=T) # 0
max(provs$abs_diff_prop_prov, na.rm=T) # 1

# plot the line #
par(mar=c(6.5, 8.1, 1.1, 1.1) + 0.5)
plot(S0$average_age_female, S0$fit, type="l", lwd=2, xlim= c(3, 23), 
     ylim=c(-0.1,1),xlab="", ylab="", axes = F, tck=-0.015, cex.axis = 1.3)

# SE for model "m1a" #
lines(S0$average_age_female,S0$upper, lty="dashed", lwd=1)
lines(S0$average_age_female,S0$lower, lty="dashed", lwd=1)

# axis #
axis(side = 2, cex.axis=1.0, las=1)
axis(side = 1, cex.axis=1.0, las=1)
box()

# data points = female age #
points(provs$female_age, provs$abs_diff_prop_prov, pch=16, cex=1.0, col="black")

# axis labels #
mtext("Absolute difference in provisioning \nbetween partners (proportion)",2,4,cex = 1.0)
mtext("Female age (years)",1,3,cex = 1.0) 

# fit of delta_age #
# transform the delta age variable to get it in the same scale as female age and average female age #
provs$delta_av_female<- provs$average_age_female + provs$delta_age_female
head(provs)
summary(provs)

# create a new model with this transformed variable #
m1b<- lmer(abs_diff_prop_prov ~ pairbond
           + average_age_female
           + delta_av_female
           + bs_fac
           + (1|birdID_female)
           + (1|year)
           , data = provs)
summary(m1b)

# fit of the delta_av_female #
S1<- effect("delta_av_female", m1b, se=T, partial.residuals=T)
S1<- as.data.frame(S1)
S1

# add the model line to existing plot #
lines(S1$delta_av,S1$fit, col = "darkgrey", lwd=2)
lines(S1$delta_av_female,S1$upper, lty="dashed", col = "darkgrey", lwd=1)
lines(S1$delta_av_female,S1$lower, lty="dashed", col = "darkgrey", lwd=1)

## Absolute difference in proportion of provisioning in males model ##

# FULL MODEL #
m2<- lmer(abs_diff_prop_prov ~ pairbond
          + I(pairbond^2) #quadratic effect 
          + average_age_male
          + delta_age_male
          + average_age_male:delta_age_male
          + abs_age_diff
          + I(abs_age_diff^2) #quadratic effect 
          + bs_fac
          + (1|birdID_male)
          + (1|year)
          , data = provs)
summary(m2)

## "FINAL" MODEL using backwards elimination ## 

m2a<- lmer(abs_diff_prop_prov ~ #pairbond
     #    + I(pairbond^2) 
          + average_age_male
          + delta_age_male
    #     + average_age_male:delta_age_male
          + abs_age_diff
          + I(abs_age_diff^2) 
          + bs_fac
          + (1|birdID_male)
          + (1|year)
          , data = provs)
summary(m2a)

# no significance #

# graph for male average and delta age #

# x= male age #
# y = absolute difference in the proportion of provisioning rate between males and females #

# get the fit for the average age effect for males #
S3<- effect("average_age_male", m2a, se=T, partial.residuals=T) # always use the final model here!
plot(S3)
S3<- as.data.frame(S3)
S3

# set graph minimum and maximum #
min(provs$average_age_male, na.rm=T) # 3
max(provs$average_age_male, na.rm=T) # 23

min(provs$abs_diff_prop_prov, na.rm=T) # 0
max(provs$abs_diff_prop_prov, na.rm=T) # 1

# plot the line #
par(mar=c(6.5, 8.1, 1.1, 1.1) + 0.5)
plot(S3$average_age_male, S3$fit, type="l", lwd=2, xlim= c(3, 23), 
     ylim=c(-0.1,1),xlab="", ylab="", axes = F, tck=-0.015, cex.axis = 1.3)

# add SE for model "m2a" #
lines(S3$average_age_male,S3$upper, lty="dashed", lwd=1)
lines(S3$average_age_male,S3$lower, lty="dashed", lwd=1)

# add axis #
axis(side = 2, cex.axis=1.0, las=1)
axis(side = 1, cex.axis=1.0, las=1)
box()

# add data points = male age #
points(provs$male_age, provs$abs_diff_prop_prov, pch=16, cex=1.0, col="black")

# add axis labels #
mtext("Absolute difference in provisioning \nbetween partners (proportion)",2,4,cex = 1.0)
mtext("Male age (years)",1,3,cex = 1.0) 

# fit of delta_age for males #
# transform the delta age variable to get it in the same scale as male age and average male age #
provs$delta_av_male<-provs$average_age_male + provs$delta_age_male
head(provs)
summary(provs)

# create a new model with this transformed variable #
mod4<- lmer(abs_diff_prop_prov ~ pairbond
              + I(pairbond^2) #quadratic effect in R 
              + abs_age_diff
              + delta_av_male
              + I(abs_age_diff^2) #quadratic effect in R 
              + bs_fac
              + (1|birdID_male)
              + (1|year)
              , data = provs)
summary(mod4)

# fit of the delta_av_male #
S4<- effect("delta_av_male", mod4, se=T, partial.residuals=T)
S4<- as.data.frame(S4)
S4

# add the model line in the existing plot #
lines(S4$delta_av,S4$fit, col = "darkgrey", lwd=2)
lines(S4$delta_av_male,S4$upper, lty="dashed", col = "darkgrey", lwd=1)
lines(S4$delta_av_male,S4$lower, lty="dashed", col = "darkgrey", lwd=1)


# ****************************************************************************************************************************************************************************************

## HYPOTHESIS II ## 

## difference in the proportion of provisioning vs absolute age difference ##
# I would expect the effect of the age difference to be significantly positive if the hypothesis holds #
# from model 1 and 2 shows age difference is not statistically significant #

# set the graph #
min(provs$abs_age_diff, na.rm=T) # 0
max(provs$abs_age_diff, na.rm=T) # 13

min(provs$diff_prop_prov, na.rm=T) # -1
max(provs$diff_prop_prov, na.rm=T) # 1

# plot #
par(mar=c(6.5, 8.1, 1.1, 1.1) + 0.5)
plot(provs$abs_age_diff, provs$abs_diff_prop_prov, 
     xlab="", ylab="", pch=16, cex=1.0, col="black",
     axes = F, tck=-0.015, cex.axis = 1.3)

# add axis #
axis(side = 2, cex.axis=1.0, las=1)
axis(side = 1, cex.axis=1.0, las=1)
box()

# add data points = age difference #
points(provs$abs_age_diff, provs$diff_prop_prov, pch=16, cex=1.0, col="black")

# add axis labels #
mtext("Difference in provisioning \nbetween partners (proportion)",2,4,cex = 1.0)
mtext("Absolute age difference",1,3,cex = 1.0) 

# the proportion of pairs who's age difference was between 0-1#
sum(provs$age_difference=="0") #57
sum(provs$age_difference=="1") #57

# age difference between 0-1 57+57/160= 71.77% #
