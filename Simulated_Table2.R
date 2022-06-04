##### 0. REQUIRED PACKAGES

#install.packages("tidyverse")

library(tidyverse)
library(survival)

##### 1. Simulate covariates, survival time, and event

set.seed(12345)

n = 5000 # number of individuals in the study

##### 1(a) Simulate covariates

AGE = round(runif(n,18,24))
GENDER = sample(c("Female", "Male", "N/A"), replace = T, prob = c(.510, .487, .003), size = n)
AFFILIATION = sample(c("On-campus", "Off-campus"), replace = T, prob = c(.323, .677), size = n)
race_eth = sample(c(1:4), replace = T, prob = c(.796, .060, .066, .078), size = n)
CONDITION = rbinom(n,size=1,prob=.053)
CONDITIONIMMUN = rbinom(n, size=1, prob=.019)
MEDICATIONS = rbinom(n, size = 1, prob = .023)
NICOTINE = rbinom(n, size = 1, prob = .059)
prev_positive = rbinom(n, size = 1, prob = .215)
manufacturer = sample(c(0:5), replace = T, prob = c(.393, .003, .003, .215, .342, .045), size = n)
manufacturer.2 = ifelse(manufacturer == 0, NA, 
                        ifelse(manufacturer == 5, 3, 
                               ifelse(manufacturer %in% c(1,3),1,2)))
vaccinated = ifelse(manufacturer %in% c(1,2), 1, ifelse(manufacturer == 0, 0, 2))
first.dose.jandj = round(rnorm(sum(manufacturer==5),-65,50))
first.dose.moderna = round(rnorm(sum(manufacturer %in% c(1,3)),-91,45))
first.dose.pfizer = round(rnorm(sum(manufacturer %in% c(2,4)),-82,50))
first.dose = rep(NA,n)
first.dose[manufacturer == 5] = first.dose.jandj
first.dose[manufacturer %in% c(1,3)] = first.dose.moderna
first.dose[manufacturer %in% c(2,4)] = first.dose.pfizer
second.dose = rep(NA,n)
second.dose[manufacturer == 3] = first.dose[manufacturer == 3] + 28
second.dose[manufacturer == 4] = first.dose[manufacturer == 4] + 28


X = data.frame(age = AGE, 
               gender1 = ifelse(GENDER == "Male", 1, 0),
               gender2 = ifelse(GENDER == "N/A", 1, 0),
               affl = ifelse(AFFILIATION == 'On-campus', 1, 0),
               race1 = ifelse(race_eth == 2, 1, 0),
               race2 = ifelse(race_eth == 3, 1, 0),
               race3 = ifelse(race_eth == 4, 1, 0),
               condition = CONDITION,
               immune = CONDITIONIMMUN,
               medic = MEDICATIONS,
               nicotine = NICOTINE,
               prev_positive = prev_positive,
               vaccinated1 = ifelse(vaccinated == 1, 1, 0),
               vaccinated2 = ifelse(vaccinated == 2, 1, 0))
X.mat = as.matrix(X)

##### 1(b) Simulate survival time and event based on regression coefficient beta.true

beta.true = c(-0.25, -0.16, -0.24, -0.27, -0.01, -0.14, -0.36, -0.34, 0.08, 0.36, 0.24, 1.47, 0.13, 1.23)

v = runif(n)
lambda = 30
rho = 1
rateC = .02
Tlat = (-log(v)/(lambda * exp(c(X.mat %*% beta.true))))^(1/rho)
Tlat = max(Tlat) - Tlat
Tlat = ifelse(Tlat == 0, 0.5, Tlat)

C = rexp(n, rate = rateC)
time = ceiling(pmin(Tlat,C))
status = as.numeric(Tlat <= C)

##### 1(c) Convert data to time-varying form with "vaccinated" as the time-varying covariate.
#####      The data frame for fitting the Cox model is survival.data

temp = data.frame(ID = 1:n, first.dose = first.dose, second.dose = second.dose, survival.time = time, 
                  status = status, manufacturer.2 = manufacturer.2)

temp = temp %>% mutate(manufacturer = ifelse(is.na(manufacturer.2) | survival.time <= first.dose, 0, 
                                             ifelse(manufacturer.2 == 3 & survival.time > first.dose, 3, 1)))
temp$manufacturer = ifelse(temp$manufacturer==1 & temp$manufacturer.2 == 2, 2, temp$manufacturer)

new.temp = tmerge(temp,temp, id = ID, status = event(survival.time,status))
new.temp = tmerge(new.temp, temp, id = ID, vaccinated = cumtdc(first.dose))
new.temp = tmerge(new.temp, temp, id = ID, vaccinated = cumtdc(second.dose))

new.temp = new.temp %>% mutate(vaccinated = ifelse(manufacturer.2 == 3 & vaccinated == 1, 2, vaccinated))
new.temp = new.temp %>% select(ID,tstart,tstop,status,vaccinated)

covariates = data.frame(ID = 1:5000, AGE = AGE, GENDER = GENDER, AFFILIATION = AFFILIATION, RACE = race_eth, CONDITION = CONDITION,
                        CONDITIONIMMUN = CONDITIONIMMUN, MEDICATIONS = MEDICATIONS, NICOTINE = NICOTINE, 
                        prev_positive = prev_positive, manufacturer = manufacturer.2)
covariates$manufacturer = ifelse(is.na(covariates$manufacturer), 0, covariates$manufacturer)
covariates2 = data.frame(covariates, vaccinated = vaccinated, status = status)


survival.data = left_join(new.temp,covariates,by = "ID")

survival.data = survival.data %>% mutate(manufacturer_v2 = case_when(vaccinated == 0 ~ 0,
                                                                     vaccinated == 1 & manufacturer == 1 ~ 1,
                                                                     vaccinated == 1 & manufacturer == 2 ~ 2,
                                                                     vaccinated == 2 & manufacturer == 1 ~ 3,
                                                                     vaccinated == 2 & manufacturer == 2 ~ 4,
                                                                     vaccinated == 2 & manufacturer == 3 ~ 5))

survival.data = survival.data %>% mutate_at(vars(c('GENDER','AFFILIATION','CONDITION','CONDITIONIMMUN',
                                                   'MEDICATIONS','NICOTINE','RACE',
                                                   'prev_positive','vaccinated','manufacturer','manufacturer_v2')), ~as.factor(.))

##################################################################################################################################

##### 2. Fit Cox model, estimate protection, format Table 2

##### 2(a) Function to extract point and interval estimates of protection

EstCI = function(output, terms, new.x){
    coef = output$coefficients[terms]
    var = vcov(output)[terms,terms]
    
    point = format(round(1 - exp(new.x%*%coef),3)*100,nsmall=1, trim=T)
    SE = sqrt(new.x%*%var%*%new.x)
    
    CI = format(round(c(1 - exp(new.x%*%coef + 1.96*SE),  1 - exp(new.x%*%coef - 1.96*SE)),3)*100,nsmall = 1,trim=T)
    est = paste0(point, '% (', CI[1],'-',CI[2],')')
    return(est)
}

##### 2(b) Headers in Table 2
#####      "top" indicates the models w/o interaction between vaccination and previous infection
#####      "bottom" indicates the models with interaction between vaccination and previous infection

top.headers = c("Vaccination Protection", "# of Individuals", "# Positive (%)", "Protection: % (95% CI)")
bottom.headers = c("Protection by Vaccination and Previous Infection History", rep(" ", 3))
status.top = c("Unvaccinated", "Fully Vaccinated", "...mRNA-1273", "...BNT162b2", "...Ad26.COV2.S")
status.bottom = c("No protection", "Fully vaccinated", 
                  "...No previous infection", "......mRNA-1273", "......BNT162b2", "......Ad26.COV2.S",
                  "...Previous infection", "......mRNA-1273", "......BNT162b2", "......Ad26.COV2.S","Previous infection only")

##### 2(c) Sample summary statistics

top.summary = table(filter(covariates2,vaccinated != 1)$manufacturer, filter(covariates2, vaccinated != 1)$status)

top.tab = matrix(NA, 5, 3)
top.tab[1,1:2] = c(top.summary[1,1]+top.summary[1,2],paste0(top.summary[1,2]," (",round(top.summary[1,2]/(top.summary[1,1]+top.summary[1,2]),3)*100,"%)"))
top.tab[2,1:2] = c(sum(top.summary[2:4,]),paste0(sum(top.summary[2:4,2])," (",round(sum(top.summary[2:4,2])/sum(top.summary[2:4,]),3)*100,"%)"))
top.tab[3,1:2] = c(top.summary[2,1]+top.summary[2,2],paste0(top.summary[2,2]," (",round(top.summary[2,2]/(top.summary[2,1]+top.summary[2,2]),3)*100,"%)"))
top.tab[4,1:2] = c(top.summary[3,1]+top.summary[3,2],paste0(top.summary[3,2]," (",round(top.summary[3,2]/(top.summary[3,1]+top.summary[3,2]),3)*100,"%)"))
top.tab[5,1:2] = c(top.summary[4,1]+top.summary[4,2],paste0(top.summary[4,2]," (",round(top.summary[4,2]/(top.summary[4,1]+top.summary[4,2]),3)*100,"%)"))

bottom.tab = matrix(NA, 11, 3)

no_prot = filter(covariates2, vaccinated == 0, prev_positive == 0)
prev_infect_only = filter(covariates2, vaccinated == 0, prev_positive == 1)
bottom.tab[1,1:2] = c(nrow(no_prot), paste0(sum(no_prot$status), " (", round(sum(no_prot$status)/nrow(no_prot),3)*100,"%)"))
bottom.tab[2,] = rep(" ", 3)

bottom.summary = table(filter(covariates2, vaccinated == 2)$manufacturer, filter(covariates2, vaccinated == 2)$status, filter(covariates2, vaccinated == 2)$prev_positive)

bottom.tab[3,1:2] = c(sum(bottom.summary[,,1]), paste0(sum(bottom.summary[,2,1]), " (", round(sum(bottom.summary[,2,1])/sum(bottom.summary[,,1]),3)*100, "%)"))
bottom.tab[4,1:2] = c(sum(bottom.summary[1,,1]), paste0(sum(bottom.summary[1,2,1]), " (", round(sum(bottom.summary[1,2,1])/sum(bottom.summary[1,,1]),3)*100, "%)"))
bottom.tab[5,1:2] = c(sum(bottom.summary[2,,1]), paste0(sum(bottom.summary[2,2,1]), " (", round(sum(bottom.summary[2,2,1])/sum(bottom.summary[2,,1]),3)*100, "%)"))
bottom.tab[6,1:2] = c(sum(bottom.summary[3,,1]), paste0(sum(bottom.summary[3,2,1]), " (", round(sum(bottom.summary[3,2,1])/sum(bottom.summary[3,,1]),3)*100, "%)"))

bottom.tab[7,1:2] = c(sum(bottom.summary[,,2]), paste0(sum(bottom.summary[,2,2]), " (", round(sum(bottom.summary[,2,2])/sum(bottom.summary[,,2]),3)*100, "%)"))
bottom.tab[8,1:2] = c(sum(bottom.summary[1,,2]), paste0(sum(bottom.summary[1,2,2]), " (", round(sum(bottom.summary[1,2,2])/sum(bottom.summary[1,,2]),3)*100, "%)"))
bottom.tab[9,1:2] = c(sum(bottom.summary[2,,2]), paste0(sum(bottom.summary[2,2,2]), " (", round(sum(bottom.summary[2,2,2])/sum(bottom.summary[2,,2]),3)*100, "%)"))
bottom.tab[10,1:2] = c(sum(bottom.summary[3,,2]), paste0(sum(bottom.summary[3,2,2]), " (", round(sum(bottom.summary[3,2,2])/sum(bottom.summary[3,,2]),3)*100, "%)"))
bottom.tab[11,1:2] = c(nrow(prev_infect_only), paste0(sum(prev_infect_only$status), " (", round(sum(prev_infect_only$status)/nrow(prev_infect_only),3)*100,"%)"))

##### 2(d) Fit Cox model, estimate protection

##### No interaction between vaccination and previous infection

out1 = coxph(Surv(tstart,tstop,status) ~ AGE + GENDER + AFFILIATION + RACE + CONDITION + CONDITIONIMMUN + MEDICATIONS + NICOTINE  + 
                 prev_positive + vaccinated, data = survival.data)
out2 = coxph(Surv(tstart,tstop,status) ~ AGE + GENDER + AFFILIATION + RACE + CONDITION + CONDITIONIMMUN + MEDICATIONS + NICOTINE  + 
                 prev_positive + manufacturer_v2, data = survival.data)

top.tab[1,3] = "Reference"
top.tab[2:5,3] = c(EstCI(out1, "vaccinated2", 1),
                   EstCI(out2, "manufacturer_v23", 1),
                   EstCI(out2, "manufacturer_v24", 1),
                   EstCI(out2, "manufacturer_v25", 1))

top = cbind(status.top, top.tab)
top = rbind(top.headers, top)


##### With interaction between vaccination and previous infection

out3 = coxph(Surv(tstart,tstop,status) ~ AGE + GENDER + AFFILIATION + RACE + CONDITION + CONDITIONIMMUN + MEDICATIONS + NICOTINE  + 
                 prev_positive * vaccinated, data = survival.data)
out4 = coxph(Surv(tstart,tstop,status) ~ AGE + GENDER + AFFILIATION + RACE + CONDITION + CONDITIONIMMUN + MEDICATIONS + NICOTINE  + 
                 prev_positive * manufacturer_v2, data = survival.data)

bottom.tab[1,3] = "Reference"

bottom.tab[3,3] = EstCI(out3, "vaccinated2", 1)
bottom.tab[4:6,3] = c(EstCI(out4,"manufacturer_v23", 1), 
                      EstCI(out4,"manufacturer_v24", 1), 
                      EstCI(out4,"manufacturer_v25", 1))
bottom.tab[7,3] = EstCI(out3, c("prev_positive1", 'vaccinated2', "prev_positive1:vaccinated2"), rep(1,3))
bottom.tab[8:10,3] = c(EstCI(out4, c("prev_positive1", "manufacturer_v23", "prev_positive1:manufacturer_v23"), rep(1,3)),
                       EstCI(out4, c("prev_positive1", "manufacturer_v24", "prev_positive1:manufacturer_v24"), rep(1,3)),
                       EstCI(out4, c("prev_positive1", "manufacturer_v25", "prev_positive1:manufacturer_v25"), rep(1,3)))
bottom.tab[11,3] = EstCI(out4, "prev_positive1", 1)

bottom = cbind(status.bottom, bottom.tab)
bottom = rbind(bottom.headers, bottom)


##### 2(e) Create Table 2

table2 = rbind(top, bottom)
colnames(table2) = NULL
rownames(table2) = NULL
table2