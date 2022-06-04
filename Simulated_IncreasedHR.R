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
time_since_infect = rep(0,n)
time_since_infect[prev_positive == 1] = round(rnorm(sum(prev_positive == 1), 270, 70))
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

covariates = data.frame(ID = 1:5000, AGE = AGE, GENDER = GENDER, AFFILIATION = AFFILIATION, RACE = race_eth, CONDITION = CONDITION,
                        CONDITIONIMMUN = CONDITIONIMMUN, MEDICATIONS = MEDICATIONS, NICOTINE = NICOTINE, 
                        prev_positive = prev_positive, manufacturer = manufacturer.2, time_since_infect = time_since_infect)
covariates$manufacturer = ifelse(is.na(covariates$manufacturer), 0, covariates$manufacturer)
covariates2 = data.frame(covariates, vaccinated = vaccinated, status = status, survival.time = time,
                         first.dose = first.dose, second.dose = second.dose)


##### 1(c) Create data where time since infection & time since vaccination are time-varying covariates

# 0. Set up the survival data matrix

survival.data = matrix(0,ncol=9,nrow=0)
colnames(survival.data) = c('ID', 'tstart', 'tstop', 'status', 'vaccinated', 
                            'time_since_p_vaccination', 'time_since_vaccination', 
                            'prev_positive', 'time_since_infect')


# 1. No vaccination (OR first.dose >= survival.time); no previous infection

temp = covariates2 %>% filter((manufacturer == 0 & prev_positive == 0) | (survival.time <= first.dose & prev_positive == 0))

for(i in 1:nrow(temp)){
    record = c(temp[i,'ID'], 0, temp[i,'survival.time'], temp[i,'status'], rep(0,5))
    survival.data = rbind(survival.data, record)
}


# 2. No vaccination (OR first.dose >= survival.time); previous infection
#
# colnames(survival.data) = c('STUDY_ID', 'tstart', 'tstop', 'status', 'vaccinated', 
#                             'time_since_p_vaccination', 'time_since_vaccination', 
#                             'prev_positive', 'time_since_infect')

temp = covariates2 %>% filter((manufacturer == 0 & prev_positive == 1) | (survival.time <= first.dose & prev_positive == 1))

for(i in 1:nrow(temp)){
    if(temp[i,'survival.time']==0){
        record = c(temp[i,'ID'], 0, 0.5, temp[i,'status'], 0, 0, 0, 1, temp[i,'time_since_infect'])
        survival.data = rbind(survival.data, record)
    } else{
        record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'])
        record[,1] = temp[i,'ID']
        record[,2] = 0:(temp[i,'survival.time']-1)
        record[,3] = 1:(temp[i,'survival.time'])
        record[nrow(record),4] = temp[i,'status']
        record[,8] = 1
        record[,9] = (temp[i,'time_since_infect'] + 1):(temp[i,'time_since_infect'] + temp[i,'survival.time'])
        survival.data = rbind(survival.data, record)
    }
}

manufacturer.0 = rep(0, nrow(survival.data))

# 3. J&J, first.dose <= follow-up start
#
# colnames(survival.data) = c('STUDY_ID', 'tstart', 'tstop', 'status', 'vaccinated', 
#                             'time_since_p_vaccination', 'time_since_vaccination', 
#                             'prev_positive', 'time_since_infect')

temp = covariates2 %>% filter(manufacturer == 3,first.dose <=0)

for(i in 1:nrow(temp)){
    if(temp[i,'prev_positive']==0){
        if(temp[i,'survival.time']==0){
            record = c(temp[i,'ID'], 0, 0.5, temp[i,'status'], 2, 0, temp[i,'survival.time']-temp[i,'first.dose'], 0, 0)
            survival.data = rbind(survival.data, record)
        } else{
            record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'])
            record[,1] = temp[i,'ID']
            record[,2] = 0:(temp[i,'survival.time']-1)
            record[,3] = 1:(temp[i,'survival.time'])
            record[nrow(record),4] = temp[i,'status']
            record[,5] = 2
            record[,6] = 0
            record[,7] = record[,2] - temp[i,'first.dose']
            survival.data = rbind(survival.data, record)
        }
    } else{
        if(temp[i,'survival.time']==0){
            record = c(temp[i,'ID'], 0, 0.5, temp[i,'status'], 2, 0, temp[i,'survival.time']-temp[i,'first.dose'], 
                       1, temp[i,'time_since_infect'])
            survival.data = rbind(survival.data, record)
        } else{
            record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'])
            record[,1] = temp[i,'ID']
            record[,2] = 0:(temp[i,'survival.time']-1)
            record[,3] = 1:(temp[i,'survival.time'])
            record[nrow(record),4] = temp[i,'status']
            record[,5] = 2
            record[,6] = 0
            record[,7] = record[,3] - temp[i,'first.dose']
            record[,8] = temp[i,'prev_positive']
            record[,9] = record[,3] + temp[i,'time_since_infect']
            survival.data = rbind(survival.data, record)
        }
    }
}

# 4. J&J, follow-up start < first.dose < survival.time
#
# colnames(survival.data) = c('ID', 'tstart', 'tstop', 'status', 'vaccinated', 
#                             'time_since_p_vaccination', 'time_since_vaccination', 
#                             'prev_positive', 'time_since_infect')

temp = covariates2 %>% filter(manufacturer == 3, first.dose > 0, survival.time > first.dose)


for(i in 1:nrow(temp)){
    if(temp[i,'prev_positive']==0){
        record = matrix(0, ncol=9, nrow = temp[i,'survival.time'] - temp[i,'first.dose'] + 1)
        record[,1] = temp[i,'ID']
        record[,3] = temp[i,'first.dose']:temp[i,'survival.time']
        record[,2] = record[,3] - 1; record[1,2] = 0
        record[nrow(record),4] = temp[i,'status']
        record[,5] = ifelse(record[,3] > temp[i,'first.dose'],2,0)
        record[record[,5]==2,7] = 1:sum(record[,5]==2)
        survival.data = rbind(survival.data, record)
    } else{
        record = matrix(0, ncol=9, nrow = temp[i,'survival.time'] - temp[i,'first.dose'] + 1)
        record[,1] = temp[i,'ID']
        record[,3] = temp[i,'first.dose']:temp[i,'survival.time']
        record[,2] = record[,3] - 1; record[1,2] = 0
        record[nrow(record),4] = temp[i,'status']
        record[,5] = ifelse(record[,3] > temp[i,'first.dose'],2,0)
        record[record[,5]==2,7] = 1:sum(record[,5]==2)
        record[,8] = temp[i,'prev_positive']
        record[,9] = record[,3] + temp[i,'time_since_infect']
        survival.data = rbind(survival.data, record)
    }
}

manufacturer.3 = rep(3, nrow(survival.data)-length(manufacturer.0))


# 5. mRNA, second.dose <= follow-up start
#
# colnames(survival.data) = c('ID', 'tstart', 'tstop', 'status', 'vaccinated', 
#                             'time_since_p_vaccination', 'time_since_vaccination', 
#                             'prev_positive', 'time_since_infect')

temp = covariates2 %>% filter(manufacturer %in% c(1,2), second.dose <= 0)

for(i in 1:nrow(temp)){
    if(temp[i,'prev_positive']==0){
        if(temp[i,'survival.time']==0){
            record = c(temp[i,'ID'], 0, 0.5, temp[i,'status'], 2, 0, temp[i,'survival.time']-temp[i,'second.dose'], 0, 0)
            survival.data = rbind(survival.data, record)
        } else{
            record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'])
            record[,1] = temp[i,'ID']
            record[,2] = 0:(temp[i,'survival.time']-1)
            record[,3] = 1:(temp[i,'survival.time'])
            record[nrow(record),4] = temp[i,'status']
            record[,5] = 2
            record[,6] = 0
            record[,7] = record[,2] - temp[i,'second.dose']
            survival.data = rbind(survival.data, record)
        }
    } else{
        if(temp[i,'survival.time']==0){
            record = c(temp[i,'ID'], 0, 0.5, temp[i,'status'], 2, 0, temp[i,'survival.time']-temp[i,'second.dose'], 
                       1, temp[i,'time_since_infect'])
            survival.data = rbind(survival.data, record)
        } else{
            record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'])
            record[,1] = temp[i,'ID']
            record[,2] = 0:(temp[i,'survival.time']-1)
            record[,3] = 1:(temp[i,'survival.time'])
            record[nrow(record),4] = temp[i,'status']
            record[,5] = 2
            record[,6] = 0
            record[,7] = record[,3] - temp[i,'second.dose']
            record[,8] = temp[i,'prev_positive']
            record[,9] = record[,3] + temp[i,'time_since_infect']
            survival.data = rbind(survival.data, record)
        }
    }
}

# 6. mRNA, first.dose <= follow-up start < second.dose <= survival.time
#
# colnames(survival.data) = c('ID', 'tstart', 'tstop', 'status', 'vaccinated', 
#                             'time_since_p_vaccination', 'time_since_vaccination', 
#                             'prev_positive', 'time_since_infect')

temp = covariates2 %>% filter(manufacturer %in% c(1,2), first.dose <= 0, second.dose > 0, second.dose <= survival.time)

for(i in 1:nrow(temp)){
    if(temp[i,'prev_positive']==0){
        record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'])
        record[,1] = temp[i,'ID']
        record[,2] = 0:(temp[i,'survival.time']-1)
        record[,3] = 1:(temp[i,'survival.time'])
        record[nrow(record),4] = temp[i,'status']
        record[,5] = 1 + (record[,3] > temp[i,'second.dose'])
        record[record[,5]==1,6] = record[record[,5]==1,3] - temp[i,'first.dose']
        record[record[,5]==2,7] = 1:sum(record[,5]==2)
        survival.data = rbind(survival.data, record)
    } else{
        record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'])
        record[,1] = temp[i,'ID']
        record[,2] = 0:(temp[i,'survival.time']-1)
        record[,3] = 1:(temp[i,'survival.time'])
        record[nrow(record),4] = temp[i,'status']
        record[,5] = 1 + (record[,3] > temp[i,'second.dose'])
        record[record[,5]==1,6] = record[record[,5]==1,3] - temp[i,'first.dose']
        record[record[,5]==2,7] = 1:sum(record[,5]==2)
        record[,8] = temp[i,'prev_positive']
        record[,9] = record[,3] + temp[i,'time_since_infect']
        survival.data = rbind(survival.data, record)
    }
}

# 7. mRNA, first.dose <= follow-up start <= survival.time < second.dose (OR no second dose)
#
# colnames(survival.data) = c('ID', 'tstart', 'tstop', 'status', 'vaccinated', 
#                             'time_since_p_vaccination', 'time_since_vaccination', 
#                             'prev_positive', 'time_since_infect')


temp = covariates2 %>% filter(manufacturer %in% c(1,2), first.dose <= 0, 
                           (second.dose > 0 & second.dose > survival.time) | is.na(second.dose))

for(i in 1:nrow(temp)){
    if(temp[i,'prev_positive']==0){
        if(temp[i,'survival.time']==0){
            record = c(temp[i,'ID'], 0, 0.5, temp[i,'status'], 2, 0, temp[i,'survival.time']-temp[i,'second.dose'], 0, 0)
            survival.data = rbind(survival.data, record)
        } else{
            record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'])
            record[,1] = temp[i,'ID']
            record[,2] = 0:(temp[i,'survival.time']-1)
            record[,3] = 1:(temp[i,'survival.time'])
            record[nrow(record),4] = temp[i,'status']
            record[,5] = 1
            record[,6] = record[,3] - temp[i,'first.dose']
            survival.data = rbind(survival.data, record)
        }
    } else{
        if(temp[i,'survival.time']==0){
            record = c(temp[i,'ID'], 0, 0.5, temp[i,'status'], 2, 0, temp[i,'survival.time']-temp[i,'second.dose'], 
                       1, temp[i,'time_since_infect'])
            survival.data = rbind(survival.data, record)
        } else{
            record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'])
            record[,1] = temp[i,'ID']
            record[,2] = 0:(temp[i,'survival.time']-1)
            record[,3] = 1:(temp[i,'survival.time'])
            record[nrow(record),4] = temp[i,'status']
            record[,5] = 1
            record[,6] = record[,3] - temp[i,'first.dose']
            record[,8] = temp[i,'prev_positive']
            record[,9] = record[,3] + temp[i,'time_since_infect']
            survival.data = rbind(survival.data, record)
        }
    }
}

# 8. mRNA, 0 < first.dose < survival.time, second.dose <= survival.time

temp = covariates2 %>% filter(manufacturer %in% c(1,2), first.dose > 0, first.dose < survival.time, 
                           second.dose <= survival.time)


for(i in 1:nrow(temp)){
    if(temp[i,'prev_positive']==0){
        record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'] - temp[i,'first.dose'] + 1)
        record[,1] = temp[i,'ID']
        record[,3] = temp[i,'first.dose']:temp[i,'survival.time']
        record[,2] = record[,3] - 1; record[1,2] = 0
        record[nrow(record),4] = temp[i,'status']
        record[,5] = (record[,3]>temp[i,'first.dose']) + (record[,3] > temp[i,'second.dose'])
        record[record[,5]==1,6] = 1:sum(record[,5]==1)
        record[record[,5]==2,7] = 1:sum(record[,5]==2)
        survival.data = rbind(survival.data, record)
    } else{
        record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'] - temp[i,'first.dose'] + 1)
        record[,1] = temp[i,'ID']
        record[,3] = temp[i,'first.dose']:temp[i,'survival.time']
        record[,2] = record[,3] - 1; record[1,2] = 0
        record[nrow(record),4] = temp[i,'status']
        record[,5] = (record[,3]>temp[i,'first.dose']) + (record[,3] > temp[i,'second.dose'])
        record[record[,5]==1,6] = 1:sum(record[,5]==1)
        record[record[,5]==2,7] = 1:sum(record[,5]==2)
        record[,8] = temp[i,'prev_positive']
        record[,9] = record[,3] + temp[i,'time_since_infect']
        survival.data = rbind(survival.data, record)
    }
}

# 9. mRNA, 0 < first.dose < survival.time, second.dose > survival.time (OR no second dose)

temp = covariates2 %>% filter(manufacturer %in% c(1,2), first.dose > 0, first.dose < survival.time, 
                           second.dose > survival.time | is.na(second.dose))


for(i in 1:nrow(temp)){
    if(temp[i,'prev_positive']==0){
        record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'] - temp[i,'first.dose'] + 1)
        record[,1] = temp[i,'ID']
        record[,3] = temp[i,'first.dose']:temp[i,'survival.time']
        record[,2] = record[,3] - 1; record[1,2] = 0
        record[nrow(record),4] = temp[i,'status']
        record[,5] = (record[,3]>temp[i,'first.dose'])
        record[record[,5]==1,6] = 1:sum(record[,5]==1)
        survival.data = rbind(survival.data, record)
    } else{
        record = matrix(0, ncol = 9, nrow = temp[i,'survival.time'] - temp[i,'first.dose'] + 1)
        record[,1] = temp[i,'ID']
        record[,3] = temp[i,'first.dose']:temp[i,'survival.time']
        record[,2] = record[,3] - 1; record[1,2] = 0
        record[nrow(record),4] = temp[i,'status']
        record[,5] = (record[,3]>temp[i,'first.dose'])
        record[record[,5]==1,6] = 1:sum(record[,5]==1)
        record[,8] = temp[i,'prev_positive']
        record[,9] = record[,3] + temp[i,'time_since_infect']
        survival.data = rbind(survival.data, record)
    }
}

manufacturer.1 = rep(1,nrow(survival.data) - length(manufacturer.0) - length(manufacturer.3))
survival.data = cbind(survival.data,c(manufacturer.0, manufacturer.3, manufacturer.1))


survival.data = data.frame(survival.data)
rownames(survival.data) = NULL
colnames(survival.data)[10] = "manufacturer"

survival.data = left_join(survival.data, select(covariates2, ID, AGE, GENDER, AFFILIATION, RACE, CONDITION, CONDITIONIMMUN, 
                                                MEDICATIONS, NICOTINE, manufacturer), by = 'ID')
survival.data$tstop = ifelse(survival.data$tstop == 0, 0.5, survival.data$tstop)
survival.data$manufacturer.x = ifelse(survival.data$manufacturer.x == 1 & survival.data$manufacturer.y == 2, 2, survival.data$manufacturer.x)

survival.data = survival.data %>% mutate(manufacturer_v2 = case_when(vaccinated == 0 ~ 0,
                                                                     vaccinated == 1 & manufacturer.x == 1 ~ 1,
                                                                     vaccinated == 1 & manufacturer.x == 2 ~ 2,
                                                                     vaccinated == 2 & manufacturer.x == 1 ~ 3,
                                                                     vaccinated == 2 & manufacturer.x == 2 ~ 4,
                                                                     vaccinated == 2 & manufacturer.x == 3 ~ 5))


survival.data = survival.data %>% mutate(t_partial_moderna = ifelse(manufacturer_v2 == 1, time_since_p_vaccination, 0),
                                         t_partial_pfizer = ifelse(manufacturer_v2 == 2, time_since_p_vaccination, 0),
                                         t_moderna = ifelse(manufacturer_v2 == 3, time_since_vaccination, 0),
                                         t_pfizer = ifelse(manufacturer_v2 == 4, time_since_vaccination, 0),
                                         t_jandj = ifelse(manufacturer_v2 == 5, time_since_vaccination, 0))

survival.data = survival.data %>% mutate_at(vars(c('GENDER','AFFILIATION','CONDITION','CONDITIONIMMUN',
                                                   'MEDICATIONS','NICOTINE','RACE',
                                                   'prev_positive','vaccinated','manufacturer.x','manufacturer_v2')), ~as.factor(.))

##### 2. Fit time-varying Cox model and estimate the monthly decline in HR

##### 2(a). Functions to estimate protection & decline in HR

EstCI.Protection = function(output, terms, new.x){
    coef = output$coefficients[terms]
    var = vcov(output)[terms,terms]
    
    point = format(round(1 - exp(new.x%*%coef),3)*100,nsmall=1, trim=T)
    SE = sqrt(new.x%*%var%*%new.x)
    
    CI = format(round(c(1 - exp(new.x%*%coef + 1.96*SE),  1 - exp(new.x%*%coef - 1.96*SE)),3)*100,nsmall = 1,trim=T)
    est = paste0(point, '% (', CI[1],'-',CI[2],')')
    return(est)
}

EstCI.Decline = function(output, terms, new.x){
    coef = output$coefficients[terms]
    var = vcov(output)[terms,terms]
    
    point = format(round(exp(new.x%*%coef),2),nsmall=2, trim=T)
    SE = sqrt(new.x %*% var %*% new.x)
    
    CI = format(round(c(exp(new.x%*%coef - 1.96*SE),  exp(new.x%*%coef + 1.96*SE)),2),nsmall = 2,trim=T)
    est = paste0(point, ' (', CI[1],'-',CI[2],')')
    return(est)
}


##### 2(b) Fit time-varying Cox model and estimate:

#####      (1) Protection at 3-months after vaccination/previous infection
#####      (2) Protection at 6-months after vaccination/previous infection
#####      (3) Monthly increase in hazard ratio (declined protection)

headers = c("Status", "3-Month Protection: % (95% CI)", "6-Month Protection: % (95% CI)", "Increase in Monthly HR")
status = c("All Vaccines", "...mRNA-1273", "...BNT162b2", "...Ad26.COV2.S", "Previous Infection")

table3 = matrix(NA, 5, 3)


# all vaccines
out = coxph(Surv(tstart, tstop, status) ~ AGE + GENDER + AFFILIATION + RACE + CONDITION + CONDITIONIMMUN + MEDICATIONS + 
                NICOTINE + vaccinated + time_since_p_vaccination + time_since_vaccination + 
                prev_positive + time_since_infect, 
            data = survival.data)

# by vaccines 
out2 = coxph(Surv(tstart, tstop, status) ~ AGE + GENDER + AFFILIATION + RACE + CONDITION + CONDITIONIMMUN + MEDICATIONS + 
                 NICOTINE + manufacturer_v2 + t_partial_moderna + t_partial_pfizer + t_moderna + t_pfizer + t_jandj + 
                 prev_positive + time_since_infect, 
             data = survival.data)

table3[1:4,1] = c(EstCI.Protection(out, c("vaccinated2", "time_since_vaccination"), c(1,90)), 
               EstCI.Protection(out2, c("manufacturer_v23", "t_moderna"), c(1,90)), 
               EstCI.Protection(out2, c("manufacturer_v24", "t_pfizer"), c(1,90)),
               EstCI.Protection(out2, c("manufacturer_v25", "t_jandj"), c(1,90)))

table3[1:4,2] = c(EstCI.Protection(out, c("vaccinated2", "time_since_vaccination"), c(1,180)), 
               EstCI.Protection(out2, c("manufacturer_v23", "t_moderna"), c(1,180)), 
               EstCI.Protection(out2, c("manufacturer_v24", "t_pfizer"), c(1,180)),
               EstCI.Protection(out2, c("manufacturer_v25", "t_jandj"), c(1,180)))

table3[1:4,3] = c(EstCI.Decline(out, "time_since_vaccination", 30),
               EstCI.Decline(out2, "t_moderna", 30),
               EstCI.Decline(out2, "t_pfizer", 30),
               EstCI.Decline(out2, "t_jandj", 30))

table3[5,] = c(EstCI.Protection(out2, c("prev_positive1","time_since_infect"), c(1, 90)), 
            EstCI.Protection(out2, c("prev_positive1","time_since_infect"), c(1, 180)), 
            EstCI.Decline(out2, c("time_since_infect"), 30))


table3 = cbind(status, table3)
table3 = rbind(headers, table3)
colnames(table3) = NULL
rownames(table3) = NULL

table3
