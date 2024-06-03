##############################################################################
##############################################################################
#######                                                                 ######
####### PROJECT: Temperature-related neonatal deaths attributable to    ######
#######          climate change in 29 low- and middle-income countries  ######                         
#######                                                                 ######
#######    CODE: Analysis - Neonatal Deaths [gswp3-w5e5]                ######
############################################################################## 
##############################################################################  
#Definition: Number of deaths during the first 28 completed days of life 

library(here)
library(survival)
library(dlnm)
library(tidyverse)

here::here()


load(here("Data/Data_1st stage/data_analysis_gswp3-w5e5.RData"))

#####################################################################
##########     Model with percentiles and with lags    ##############
#####################################################################
library(zoo)

data <- neonat_mort  %>%
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,tmean_pct_lag2, tmean_pct_lag3, 
            tmean_pct_lag4, tmean_pct_lag5, tmean_pct_lag6, tmean_pct_lag7))   

data_neonat <- data %>%   # Number of neonatal mortality observations
  summarise(neonat_mort = sum(neonat_mort))

data <- data %>% 
  dplyr::select(c(CountryName, neonat_mort, strata_id,tmean_pct_lag0, tmean_pct_lag1,
                  tmean_pct_lag2, tmean_pct_lag3, tmean_pct_lag4, tmean_pct_lag5,
                  tmean_pct_lag6, tmean_pct_lag7)) 

################################################################################
#####            USE AIC to select best model for lag 7 days               #####
################################################################################
#Define a matrix of exposure histories for each observation (7 days lag)
Qneonat <- data  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag7) 

Qneonat <- as.matrix(Qneonat) 
Qneonat <- unname(Qneonat)
arglag <- list(fun="ns", knots=logknots(3,2), intercept=F)  

cb.temp1 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
        knots = quantile(data$tmean_pct_lag0, c(10/100), na.rm=T),  
        Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp2 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
        knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),   
        Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp3 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
        knots = quantile(data$tmean_pct_lag0, c(30/100), na.rm=T),   
        Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp4 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
        knots = quantile(data$tmean_pct_lag0, c(40/100), na.rm=T),  
        Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp5 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
       knots = quantile(data$tmean_pct_lag0, c(50/100), na.rm=T),  
       Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp6 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
       knots = quantile(data$tmean_pct_lag0, c(60/100), na.rm=T),   
       Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp7 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(70/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp8 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(80/100), na.rm=T),  
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp9 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(90/100), na.rm=T),  
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp10 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(95/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp11 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(98/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp12 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(25/100, 50/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp13 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(50/100, 75/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp14 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(50/100, 90/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp15 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(75/100, 95/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp16 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
      knots = quantile(data$tmean_pct_lag0, c(10/100, 50/100, 75/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp17 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(10/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp18 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
       knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),   
       Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp19 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
       knots = quantile(data$tmean_pct_lag0, c(30/100), na.rm=T),  
       Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp20 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
       knots = quantile(data$tmean_pct_lag0, c(40/100), na.rm=T),   
       Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp21 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
       knots = quantile(data$tmean_pct_lag0, c(50/100), na.rm=T),   
       Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp22 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
       knots = quantile(data$tmean_pct_lag0, c(60/100), na.rm=T),   
       Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp23 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
       knots = quantile(data$tmean_pct_lag0, c(70/100), na.rm=T),  
       Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp24 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(80/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp25 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(90/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp26 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(95/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp27 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(98/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp28 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(25/100, 50/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp29 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(50/100, 75/100), na.rm=T),  
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp30 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(50/100, 90/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp31 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(75/100, 95/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp32 <- crossbasis(Qneonat, lag=7, argvar=list(fun="bs", 
      knots = quantile(data$tmean_pct_lag0, c(10/100, 50/100, 75/100), na.rm=T),   
      Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)


model1 <- clogit(neonat_mort~cb.temp1 + strata(strata_id), data=data)
model2 <- clogit(neonat_mort~cb.temp2 + strata(strata_id), data=data)
model3 <- clogit(neonat_mort~cb.temp3 + strata(strata_id), data=data)
model4 <- clogit(neonat_mort~cb.temp4 + strata(strata_id), data=data)
model5 <- clogit(neonat_mort~cb.temp5 + strata(strata_id), data=data)
model6 <- clogit(neonat_mort~cb.temp6 + strata(strata_id), data=data)
model7 <- clogit(neonat_mort~cb.temp7 + strata(strata_id), data=data)
model8 <- clogit(neonat_mort~cb.temp8 + strata(strata_id), data=data)
model9 <- clogit(neonat_mort~cb.temp9 + strata(strata_id), data=data)
model10 <- clogit(neonat_mort~cb.temp10 + strata(strata_id), data=data)
model11 <- clogit(neonat_mort~cb.temp11 + strata(strata_id), data=data)
model12 <- clogit(neonat_mort~cb.temp12 + strata(strata_id), data=data)
model13 <- clogit(neonat_mort~cb.temp13 + strata(strata_id), data=data)
model14 <- clogit(neonat_mort~cb.temp14 + strata(strata_id), data=data)
model15 <- clogit(neonat_mort~cb.temp15 + strata(strata_id), data=data)
model16 <- clogit(neonat_mort~cb.temp16 + strata(strata_id), data=data)
model17 <- clogit(neonat_mort~cb.temp17 + strata(strata_id), data=data)
model18 <- clogit(neonat_mort~cb.temp18 + strata(strata_id), data=data)
model19 <- clogit(neonat_mort~cb.temp19 + strata(strata_id), data=data)
model20 <- clogit(neonat_mort~cb.temp20 + strata(strata_id), data=data)
model21 <- clogit(neonat_mort~cb.temp21 + strata(strata_id), data=data)
model22 <- clogit(neonat_mort~cb.temp22 + strata(strata_id), data=data)
model23 <- clogit(neonat_mort~cb.temp23 + strata(strata_id), data=data)
model24 <- clogit(neonat_mort~cb.temp24 + strata(strata_id), data=data)
model25 <- clogit(neonat_mort~cb.temp25 + strata(strata_id), data=data)
model26 <- clogit(neonat_mort~cb.temp26 + strata(strata_id), data=data)
model27 <- clogit(neonat_mort~cb.temp27 + strata(strata_id), data=data)
model28 <- clogit(neonat_mort~cb.temp28 + strata(strata_id), data=data)
model29 <- clogit(neonat_mort~cb.temp29 + strata(strata_id), data=data)
model30 <- clogit(neonat_mort~cb.temp30 + strata(strata_id), data=data)
model31 <- clogit(neonat_mort~cb.temp31 + strata(strata_id), data=data)
model32 <- clogit(neonat_mort~cb.temp32 + strata(strata_id), data=data)
#install.packages("AICcmodavg")
library(AICcmodavg)
models_ns <- list(model1, model2, model3, model4, model5, model6,
                  model7, model8, model9, model10, model11, model12,
                  model13, model14, model15, model16)
models_bs <- list(model17, model18,model19, model20, model21, model22, model23, 
                  model24, model25, model26, model27, model28, model29, model30,
                  model31, model32)

model.names_ns <- c('model1', 'model2', 'model3', 'model4', 'model5', 'model6',
                    'model7', 'model8', 'model9', 'model10', 'model11',
                    'model12', 'model13', 'model14', 'model15', 'model16')

model.names_bs <- c('model17','model18', 'model19', 'model20', 'model21', 'model22', 
                    'model23','model24', 'model25', 'model26', 'model27', 'model28', 
                    'model29', 'model30', 'model31', 'model32')

aictab(cand.set = models_ns, modnames = model.names_ns)
aictab(cand.set = models_bs, modnames = model.names_bs)
#Best fitting model is model 1
str(model1)

cb.temp1 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
             knots = quantile(data$tmean_pct_lag0, c(10/100), na.rm=T),   
             Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
model1 <- clogit(neonat_mort ~ cb.temp1 + strata(strata_id), data=data)

# Find the minimum mortality temperature
cenvar = 26.9
cp <- crosspred(cb.temp1, model1, cen=cenvar, by=0.1, cumul=TRUE)

cen <- cp$predvar[which.min(cp$allRRfit)] 
pred <- crosspred(cb.temp1, model1, cen=cen, by=0.1, cumul=TRUE)

str(RR)


#Sensitivity analysis
cb.temp6 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
          knots = quantile(data$tmean_pct_lag0, c(60/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
model6 <- clogit(neonat_mort ~ cb.temp1 + strata(strata_id), data=data)


cenvar = 26.9
cp <- crosspred(cb.temp1, model1, cen=cenvar, by=0.1, cumul=TRUE)

cen <- cp$predvar[which.min(cp$allRRfit)] 
pred <- crosspred(cb.temp1, model1, cen=cen, by=0.1, cumul=TRUE)

str(RR)



#################################################################
#####        USE AIC to select best model for 2 days lag    #####
#################################################################
#Define a matrix of exposure histories for each observation (7 days lag)
Qneonat <- data  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag2) 

Qneonat <- as.matrix(Qneonat) 
Qneonat <- unname(Qneonat)

arglag <- list(fun="ns",knots=logknots(2,nk=1))  #For lag 2 days

cb.temp1 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(10/100), na.rm=T),  
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp2 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),  
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp3 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(30/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp4 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(40/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp5 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(50/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp6 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(60/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp7 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(70/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp8 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(80/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp9 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(90/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp10 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(95/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp11 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(98/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp12 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(25/100, 50/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp13 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(50/100, 75/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp14 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(50/100, 90/100), na.rm=T),   
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp15 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
           knots = quantile(data$tmean_pct_lag0, c(75/100, 95/100), na.rm=T),   
           Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp16 <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
           knots = quantile(data$tmean_pct_lag0, c(10/100, 50/100, 75/100), na.rm=T),   
           Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp17 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
           knots = quantile(data$tmean_pct_lag0, c(10/100), na.rm=T), 
           Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp18 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
           knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),   
           Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp19 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(30/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp20 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(40/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp21 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(50/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp22 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(60/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp23 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(70/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp24 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(80/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp25 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(90/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp26 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(95/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp27 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(98/100), na.rm=T),  
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp28 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(25/100, 50/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp29 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(50/100, 75/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

cb.temp30 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(50/100, 90/100), na.rm=T),  
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp31 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(75/100, 95/100), na.rm=T),  
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
cb.temp32 <- crossbasis(Qneonat, lag=2, argvar=list(fun="bs", 
          knots = quantile(data$tmean_pct_lag0, c(10/100, 50/100, 75/100), na.rm=T),   
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

model1 <- clogit(neonat_mort~cb.temp1 + strata(strata_id), data=data)
model2 <- clogit(neonat_mort~cb.temp2 + strata(strata_id), data=data)
model3 <- clogit(neonat_mort~cb.temp3 + strata(strata_id), data=data)
model4 <- clogit(neonat_mort~cb.temp4 + strata(strata_id), data=data)
model5 <- clogit(neonat_mort~cb.temp5 + strata(strata_id), data=data)
model6 <- clogit(neonat_mort~cb.temp6 + strata(strata_id), data=data)
model7 <- clogit(neonat_mort~cb.temp7 + strata(strata_id), data=data)
model8 <- clogit(neonat_mort~cb.temp8 + strata(strata_id), data=data)
model9 <- clogit(neonat_mort~cb.temp9 + strata(strata_id), data=data)
model10 <- clogit(neonat_mort~cb.temp10 + strata(strata_id), data=data)
model11 <- clogit(neonat_mort~cb.temp11 + strata(strata_id), data=data)
model12 <- clogit(neonat_mort~cb.temp12 + strata(strata_id), data=data)
model13 <- clogit(neonat_mort~cb.temp13 + strata(strata_id), data=data)
model14 <- clogit(neonat_mort~cb.temp14 + strata(strata_id), data=data)
model15 <- clogit(neonat_mort~cb.temp15 + strata(strata_id), data=data)
model16 <- clogit(neonat_mort~cb.temp16 + strata(strata_id), data=data)
model17 <- clogit(neonat_mort~cb.temp17 + strata(strata_id), data=data)
model18 <- clogit(neonat_mort~cb.temp18 + strata(strata_id), data=data)
model19 <- clogit(neonat_mort~cb.temp19 + strata(strata_id), data=data)
model20 <- clogit(neonat_mort~cb.temp20 + strata(strata_id), data=data)
model21 <- clogit(neonat_mort~cb.temp21 + strata(strata_id), data=data)
model22 <- clogit(neonat_mort~cb.temp22 + strata(strata_id), data=data)
model23 <- clogit(neonat_mort~cb.temp23 + strata(strata_id), data=data)
model24 <- clogit(neonat_mort~cb.temp24 + strata(strata_id), data=data)
model25 <- clogit(neonat_mort~cb.temp25 + strata(strata_id), data=data)
model26 <- clogit(neonat_mort~cb.temp26 + strata(strata_id), data=data)
model27 <- clogit(neonat_mort~cb.temp27 + strata(strata_id), data=data)
model28 <- clogit(neonat_mort~cb.temp28 + strata(strata_id), data=data)
model29 <- clogit(neonat_mort~cb.temp29 + strata(strata_id), data=data)
model30 <- clogit(neonat_mort~cb.temp30 + strata(strata_id), data=data)
model31 <- clogit(neonat_mort~cb.temp31 + strata(strata_id), data=data)
model32 <- clogit(neonat_mort~cb.temp32 + strata(strata_id), data=data)
#install.packages("AICcmodavg")
library(AICcmodavg)
models_ns <- list(model1, model2, model3, model4, model5, model6,
                  model7, model8, model9, model10, model11, model12,
                  model13, model14, model15, model16)
models_bs <- list(model17, model18,model19, model20, model21, model22, model23, 
                  model24, model25, model26, model27, model28, model29, model30,
                  model31, model32)

model.names_ns <- c('model1', 'model2', 'model3', 'model4', 'model5', 'model6',
                    'model7', 'model8', 'model9', 'model10', 'model11',
                    'model12', 'model13', 'model14', 'model15', 'model16')

model.names_bs <- c('model17','model18', 'model19', 'model20', 'model21', 'model22', 
                    'model23','model24', 'model25', 'model26', 'model27', 'model28', 
                    'model29', 'model30', 'model31', 'model32')

aictab(cand.set = models_ns, modnames = model.names_ns)
aictab(cand.set = models_bs, modnames = model.names_bs)



cb.temp1 <- crossbasis(Qneonat, lag=7, argvar=list(fun="ns", 
           knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),   
           Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)


# # # # # Lag 2
Qneonat <- data  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag2) 

Qneonat <- as.matrix(Qneonat) 
Qneonat <- unname(Qneonat)

arglag <- list(fun="ns",knots=logknots(2,nk=1))  #For lag 2 days

cb.temp <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),  
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)
cen <- cp$predvar[which.min(cp$allRRfit)] 

pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)


jpeg(here("Figures","Neonatal mortality", "Stage I", "gswp3-w5e5", "Sensitivity analysis", "pt_Lag0-2_gswp3-w5e5_20knot.jpg"),
        units = 'in', res = 300, width = 10, height = 6)
plot(pred,"overall",  main= ("Neonatal deaths reported last 15 years, Lag 0-2"), 
           ylab="OR", xlab="Mean temperature (percentile)",
           xlim=c(0,100), ylim=c(0.9, 1.35), lwd=1.5)
 dev.off()
 ### Sensitivity analysis ---> different knot placements
 
cenvar = 26.9
 cp12 <- crosspred(cb.temp12, model12, cen=cenvar, by=0.1, cumul=TRUE)
 cen12 <- cp$predvar[which.min(cp12$allRRfit)] 
 
 pred12 <- crosspred(cb.temp12, model12, cen=cen12, by=0.1, cumul=TRUE)
 

 jpeg(here("Figures","Neonatal mortality", "Stage I", "gswp3-w5e5", "Sensitivity analysis", "knot_sensetivity_25_50knot.jpg"),
      units = 'in', res = 300, width = 10, height = 6)
 plot(pred12,"overall",  main= ("Two knots at the 25th and 50th percentiles"), 
      ylab="OR", xlab="Mean temperature (percentile)",
      xlim=c(0,100), ylim=c(0.9, 1.35), lwd=1.5)
 dev.off()
 
 
 
 ### Sensitivity analysis ---> deaths reported last 10 years
 #Model 12 is with lowest AIC
 Qneonat <- data  %>% 
   dplyr::select(tmean_pct_lag0:tmean_pct_lag2) 
 
 Qneonat <- as.matrix(Qneonat) 
 Qneonat <- unname(Qneonat)
 
 arglag <- list(fun="ns",knots=logknots(2,nk=1))  #For lag 2 days
 cb.temp <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
                knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),   
                Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
 
 model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)
 cenvar = 26.9
 cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)
 cen <- cp$predvar[which.min(cp$allRRfit)] 
 
 pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)
 
 jpeg(here("Figures","Neonatal mortality", "Stage I", "gswp3-w5e5", 
           "Sensitivity analysis", "pt_Lag0-2_gswp3-w5e5_10_years_20knot.jpg"),
            units = 'in', res = 300, width = 10, height = 6)
 plot(pred,"overall", main= ("Neonatal deaths  reported last 10 years, Lag 0-2"), 
      ylab="OR", xlab="Mean temperature (percentile)",
      xlim=c(0,100), ylim=c(0.9, 1.35), lwd=1.5)
 dev.off()
 
 ### Sensitivity analysis ---> deaths reported last 5 years
 #Model 3 is with lowest AIC
 Qneonat <- data  %>% 
   dplyr::select(tmean_pct_lag0:tmean_pct_lag2) 
 
 Qneonat <- as.matrix(Qneonat) 
 Qneonat <- unname(Qneonat)
 
 arglag <- list(fun="ns",knots=logknots(2,nk=1))  #For lag 2 days
 
 cb.temp <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
          knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),  
          Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
 model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)
 cenvar = 26.9
 cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)
 cen <- cp$predvar[which.min(cp$allRRfit)] 
 
 pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)
 
 jpeg(here("Figures","Neonatal mortality", "Stage I", "gswp3-w5e5", 
           "Sensitivity analysis", "pt_Lag0-2_gswp3-w5e5_5_years_20knot.jpg"),
      units = 'in', res = 300, width = 10, height = 6)
 plot(pred,"overall", main= ("Neonatal deaths reported last 5 years, Lag 0-2"), 
      ylab="OR", xlab="Mean temperature (percentile)",
      xlim=c(0,100), ylim=c(0.9, 1.35), lwd=1.5)
 dev.off()

######## For lag 3 days - sensitivity analysis as results for lag 3 fall shortly of significance


Qneonat <- data  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag3) 
Qneonat <- as.matrix(Qneonat) 
Qneonat <- unname(Qneonat)
arglag <- list(fun="ns", knots=logknots(3,2), intercept=F)  #For lag 3 days
cb.temp <- crossbasis(Qneonat, lag=3, argvar=list(fun="ns", 
             knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),  
             Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)


model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)
cen <- cp$predvar[which.min(cp$allRRfit)] 

pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)


#Cummulative ERF 3 DAYS LAG
 jpeg(here("Figures","Neonatal mortality", "Stage I", "gswp3-w5e5", "pt_Lag0-3_gswp3-w5e5.jpg"),
      units = 'in', res = 300, width = 10, height = 6)
      plot(pred,"overall",  main= ("Neonatal mortality, Lag 0-3"), 
     ylab="OR", xlab="Mean temperature (percentile)",
     xlim=c(0,100), ylim=c(0.9, 1.35), lwd=1.5)
 dev.off()







