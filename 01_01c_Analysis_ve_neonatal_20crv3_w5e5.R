##############################################################################
##############################################################################
#######                                                                 ######
####### PROJECT: Temperature-related neonatal deaths attributable to    ######
#######          climate change in 29 low- and middle-income countries  ######                         
#######                                                                 ######
#######    CODE: Analysis - Very Early Neonatal Deaths [20crv3_w5e5]     ######
############################################################################## 
##############################################################################  
library(here)
library(survival)
library(dlnm)
library(tidyverse)
options(scipen=999)

here::here()


load(here("Data/Data_1st stage/data_analysis_20crv3-w5e5.RData"))

#####################################################################
##########     Model with percentiles and with lags    ##############
#####################################################################


data <- ve_neonat %>%
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,
            tmean_pct_lag2, tmean_pct_lag3, tmean_pct_lag4, tmean_pct_lag5, 
            tmean_pct_lag6,tmean_pct_lag7))  



data <- data %>% 
  dplyr::select(c(CountryName, ve_neonat, strata_id,tmean_pct_lag0, tmean_pct_lag1,
                  tmean_pct_lag2, tmean_pct_lag3, tmean_pct_lag4, tmean_pct_lag5,
                  tmean_pct_lag6,tmean_pct_lag7))



#Define a matrix of exposure histories for each observation (2 days lag)
Qneonat <- data  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag2) 

Qneonat <- as.matrix(Qneonat) 
Qneonat <- unname(Qneonat)
arglag <- list(fun="ns",knots=logknots(2,nk=1))

cb.temp <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
                knots = quantile(data$tmean_pct_lag0, c(10/100), na.rm=T),
                Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)



model <- clogit(Qneonat~cb.temp + strata(strata_id), data=data)

#Find the minimum mortality temperature
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)

cen <- cp$predvar[which.min(cp$allRRfit)] 
pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)


#Cummulative ERF with 2 DAYS LAG
jpeg(here("Figures","Very_early_neonat_mort", "Stage I", "20crv3-w5e5", "pt_Lag0-2_20crv3-w5e5.jpg"),
     units = 'in', res = 300, width = 10, height = 6)
plot(pred,"overall",  main= paste("Very early neonatal mortality 
  [dataset: 20crv3-w5e5], Lag 0-2, cen =", cen, "pt"),  
     ylab="OR",xlab="Mean temperature (percentile)",
     xlim=c(0, 100), ylim=c(0.9, 1.35), lwd=1.5)
dev.off()


pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)


#Cummulative ERF with 2 DAYS LAG (colour)
cp1cold <- crosspred(cb.temp,model,at=seq(min(data$tmean_pct_lag0,na.rm=T),
                                          cen,length=20),cen=cen)
cp1hot <- crosspred(cb.temp,model,at=seq(cen,max(data$tmean_pct_lag0,na.rm=T),
                                         length=20),cen=cen)


jpeg(here("Figures","Very_early_neonat_mort", "Stage I", "20crv3-w5e5", "ERF_20crv3-w5e5.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

plot(cp1cold,"overall",  col=4,  
     ylab="OR",xlab="Mean temperature (percentile)",
     xlim=c(0,100), ylim=c(0.9, 1.35), lwd=1.5)
lines(cp1hot,"overall",ci="area",col=2,lwd=1.5)
abline(v=cen, lty=2)

dev.off()


