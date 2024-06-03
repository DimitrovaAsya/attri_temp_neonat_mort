##############################################################################
##############################################################################
#######                                                                 ######
####### PROJECT: Temperature-related neonatal deaths attributable to    ######
#######          climate change in 29 low- and middle-income countries  ######                         
#######                                                                 ######
#######    CODE: Analysis - Neonatal Deaths [20crv3-era5]              ######
############################################################################## 
##############################################################################  
#Definition: Number of deaths during the first 28 completed days of life 

library(here)
library(survival)
library(dlnm)
library(tidyverse)

here::here()

load(here("Data/Data_1st stage/data_analysis_20crv3-era5.RData"))


#####################################################################
##########     Model with percentiles and with lags    ##############
#####################################################################

data <- neonat_mort  %>%
  mutate(year  = format(date, format = "%Y"),
         year=as.numeric(year),
         month = format(date, format = "%m"),
         month=as.numeric(month)) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,
            tmean_pct_lag2, tmean_pct_lag3, tmean_pct_lag4, tmean_pct_lag5, 
            tmean_pct_lag6,tmean_pct_lag7))  



data <- data %>% 
  dplyr::select(c(CountryName, neonat_mort, strata_id,tmean_pct_lag0, tmean_pct_lag1,
                  tmean_pct_lag2, tmean_pct_lag3, tmean_pct_lag4, tmean_pct_lag5,
                  tmean_pct_lag6, tmean_pct_lag7)) 


#Define a matrix of exposure histories for each observation (2 days lag)
Qneonat <- data  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag2) 

Qneonat <- as.matrix(Qneonat) 
Qneonat <- unname(Qneonat)
arglag <- list(fun="ns",knots=logknots(2,nk=1))  #For lag 2 days

cb.temp <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
                              knots = quantile(data$tmean_pct_lag0, c(20/100), na.rm=T),  
                              Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)
#Best fitting is model 2
model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)

# Find the minimum mortality tmp 
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)


cen <- cp$predvar[which.min(cp$allRRfit)] 

pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)


#Cumulative ERF - 2 DAYS LAG
jpeg(here("Figures","Neonatal mortality", "Stage I", "20crv3-era5", "pt_Lag0-2_20crv3-era5.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

plot(pred,"overall",  main= paste("Neonatal mortality 
[dataset: 20crv3-era5], Lag 0-2, cen =", cen, "pt"),  
     ylab="OR",xlab="Mean temperature (percentile)",
     xlim=c(0,100), ylim=c(0.9, 1.35), lwd=1.5)
dev.off()


# 2 DAYS LAGS COLOUR
cp1cold <- crosspred(cb.temp,model,at=seq(min(data$tmean_pct_lag0,na.rm=T),
                                          cen,length=20),cen=cen)
cp1hot <- crosspred(cb.temp,model,at=seq(cen,max(data$tmean_pct_lag0,na.rm=T),
                                         length=20),cen=cen)


jpeg(here("Figures","Neonatal mortality", "Stage I", "20crv3-era5", "ERF_2daylag_20crv3-era5.jpg"),
     units = 'in', res = 300, width = 10, height = 6)
plot(cp1cold,"overall",  col=4,  
     ylab="OR",xlab="Mean temperature (percentile)",
     xlim=c(0,100), ylim=c(0.9, 1.35), lwd=1.5)
lines(cp1hot,"overall",ci="area",col=2,lwd=1.5)
abline(v=cen, lty=2)
# title("Overall cumulative temperature-neonatal mortality association")
dev.off()


