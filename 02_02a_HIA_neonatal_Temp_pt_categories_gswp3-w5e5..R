#################################################################################
#################################################################################
#######                                                                    ######
####### PROJECT: Temperature-related neonatal deaths attributable to       ######
#######          climate change in 29 low- and middle-income countries     ######                         
#######                                                                    ######
####### CODE: HIA [CATEGORIES] -  Early Neonatal Deaths; [gswp3-w5e5]      ######
################################################################################# 
#################################################################################  


library(reshape2)
library(here)
library(survival)
library(dlnm)
library(tidyverse)
options(scipen=999)
source(here("Codes/attrdl.R"))

here::here()

load(here("Data/Data_1st stage/data_analysis_gswp3-w5e5.RData"))

########################################################
######           0 FACTUAL SCENARIO               ###### 
########################################################

########################################################
######           0-1  EXTREMELY HOT               ###### 
########################################################
data <- neonat_mort  %>%
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,tmean_pct_lag2, tmean_pct_lag3, 
            tmean_pct_lag4, tmean_pct_lag5, tmean_pct_lag6, tmean_pct_lag7))   

data_neonat <- data %>%   # Number of neonatal deaths
  summarise(neonat_mort = sum(neonat_mort))

data <- data %>% 
  dplyr::select(c(CountryName, year, neonat_mort, strata_id,tmean_pct_lag0, tmean_pct_lag1,
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


model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)

# Find the minimum mortality temperature 
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)


cen <- cp$predvar[which.min(cp$allRRfit)] 

pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)


countries <- unique(data$CountryName)

df_0 = NULL
df_1 = NULL


for(j in countries){
  country <- j
  df_2 <- data %>%
    dplyr::filter(CountryName == j)   
  
  if (nrow(df_2) == 0) next
  Qneonat <- df_2  %>% 
    dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
  Qneonat <- as.matrix(Qneonat) 
  Qneonat <- unname(Qneonat)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
  #95% CIs calculated
  
  attr_f_hot <- attrdl(Qneonat, cb.temp,df_2$neonat_mort, model, type="af", 
                    cen=cen, range=c(97.5,100))   
  attr_f_hot <- as.data.frame(attr_f_hot)
  set.seed(42)
  CIs_hot <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
                             type="af", cen=cen, range=c(97.5,100), sim=TRUE, 
                             nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)  
  CIs_hot <- as.data.frame(t(CIs_hot))
  CIs_hot <- CIs_hot %>% 
    dplyr::rename('attr_f_hl' = '2.5%',
                  'attr_f_hh' = '97.5%') 
  
  df_attr <- cbind (attr_f_hot, CIs_hot) %>%
    mutate(country = j)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
  #95% CIs calculated
  set.seed(42)
  sim_hot <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                    type="af",cen=cen, range=c(97.5,100), sim=TRUE, nsim=10000)  
  sim_hot <- as.data.frame(sim_hot)
  sim_hot$country <- j
  sim_hot$scenario <- "Factual"
  sim_hot$dataset <- "gswp3-w5e5"
  
  sim <- sim_hot 
  
  df_0 <- bind_rows(df_0, df_attr)
  df_1 <- bind_rows(df_1, sim)
  
}

attr_frac_df <- df_0 %>% 
  mutate(scen="obsclim",
         sim="gswp3-w5e5",
         range ="extremely hot") %>% 
  write.csv(file= here("Results", "Attributable_neonatal", "gswp3-w5e5", 
          "obsclim_gswp3-w5e5_extr_hot.csv"), row.names=FALSE)


df_1 <- df_1 %>% 
  mutate(range ="extremely hot") %>% 
  write.csv(file= here("Results", "Attributable_neonatal", "Monte_Carlo", 
                       "Factual_gswp3-w5e5_extr_hot.csv"))

########################################################
######           0-2  MODERATELY  HOT             ###### 
########################################################
library(reshape2)
library(here)
library(survival)
library(dlnm)
library(tidyverse)
options(scipen=999)
source(here("Codes/attrdl.R"))

here::here()

load(here("Data/Data_1st stage/data_analysis_gswp3-w5e5.RData"))

data <- neonat_mort  %>%
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,tmean_pct_lag2, tmean_pct_lag3, 
            tmean_pct_lag4, tmean_pct_lag5, tmean_pct_lag6, tmean_pct_lag7))   

data_neonat <- data %>%   
  summarise(neonat_mort = sum(neonat_mort))

data <- data %>% 
  dplyr::select(c(CountryName, year, neonat_mort, strata_id,tmean_pct_lag0, tmean_pct_lag1,
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


model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)

# Find the minimum mortality temperature 
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)


cen <- cp$predvar[which.min(cp$allRRfit)] 

pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)


countries <- unique(data$CountryName)

df_0 = NULL
df_1 = NULL



for(j in countries){
  country <- j
  df_2 <- data %>%
    dplyr::filter(CountryName == j)
  
  if (nrow(df_2) == 0) next
  Qneonat <- df_2  %>% 
    dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
  Qneonat <- as.matrix(Qneonat) 
  Qneonat <- unname(Qneonat)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
  #95% CIs calculated
  attr_f_hot <- attrdl(Qneonat, cb.temp,df_2$neonat_mort, model, type="af", cen=cen,
                       range=c(75, 97.5))   
  attr_f_hot <- as.data.frame(attr_f_hot)
  set.seed(42)
  CIs_hot <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
                             type="af", cen=cen, range=c(75, 97.5), 
                             sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE) 
  CIs_hot <- as.data.frame(t(CIs_hot))
  CIs_hot <- CIs_hot %>% 
    dplyr::rename('attr_f_hl' = '2.5%',
                  'attr_f_hh' = '97.5%') 
  
  df_attr <- cbind (attr_f_hot, CIs_hot) %>%
    mutate(country = j)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
  #95% CIs calculated
  set.seed(42)
  sim_hot <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                 type="af",cen=cen, range=c(75, 97.5), sim=TRUE, nsim=10000)  
  sim_hot <- as.data.frame(sim_hot)
  sim_hot$country <- j
  sim_hot$scenario <- "Factual"
  sim_hot$dataset <- "gswp3-w5e5"
  sim <- sim_hot 
  
  df_0 <- bind_rows(df_0, df_attr)
  df_1 <- bind_rows(df_1, sim)
  
}


attr_frac_df <- df_0 %>% 
  mutate(scen="obsclim",
         sim="gswp3-w5e5",
         range ="moderately hot") %>% 
  write.csv(file= here("Results", "Attributable_neonatal", "gswp3-w5e5", 
          "obsclim_gswp3-w5e5_mod_hot.csv"), row.names=FALSE)

df_1 <- df_1 %>% 
  mutate(range ="moderately hot") %>% 
  write.csv(file= here("Results", "Attributable_neonatal", "Monte_Carlo", 
                       "Factual_gswp3-w5e5_mod_hot.csv"))

########################################################
######           0-3  MILDLY  HOT                 ###### 
########################################################
library(reshape2)
library(here)
library(survival)
library(dlnm)
library(tidyverse)
options(scipen=999)
source(here("Codes/attrdl.R"))

here::here()

load(here("Data/Data_1st stage/data_analysis_gswp3-w5e5.RData"))

here::here()

data <- neonat_mort  %>%
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,tmean_pct_lag2, tmean_pct_lag3, 
            tmean_pct_lag4, tmean_pct_lag5, tmean_pct_lag6, tmean_pct_lag7))   

data_neonat <- data %>%   # Number of neonat observations
  summarise(neonat_mort = sum(neonat_mort))

data <- data %>% 
  dplyr::select(c(CountryName, year, neonat_mort, strata_id,tmean_pct_lag0, tmean_pct_lag1,
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


model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)

#Find the minimum mortality temperature 
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)


cen <- cp$predvar[which.min(cp$allRRfit)] 

pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)


countries <- unique(data$CountryName)
df_0 = NULL
df_1 = NULL


for(j in countries){
  country <- j
  df_2 <- data %>%
    dplyr::filter(CountryName == j)
  
  if (nrow(df_2) == 0) next
  Qneonat <- df_2  %>% 
    dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
  Qneonat <- as.matrix(Qneonat) 
  Qneonat <- unname(Qneonat)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
  #95% CIs calculated
  attr_f_hot <- attrdl(Qneonat, cb.temp,df_2$neonat_mort, model, type="af", cen=cen,
                       range=c(cen ,75))  
  attr_f_hot <- as.data.frame(attr_f_hot)
  set.seed(42)
  CIs_hot <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
                             type="af", cen=cen, range=c(cen ,75), sim=TRUE, 
                             nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)  
  CIs_hot <- as.data.frame(t(CIs_hot))
  CIs_hot <- CIs_hot %>% 
    dplyr::rename('attr_f_hl' = '2.5%',
                  'attr_f_hh' = '97.5%') 
  
  df_attr <- cbind (attr_f_hot, CIs_hot) %>%
    mutate(country = j)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
  #95% CIs calculated
  set.seed(42)
  sim_hot <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                 type="af",cen=cen, range=c(cen ,75), sim=TRUE, nsim=10000)  
  sim_hot <- as.data.frame(sim_hot)
  sim_hot$country <- j
  sim_hot$scenario <- "Factual"
  sim_hot$dataset <- "gswp3-w5e5"
  
  sim <- sim_hot 
  
  df_0 <- bind_rows(df_0, df_attr)
  df_1 <- bind_rows(df_1, sim)
  
}

attr_frac_df <- df_0 %>% 
  mutate(scen="obsclim",
         sim="gswp3-w5e5",
         range ="mildly hot") %>% 
  write.csv(file= here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                       "obsclim_gswp3-w5e5_mild_hot.csv"),
            row.names=FALSE)

df_1 <- df_1 %>% 
  mutate(range ="mildly hot") %>% 
  write.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                 "Factual_gswp3-w5e5_mild_hot.csv"))


########################################################
######           0-4  MILDLY  COLD                ###### 
########################################################
library(reshape2)
library(here)
library(survival)
library(dlnm)
library(tidyverse)
options(scipen=999)
source(here("Codes/attrdl.R"))

here::here()

load(here("Data/Data_1st stage/data_analysis_gswp3-w5e5.RData"))

data <- neonat_mort  %>%
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,tmean_pct_lag2, tmean_pct_lag3, 
            tmean_pct_lag4, tmean_pct_lag5, tmean_pct_lag6, tmean_pct_lag7))   

data_neonat <- data %>%   # Number of neonat observations
  summarise(neonat_mort = sum(neonat_mort))

data <- data %>% 
  dplyr::select(c(CountryName, year, neonat_mort, strata_id,tmean_pct_lag0, tmean_pct_lag1,
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


model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)

#Find the minimum mortality temperature
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)
cen <- cp$predvar[which.min(cp$allRRfit)] 
pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)
countries <- unique(data$CountryName)


df_0 = NULL
df_1 = NULL

for(j in countries){
  country <- j
  df_2 <- data %>%
    dplyr::filter(CountryName == j)   
  if (nrow(df_2) == 0) next
  Qneonat <- df_2  %>% 
    dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
  Qneonat <- as.matrix(Qneonat) 
  Qneonat <- unname(Qneonat)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
  #95% CIs calculated
  
  attr_f_cold <- attrdl(Qneonat, cb.temp,df_2$neonat_mort, model, type="af", 
                        cen=cen, range=c(25, cen))  
  attr_f_cold <- as.data.frame(attr_f_cold)
  set.seed(42)
  CIs_cold <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
             type="af", cen=cen,range=c(25, cen), sim=TRUE, nsim=10000),
             c(2.5,97.5)/100, na.rm=TRUE) 
  
  CIs_cold <- as.data.frame(t(CIs_cold))
  CIs_cold <- CIs_cold %>% 
    dplyr::rename('attr_f_cl' = '2.5%',
                  'attr_f_ch' = '97.5%') 
  df_attr <- cbind (attr_f_cold, CIs_cold) %>%
    mutate(country = j)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
  #95% CIs calculated
  
  set.seed(42)
  sim_cold <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                     type="af",cen=cen, range=c(25, cen), sim=TRUE, nsim=10000) 
  sim_cold <- as.data.frame(sim_cold)
  
  sim <- sim_cold 
  sim$country <- j
  sim$scenario <- "Factual"
  sim$dataset <- "gswp3-w5e5"
  
  df_0 <- bind_rows(df_0, df_attr)
  df_1 <- bind_rows(df_1, sim)
  
}

attr_frac_df <- df_0 %>% 
  mutate(scen="obsclim",
         sim="gswp3-w5e5",
         range ="mildly cold") %>% 
  write.csv(file= here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                   "obsclim_gswp3-w5e5_mild_cold.csv"), row.names=FALSE)

df_1 <- df_1 %>% 
  mutate(range ="mildly cold") %>% 
  write.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                 "Factual_gswp3-w5e5_mild_cold.csv"))

########################################################
######           0-5  MODERATELY  COLD            ###### 
########################################################
library(reshape2)
library(here)
library(survival)
library(dlnm)
library(tidyverse)
options(scipen=999)
source(here("Codes/attrdl.R"))

here::here()

load(here("Data/Data_1st stage/data_analysis_gswp3-w5e5.RData"))

data <- neonat_mort  %>%
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,tmean_pct_lag2, tmean_pct_lag3, 
            tmean_pct_lag4, tmean_pct_lag5, tmean_pct_lag6, tmean_pct_lag7))  

data_neonat <- data %>%   # Number of neonat observations
  summarise(neonat_mort = sum(neonat_mort))

data <- data %>% 
  dplyr::select(c(CountryName, year, neonat_mort, strata_id,tmean_pct_lag0, tmean_pct_lag1,
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


model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)

#Find the minimum mortality temperature 
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)


cen <- cp$predvar[which.min(cp$allRRfit)] 

pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)


countries <- unique(data$CountryName)

df_0 = NULL
df_1 = NULL


for(j in countries){
  country <- j
  df_2 <- data %>%
    dplyr::filter(CountryName == j) 
  if (nrow(df_2) == 0) next
  Qneonat <- df_2  %>% 
    dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
  Qneonat <- as.matrix(Qneonat) 
  Qneonat <- unname(Qneonat)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
  
  attr_f_cold <- attrdl(Qneonat, cb.temp,df_2$neonat_mort, model, type="af", cen=cen,
                        range=c(2.5, 25))   
  attr_f_cold <- as.data.frame(attr_f_cold)
  set.seed(42)
  CIs_cold <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
                    type="af", cen=cen,range=c(2.5, 25), sim=TRUE,nsim=10000),
                    c(2.5,97.5)/100, na.rm=TRUE) 
  
  CIs_cold <- as.data.frame(t(CIs_cold))
  CIs_cold <- CIs_cold %>% 
    dplyr::rename('attr_f_cl' = '2.5%',
                  'attr_f_ch' = '97.5%') 
  df_attr <- cbind (attr_f_cold, CIs_cold) %>%
    mutate(country = j)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
  #95% CIs calculated
  
  set.seed(42)
  sim_cold <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                     type="af",cen=cen, range=c(2.5, 25),
                     sim=TRUE, nsim=10000) 
  sim_cold <- as.data.frame(sim_cold)
  
  sim <-sim_cold 
  sim$country <- j
  sim$scenario <- "Factual"
  sim$dataset <- "gswp3-w5e5"
  
  df_0 <- bind_rows(df_0, df_attr)
  df_1 <- bind_rows(df_1, sim)
  
}

attr_frac_df <- df_0 %>% 
  mutate(scen="obsclim",
         sim="gswp3-w5e5",
         range ="moderately cold") %>% 
  write.csv(file= here("Results", "Attributable_neonatal", "gswp3-w5e5",
              "obsclim_gswp3-w5e5_mod_cold.csv"), row.names=FALSE)

df_1 <- df_1 %>% 
  write.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                 "Factual_gswp3-w5e5_mod_cold.csv"))

########################################################
######           0-6  EXTREMELY  COLD             ###### 
########################################################
library(reshape2)
library(here)
library(survival)
library(dlnm)
library(tidyverse)
options(scipen=999)
source(here("Codes/attrdl.R"))

here::here()

load(here("Data/Data_1st stage/data_analysis_gswp3-w5e5.RData"))


data <- neonat_mort  %>%
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,tmean_pct_lag2, tmean_pct_lag3, 
            tmean_pct_lag4, tmean_pct_lag5, tmean_pct_lag6, tmean_pct_lag7))   

data_neonat <- data %>%   # Number of neonat observations
  summarise(neonat_mort = sum(neonat_mort))

data <- data %>% 
  dplyr::select(c(CountryName, year, neonat_mort, strata_id,tmean_pct_lag0, tmean_pct_lag1,
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


model <- clogit(neonat_mort ~ cb.temp + strata(strata_id), data=data)

#Find the minimum mortality temperature 
cenvar = 26.9
cp <- crosspred(cb.temp, model, cen=cenvar, by=0.1, cumul=TRUE)


cen <- cp$predvar[which.min(cp$allRRfit)] 
pred <- crosspred(cb.temp, model, cen=cen, by=0.1, cumul=TRUE)
countries <- unique(data$CountryName)

df_0 = NULL
df_1 = NULL


for(j in countries){
  country <- j
  df_2 <- data %>%
    dplyr::filter(CountryName == j)   
  
  if (nrow(df_2) == 0) next
  Qneonat <- df_2  %>% 
    dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
  Qneonat <- as.matrix(Qneonat) 
  Qneonat <- unname(Qneonat)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
  #95% CIs calculated
  
  attr_f_cold <- attrdl(Qneonat, cb.temp,df_2$neonat_mort, model, type="af", cen=cen,
                        range=c(0, 2.5))   
  attr_f_cold <- as.data.frame(attr_f_cold)
  set.seed(42)
  CIs_cold <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
                              type="af", cen=cen,range=c(0, 2.5), sim=TRUE,
                              nsim=10000),c(2.5,97.5)/100, na.rm=TRUE) 
  
  CIs_cold <- as.data.frame(t(CIs_cold))
  CIs_cold <- CIs_cold %>% 
    dplyr::rename('attr_f_cl' = '2.5%',
                  'attr_f_ch' = '97.5%') 
  df_attr <- cbind (attr_f_cold, CIs_cold) %>%
    mutate(country = j)
  
  # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
  #95% CIs calculated
  
  set.seed(42)
  sim_cold <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                     type="af",cen=cen, range=c(0, 2.5), sim=TRUE,nsim=10000) 
  sim_cold <- as.data.frame(sim_cold)
  
  sim <- sim_cold  
  sim$country <- j
  sim$scenario <- "Factual"
  sim$dataset <- "gswp3-w5e5"
  
  df_0 <- bind_rows(df_0, df_attr)
  df_1 <- bind_rows(df_1, sim)
  
}

attr_frac_df <- df_0 %>% 
  mutate(scen="obsclim",
         sim="gswp3-w5e5",
         range = "extremely cold") %>% 
  write.csv(file= here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                       "obsclim_gswp3-w5e5_extr_cold.csv"),
            row.names=FALSE)

df_1 <- df_1 %>% 
  write.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                 "Factual_gswp3-w5e5_extr_cold.csv"))


########################################################
######              FIGURE                        ###### 
########################################################




extr_h <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                        "obsclim_gswp3-w5e5_extr_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))

mod_h <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                       "obsclim_gswp3-w5e5_mod_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))

mild_h <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                        "obsclim_gswp3-w5e5_mild_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))


hot_attr <- bind_rows(extr_h, mild_h, mod_h) %>% 
  dplyr::rename("attr_f_total" = "attr_f_hot")

extr_c <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                        "obsclim_gswp3-w5e5_extr_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))

mod_c <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                       "obsclim_gswp3-w5e5_mod_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))

mild_c <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                        "obsclim_gswp3-w5e5_mild_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))

cold_attr <- bind_rows(extr_c, mod_c, mild_c) %>% 
  dplyr::rename("attr_f_total" = "attr_f_cold")
library("readxl")
library("readxl")
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 


neonat_data <- read_excel(here("Data/Neonatal_deaths_data_unicef.xlsx"))

#Calculate births by country for the whole period (2001-2019)
births_data <- read_excel(here("Data/Births_data_unicef.xlsx")) %>% 
  dplyr::select(country, Year, births_total) %>% 
  dplyr::mutate(births_total= as.numeric(births_total)) %>% 
  group_by(country) %>% 
  summarise(births_total = sum(births_total))
#Calculate neonats by country for the whole period (2001-2019)
neonat_data <- neonat_data %>% 
  dplyr::select(country, Year, total_neonat_deaths) %>% 
  dplyr::mutate(total_neonat_deaths= as.numeric(total_neonat_deaths)) %>% 
  group_by(country) %>% 
  summarise(total_neonat_deaths = sum(total_neonat_deaths))


total_attr <- bind_rows(hot_attr,cold_attr) %>% 
  left_join(neonat_data) %>%
  left_join(births_data) %>% 
  mutate(attr_f_total = attr_f_total*total_neonat_deaths,
         rate_f_total = (attr_f_total/(births_total*1000))*100000) %>% 
  dplyr::select(country, rate_f_total, range) 


attr_frac_1 <- melt(total_attr[,c("rate_f_total", 
                           "country", "range")], id=c("country", "range")) %>% 
  dplyr::mutate(scen = "obsclim") %>% 
  dplyr::filter(variable =="rate_f_total")

attr_frac_1 <- attr_frac_1 %>% 
  mutate(range= as.factor(range),
         range=factor(range, c("extremely hot", "moderately hot", "mildly hot", 
                          "mildly cold", "moderately cold", "extremely cold"))) 

attr_frac_1 <- left_join(attr_frac_1, continent) 

#Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "gswp3-w5e5", 
          "heat_cold_rate_gswp3-w5e5_6categories.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(attr_frac_1, aes(x=country, y=value, fill=range)) +
  geom_bar(position = "stack", stat="identity", color="black") +  
  scale_fill_manual(values=c('#e85b54', '#eb7c60', '#eba160','#9bd2f2',  '#4193f0', '#415ef0'))+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  labs(fill = "Temperature range")+
  labs(x ="Country", y = "Heat- and cold-related neonatal mortality rate (per 100,000)")+  
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))+
  theme(legend.key.size = unit(0.4, "cm"), legend.title=element_text(size=10))+ 
  coord_flip()

dev.off()


