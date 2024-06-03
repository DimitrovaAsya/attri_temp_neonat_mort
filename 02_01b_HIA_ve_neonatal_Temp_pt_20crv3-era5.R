#################################################################################
#################################################################################
#######                                                                    ######
####### PROJECT: Temperature-related neonatal deaths attributable to       ######
#######          climate change in 29 low- and middle-income countries     ######                         
#######                                                                    ######
#######   CODE: HIA - Very Early Neonatal Deaths [20crv3-era5]             ######
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


load(here("Data/Data_1st stage/data_analysis_20crv3-era5.RData"))

########################################################
######           0-1 FACTUAL SCENARIO             ###### 
########################################################
data <- ve_neonat %>%
  mutate(year  = format(date, format = "%Y"),
         year=as.numeric(year),
         month = format(date, format = "%m"),
         month=as.numeric(month)) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,
            tmean_pct_lag2))   

data_neonat <- data %>%   # Number of ve_neonat observations
  summarise(ve_neonat = sum(ve_neonat))
data <- data %>% 
  dplyr::select(c(CountryName, year, clim.zone, ve_neonat, strata_id,tmean_pct_lag0, tmean_pct_lag1,
                  tmean_pct_lag2))


#Define a matrix of exposure histories for each observation (2 days lag)
Qve_neonat <- data  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag2) 

Qve_neonat <- as.matrix(Qve_neonat) 
Qve_neonat <- unname(Qve_neonat)
arglag <- list(fun="ns", knots=logknots(2,nk=1))  

cb.temp <- crossbasis(Qve_neonat, lag=2, argvar=list(fun="ns", 
                knots = quantile(data$tmean_pct_lag0, c(10/100), na.rm=T),   
                Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

model <- clogit(ve_neonat~cb.temp + strata(strata_id), data=data)

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
       dplyr::filter(CountryName ==j)  
               
    if (nrow(df_2) == 0) next
    Qve_neonat <- df_2  %>% 
      dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
    Qve_neonat <- as.matrix(Qve_neonat) 
    Qve_neonat <- unname(Qve_neonat)

    # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
    #95% CIs calculated
    attr_f_hot <- attrdl(Qve_neonat, cb.temp,df_2$ve_neonat, model, type="af", 
                         cen=cen, range=c(cen,100))
    attr_f_hot <- as.data.frame(attr_f_hot)
    set.seed(42)
    CIs_hot <- quantile(attrdl(Qve_neonat,cb.temp,df_2$ve_neonat, model, 
                      type="af", cen=cen,range=c(cen,100), sim=TRUE,nsim=10000),
                      c(2.5,97.5)/100, na.rm=TRUE)
    CIs_hot <- as.data.frame(t(CIs_hot))
    CIs_hot <- CIs_hot %>% 
      dplyr::rename('attr_f_hl' = '2.5%',
                    'attr_f_hh' = '97.5%') 
    
    attr_f_cold <- attrdl(Qve_neonat, cb.temp,df_2$ve_neonat, model, type="af", 
                          cen=cen, range=c(0,cen))
    attr_f_cold <- as.data.frame(attr_f_cold)
    set.seed(42)
    CIs_cold <- quantile(attrdl(Qve_neonat,cb.temp,df_2$ve_neonat, model, 
                                type="af", cen=cen,range=c(0,cen), sim=TRUE, 
                                nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)
    
    CIs_cold <- as.data.frame(t(CIs_cold))
    CIs_cold <- CIs_cold %>% 
      dplyr::rename('attr_f_cl' = '2.5%',
                    'attr_f_ch' = '97.5%') 
    df_attr <- cbind (attr_f_hot, CIs_hot, attr_f_cold, CIs_cold) %>%
      mutate(country = j)
    
    # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
    #95% CIs calculated
    set.seed(42)
    sim_hot <- attrdl(Qve_neonat ,cb.temp,df_2$ve_neonat, model, 
                      type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)
    sim_hot <- as.data.frame(sim_hot)
    sim_hot$country <- j
    sim_hot$scenario <- "Factual"
    sim_hot$dataset <- "20crv3-era5"
    
    set.seed(42)
    sim_cold <- attrdl(Qve_neonat ,cb.temp,df_2$ve_neonat, model, 
                       type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)
    sim_cold <- as.data.frame(sim_cold)
    
    sim <- cbind(sim_hot,sim_cold)  

  
    df_0 <- bind_rows(df_0, df_attr)
    df_1 <- bind_rows(df_1, sim)
}

attr_frac_df <- df_0 %>% 
  mutate(scen="obsclim",
         sim="20crv3-era5") %>% 
  write.csv(file= here("Results", "Attributable_ve_neonats", "20crv3-era5", "obsclim_20crv3-era5.csv"),
            row.names=FALSE)

write.csv(df_1, here("Results", "Attributable_ve_neonats", "Monte_Carlo", "Factual_20crv3-era5.csv"))


#Results for all countries
set.seed(42)
sim_hot <- attrdl(Qve_neonat ,cb.temp,data$ve_neonat, model, 
                  type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)
sim_hot <- as.data.frame(sim_hot)
# sim_hot$country <- j
sim_hot$scenario <- "Factual"
sim_hot$dataset <- "20crv3-era5"
set.seed(42)
sim_cold <- attrdl(Qve_neonat ,cb.temp,data$ve_neonat, model, 
                   type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)
sim_cold <- as.data.frame(sim_cold)

sim <- cbind(sim_hot,sim_cold) 
write.csv(sim, here("Results", "Attributable_ve_neonats", "Monte_Carlo", 
                    "Factual_20crv3-era5_all_countries.csv"))

########################################################
######    0-2 COUNTERFACTUAL SCENARIO             ###### 
########################################################

data_f2 <- ve_neonat_counterclim %>% 
  mutate(date=as.Date(date)) %>% 
  mutate(year= format(date, format = "%Y")) %>% 
  mutate(year=as.numeric(year)) %>%
  filter(year_int - year <= 14) %>%
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,
            tmean_pct_lag2)) %>% 
  dplyr::select(c(ve_neonat, clim.zone, strata_id, CountryName, year, tmean_pct_lag0, tmean_pct_lag1,
                  tmean_pct_lag2))

countries <- unique(data_f2$CountryName)


df_0 = NULL
df_1 = NULL

for(j in countries){
    country <- j
    df_2 <- data_f2 %>%
       dplyr::filter(CountryName ==j)              
    if (nrow(df_2) == 0) next
    
    Qve_neonat <- df_2  %>% 
      dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
    Qve_neonat <- as.matrix(Qve_neonat) 
    Qve_neonat <- unname(Qve_neonat)
    
    
    attr_f_hot <- attrdl(Qve_neonat, cb.temp,df_2$ve_neonat, model, type="af", cen=cen,
                         range=c(cen,100))
    attr_f_hot <- as.data.frame(attr_f_hot)
    set.seed(42)
    CIs_hot <- quantile(attrdl(Qve_neonat,cb.temp,df_2$ve_neonat, model, 
                        type="af", cen=cen,range=c(cen,100), 
                        sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)
    CIs_hot <- as.data.frame(t(CIs_hot))
    CIs_hot <- CIs_hot %>% 
      dplyr::rename('attr_f_hl' = '2.5%',
                    'attr_f_hh' = '97.5%') 
    attr_f_cold <- attrdl(Qve_neonat, cb.temp,df_2$ve_neonat, model, type="af", 
                          cen=cen, range=c(0,cen))
    attr_f_cold <- as.data.frame(attr_f_cold)
    set.seed(42)
    CIs_cold <- quantile(attrdl(Qve_neonat,cb.temp,df_2$ve_neonat, model, 
                            type="af", cen=cen,range=c(0,cen), 
                            sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)
    
    CIs_cold <- as.data.frame(t(CIs_cold))
    CIs_cold <- CIs_cold %>% 
      dplyr::rename('attr_f_cl' = '2.5%',
                    'attr_f_ch' = '97.5%') 
    df_attr <- cbind (attr_f_hot, CIs_hot, attr_f_cold, CIs_cold) %>%
      mutate(country = j)
    
    set.seed(42)
    sim_hot <- attrdl(Qve_neonat ,cb.temp,df_2$ve_neonat, model, 
                      type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)
    sim_hot <- as.data.frame(sim_hot)
    sim_hot$country <- j
    sim_hot$scenario <- "Counterfactual"
    sim_hot$dataset <- "20crv3-era5"
    
    set.seed(42)
    sim_cold <- attrdl(Qve_neonat ,cb.temp,df_2$ve_neonat, model, 
                       type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)
    sim_cold <- as.data.frame(sim_cold)
    
    sim <- cbind(sim_hot,sim_cold)  
  
    df_0 <- bind_rows(df_0, df_attr)
    df_1 <- bind_rows(df_1, sim)
}

attr_frac_df <- df_0 %>% 
  mutate(scen="counterclim",
         sim="20crv3-era5") %>% 
  write.csv(file= here("Results", "Attributable_ve_neonats", "20crv3-era5", 
                       "counterclim_20crv3-era5.csv"), row.names=FALSE)


write.csv(df_1, here("Results", "Attributable_ve_neonats", "Monte_Carlo", 
                     "Counterfactual_20crv3-era5.csv"))



# All countries
Qve_neonat <- data_f2  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
Qve_neonat <- as.matrix(Qve_neonat) 
Qve_neonat <- unname(Qve_neonat)
set.seed(42)
sim_hot <- attrdl(Qve_neonat ,cb.temp, data_f2$ve_neonat, model, 
                  type="af",cen=cen, range=c(cen,100), sim=TRUE, nsim=10000)
sim_hot <- as.data.frame(sim_hot)
sim_hot$scenario <- "Counterfactual"
sim_hot$dataset <- "20crv3-era5"
set.seed(42)
sim_cold <- attrdl(Qve_neonat ,cb.temp, data_f2$ve_neonat, model, 
                   type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)
sim_cold <- as.data.frame(sim_cold)

sim <- cbind(sim_hot,sim_cold) 
write.csv(sim, here("Results", "Attributable_ve_neonats", "Monte_Carlo", 
                    "Counterfactual_20crv3-era5_all_countries.csv"))



