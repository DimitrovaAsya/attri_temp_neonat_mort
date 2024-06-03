
#################################################################################
#################################################################################
#######                                                                    ######
####### PROJECT: Temperature-related neonatal deaths attributable to       ######
#######          climate change in 29 low- and middle-income countries     ######                         
#######                                                                    ######
#######        CODE: HIA - Neonatal Deaths [gswp3-w5e5]                    ######
################################################################################# 
#################################################################################  


library(reshape2)
library(here)
library(survival)
library(dlnm)
library(tidyverse)
options(scipen=999)


here::here()

load(here("Data/Data_1st stage/data_analysis_gswp3-w5e5.RData"))


########################################################
######           0-1 FACTUAL SCENARIO             ###### 
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
      dplyr::filter(CountryName ==j)
    
    if (nrow(df_2) == 0) next
    Qneonat <- df_2  %>% 
      dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
    Qneonat <- as.matrix(Qneonat) 
    Qneonat <- unname(Qneonat)
    
    # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
    #95% CIs calculated

    attr_f_hot <- attrdl(Qneonat, cb.temp,df_2$neonat_mort, model, type="af", cen=cen,
                           range=c(cen,100))  
     ############################################################# 
    attr_f_hot <- as.data.frame(attr_f_hot)
    set.seed(42)
    CIs_hot <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
                  type="af", cen=cen,range=c(cen,100), sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)
      ####################################
 
    CIs_hot <- as.data.frame(t(CIs_hot))
    CIs_hot <- CIs_hot %>% 
      dplyr::rename('attr_f_hl' = '2.5%',
                    'attr_f_hh' = '97.5%') 

    attr_f_cold <- attrdl(Qneonat, cb.temp, df_2$neonat_mort, model, type="af", cen=cen,
                   range=c(0, cen))
     #################################################
     attr_f_cold <- as.data.frame(attr_f_cold)
     set.seed(42)
     CIs_cold <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
              type="af", cen=cen,range=c(0,cen), sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)
    #################################### 
    
    CIs_cold <- as.data.frame(t(CIs_cold))
    CIs_cold <- CIs_cold %>% 
      dplyr::rename('attr_f_cl' = '2.5%',
                    'attr_f_ch' = '97.5%') 
    df_attr <- cbind (attr_f_hot, CIs_hot, attr_f_cold, CIs_cold) %>%
      mutate(country = j)
    
    # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
    #95% CIs calculated
    set.seed(42)
    sim_hot <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                        type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)
     ##############################################
    sim_hot <- as.data.frame(sim_hot)
    sim_hot$country <- j
    sim_hot$scenario <- "Factual"
    sim_hot$dataset <- "gswp3-w5e5"
    
    set.seed(42)
    sim_cold <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                         type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)
     ##########################################
    sim_cold <- as.data.frame(sim_cold)
    sim <- cbind(sim_hot,sim_cold)  
    
    df_0 <- bind_rows(df_0, df_attr)
    df_1 <- bind_rows(df_1, sim)
  
}

attr_frac_df <- df_0 %>% 
  mutate(scen="obsclim",
         sim="gswp3-w5e5") %>%

write.csv(file= here("Results", "Attributable_neonatal", "gswp3-w5e5", "obsclim_gswp3-w5e5.csv"),
              row.names=FALSE)
  

  write.csv(df_1, here("Results", "Attributable_neonatal", "Monte_Carlo", "Factual_gswp3-w5e5.csv"))

 #################################################
 #Results for all countries
 set.seed(42)
 sim_hot <- attrdl(Qneonat ,cb.temp,data$neonat_mort, model, 
                   type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)
 sim_hot <- as.data.frame(sim_hot)
 # sim_hot$country <- j
 sim_hot$scenario <- "Factual"
 sim_hot$dataset <- "gswp3-w5e5"
 set.seed(42)
 sim_cold <- attrdl(Qneonat ,cb.temp,data$neonat_mort, model, 
                    type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)
 sim_cold <- as.data.frame(sim_cold)
 
 sim <- cbind(sim_hot,sim_cold) 
 write.csv(sim, here("Results", "Attributable_neonatal", "Monte_Carlo", 
                     "Factual_gswp3-w5e5_all_countries.csv"))
 
 


########################################################
######    0-2 COUNTERFACTUAL SCENARIO             ###### 
########################################################

data_f2 <- neonat_mort_counterclim %>% 
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,tmean_pct_lag2, tmean_pct_lag3, 
            tmean_pct_lag4, tmean_pct_lag5, tmean_pct_lag6, tmean_pct_lag7))   

 
countries <- unique(data_f2$CountryName)


df_0 = NULL
df_1 = NULL


for(j in countries){
    country <- j
    df_2 <- data_f2 %>%
      dplyr::filter(CountryName ==j)
   
    if (nrow(df_2) == 0) next
    Qneonat <- df_2  %>% 
      dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
    Qneonat <- as.matrix(Qneonat) 
    Qneonat <- unname(Qneonat)
    
    # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
    #95% CIs calculated
    
     attr_f_hot <- attrdl(Qneonat, cb.temp,df_2$neonat_mort, model, type="af", cen=cen,
                          range=c(cen,100))  
    ############################################################# 
    attr_f_hot <- as.data.frame(attr_f_hot)
     set.seed(42) 
    CIs_hot <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
             type="af", cen=cen,range=c(cen,100), sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)
    ####################################
    
    CIs_hot <- as.data.frame(t(CIs_hot))
    CIs_hot <- CIs_hot %>% 
      dplyr::rename('attr_f_hl' = '2.5%',
                    'attr_f_hh' = '97.5%') 

    attr_f_cold <- attrdl(Qneonat, cb.temp, df_2$neonat_mort, model, type="af", cen=cen,
                          range=c(0, cen))
    #################################################
    attr_f_cold <- as.data.frame(attr_f_cold)
    set.seed(42)
    CIs_cold <- quantile(attrdl(Qneonat,cb.temp,df_2$neonat_mort, model, 
              type="af", cen=cen,range=c(0,cen), sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)
    #################################### 
    
    CIs_cold <- as.data.frame(t(CIs_cold))
    CIs_cold <- CIs_cold %>% 
      dplyr::rename('attr_f_cl' = '2.5%',
                    'attr_f_ch' = '97.5%') 
    df_attr <- cbind (attr_f_hot, CIs_hot, attr_f_cold, CIs_cold) %>%
      mutate(country = j)
    
    # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
    #95% CIs calculated
     set.seed(42)
     sim_hot <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                       type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)
    ##############################################
    sim_hot <- as.data.frame(sim_hot)
    sim_hot$country <- j
    sim_hot$scenario <- "Counterfactual"
    sim_hot$dataset <- "gswp3-w5e5"
    set.seed(42)
    sim_cold <- attrdl(Qneonat ,cb.temp,df_2$neonat_mort, model, 
                        type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)

    sim_cold <- as.data.frame(sim_cold)
    sim <- cbind(sim_hot,sim_cold)  
    
    df_0 <- bind_rows(df_0, df_attr)
    df_1 <- bind_rows(df_1, sim)
    
}
attr_frac_df <- df_0 %>% 
  mutate(scen="counterclim",
         sim="gswp3-w5e5") %>%

 write.csv(file= here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                      "counterclim_gswp3-w5e5.csv"),
          row.names=FALSE)


 write.csv(df_1, here("Results", "Attributable_neonatal", "Monte_Carlo", 
                      "Counterfactual_gswp3-w5e5.csv"))



####################################################################
# All countries

Qneonat <- data_f2  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
Qneonat <- as.matrix(Qneonat) 
Qneonat <- unname(Qneonat)
set.seed(42)
sim_hot <- attrdl(Qneonat ,cb.temp, data_f2$neonat_mort, model, 
                  type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)
sim_hot <- as.data.frame(sim_hot)
# sim_hot$country <- j
sim_hot$scenario <- "Counterfactual"
sim_hot$dataset <- "gswp3-w5e5"
set.seed(42)
sim_cold <- attrdl(Qneonat ,cb.temp, data_f2$neonat_mort, model, 
                   type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)
sim_cold <- as.data.frame(sim_cold)

sim <- cbind(sim_hot,sim_cold) 
write.csv(sim, here("Results", "Attributable_neonatal", "Monte_Carlo", 
                    "Counterfactual_gswp3-w5e5_all_countries.csv"))






