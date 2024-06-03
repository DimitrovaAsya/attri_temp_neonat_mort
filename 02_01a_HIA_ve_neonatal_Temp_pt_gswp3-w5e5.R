#################################################################################
#################################################################################
#######                                                                    ######
####### PROJECT: Temperature-related neonatal deaths attributable to       ######
#######          climate change in 29 low- and middle-income countries     ######                         
#######                                                                    ######
#######   CODE: HIA - Very Early Neonatal Deaths [gswp3-w5e5]              ######
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
######           0-1 FACTUAL SCENARIO             ###### 
########################################################

data <- ve_neonat %>%
  mutate(year  = format(date, format = "%Y"),
         month = format(date, format = "%m")) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_pct_lag0, tmean_pct_lag1,
            tmean_pct_lag2))   

data <- data %>% 
  dplyr::select(c(ve_neonat, clim.zone, strata_id, CountryName, year, tmean_pct_lag0, tmean_pct_lag1,
                  tmean_pct_lag2)) 

#Define a matrix of exposure histories for each observation (7 days lag)
Qneonat <- data  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag2) 

Qneonat <- as.matrix(Qneonat) 
Qneonat <- unname(Qneonat)
arglag <- list(fun="ns",knots=logknots(2,nk=1))

cb.temp <- crossbasis(Qneonat, lag=2, argvar=list(fun="ns", 
            knots = quantile(data$tmean_pct_lag0, c(10/100), na.rm=T),  
            Bound=range(data$tmean_pct_lag0, na.rm=T)), arglag)

model <- clogit(ve_neonat~cb.temp + strata(strata_id), data=data)

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
       dplyr::filter(CountryName == j)   #Evaluate impact of extreme heat only

    if (nrow(df_2) == 0) next
    Qneonat <- df_2  %>% 
      dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
    Qneonat <- as.matrix(Qneonat) 
    Qneonat <- unname(Qneonat)
    
    # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD
    #95% CIs calculated
    
    attr_f_hot <- attrdl(Qneonat, cb.temp,df_2$ve_neonat, model, type="af", cen=cen,
                          range=c(cen,100))  
    attr_f_hot <- as.data.frame(attr_f_hot)
    set.seed(42)
    CIs_hot <- quantile(attrdl(Qneonat,cb.temp,df_2$ve_neonat, model, 
                  type="af", cen=cen, range=c(cen,100), sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE)  
    CIs_hot <- as.data.frame(t(CIs_hot))
    CIs_hot <- CIs_hot %>% 
       dplyr::rename('attr_f_hl' = '2.5%',
                     'attr_f_hh' = '97.5%')
    
    attr_f_cold <- attrdl(Qneonat, cb.temp,df_2$ve_neonat, model, type="af", cen=cen,
                         range=c(0, cen))   
    attr_f_cold <- as.data.frame(attr_f_cold)
    set.seed(42)
    CIs_cold <- quantile(attrdl(Qneonat,cb.temp,df_2$ve_neonat, model, 
                  type="af", cen=cen,range=c(0,cen), sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE) 
    
    CIs_cold <- as.data.frame(t(CIs_cold))
    CIs_cold <- CIs_cold %>% 
      dplyr::rename('attr_f_cl' = '2.5%',
                    'attr_f_ch' = '97.5%') 
    df_attr <- cbind (attr_f_hot, CIs_hot, attr_f_cold, CIs_cold) %>%
         mutate(country = j)
    
    # ATTRIBUTABLE FRACTION COMPONENT DUE TO HEAT AND COLD - 1000 MC simulations
    #95% CIs calculated
    
    set.seed(42)
    sim_hot <- attrdl(Qneonat ,cb.temp,df_2$ve_neonat, model, 
             type="af",cen=cen, range=c(cen,100), sim=TRUE, nsim=10000) 
  
    sim_hot <- as.data.frame(sim_hot)
    sim_hot$country <- j
    sim_hot$scenario <- "Factual"
    sim_hot$dataset <- "gswp3-w5e5"
    set.seed(42)
    sim_cold <- attrdl(Qneonat ,cb.temp,df_2$ve_neonat, model, 
                  type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)  
    sim_cold <- as.data.frame(sim_cold)
    
    sim <- cbind(sim_hot,sim_cold)  
    
    df_0 <- bind_rows(df_0, df_attr)
    df_1 <- bind_rows(df_1, sim)
    
  }
  
attr_frac_df <- df_0 %>% 
  mutate(scen="obsclim",
         sim="gswp3-w5e5") %>% 
  write.csv(file= here("Results", "Attributable_ve_neonats", "gswp3-w5e5", "obsclim_gswp3-w5e5.csv"),
            row.names=FALSE)


write.csv(df_1, here("Results", "Attributable_ve_neonats", "Monte_Carlo", "Factual_gswp3-w5e5.csv"))


#################################################
#Results for all countries
set.seed(42)
sim_hot <- attrdl(Qneonat ,cb.temp,data$ve_neonat, model, 
                  type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)
sim_hot <- as.data.frame(sim_hot)

sim_hot$scenario <- "Factual"
sim_hot$dataset <- "gswp3-w5e5"
set.seed(42)
sim_cold <- attrdl(Qneonat ,cb.temp,data$ve_neonat, model, 
                   type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)
sim_cold <- as.data.frame(sim_cold)

sim <- cbind(sim_hot,sim_cold) 
write.csv(sim, here("Results", "Attributable_ve_neonats", "Monte_Carlo", "Factual_gswp3-w5e5_all_countries.csv"))

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
    Qneonat <- df_2  %>% 
      dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
    Qneonat <- as.matrix(Qneonat) 
    Qneonat <- unname(Qneonat)
  
    attr_f_hot <- attrdl(Qneonat, cb.temp,df_2$ve_neonat, model, type="af", cen=cen,
                         range=c(cen,100))   
    attr_f_hot <- as.data.frame(attr_f_hot)
    
    set.seed(42)
    CIs_hot <- quantile(attrdl(Qneonat,cb.temp,df_2$ve_neonat, model, 
                               type="af", cen=cen,range=c(cen,100), sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE) 
    CIs_hot <- as.data.frame(t(CIs_hot))
    CIs_hot <- CIs_hot %>% 
      dplyr::rename('attr_f_hl' = '2.5%',
                    'attr_f_hh' = '97.5%') 
    attr_f_cold <- attrdl(Qneonat, cb.temp,df_2$ve_neonat, model, type="af", cen=cen,
                          range=c(0,cen))  
    attr_f_cold <- as.data.frame(attr_f_cold)
    
    set.seed(42)
    CIs_cold <- quantile(attrdl(Qneonat,cb.temp,df_2$ve_neonat, model, 
                                type="af", cen=cen,range=c(0,cen), sim=TRUE,nsim=10000),c(2.5,97.5)/100, na.rm=TRUE) 
    
    CIs_cold <- as.data.frame(t(CIs_cold))
    CIs_cold <- CIs_cold %>% 
      dplyr::rename('attr_f_cl' = '2.5%',
                    'attr_f_ch' = '97.5%') 
    df_attr <- cbind (attr_f_hot, CIs_hot, attr_f_cold, CIs_cold) %>%
      mutate(country = j)
    

    set.seed(42)
    sim_hot <- attrdl(Qneonat ,cb.temp,df_2$ve_neonat, model, 
                      type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)  

    sim_hot <- as.data.frame(sim_hot)
    sim_hot$country <- j
    sim_hot$scenario <- "Counterfactual"
    sim_hot$dataset <- "gswp3-w5e5"
    
    set.seed(42)
    sim_cold <- attrdl(Qneonat ,cb.temp,df_2$ve_neonat, model, 
                       type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)  
    sim_cold <- as.data.frame(sim_cold)
    
    sim <- cbind(sim_hot,sim_cold) 

  df_0 <- bind_rows(df_0, df_attr)
  df_1 <- bind_rows(df_1, sim)
}



attr_frac_df <- df_0 %>% 
  mutate(scen="counterclim",
         sim="gswp3-w5e5") %>% 
write.csv(file= here("Results", "Attributable_ve_neonats", "gswp3-w5e5", "counterclim_gswp3-w5e5.csv"),
            row.names=FALSE)


write.csv(df_1, here("Results", "Attributable_ve_neonats", "Monte_Carlo", "Counterfactual_gswp3-w5e5.csv"))

####################################################################
# All countries

Qneonat <- data_f2  %>% 
  dplyr::select(tmean_pct_lag0:tmean_pct_lag2)  
Qneonat <- as.matrix(Qneonat) 
Qneonat <- unname(Qneonat)
set.seed(42)
sim_hot <- attrdl(Qneonat ,cb.temp,data_f2$ve_neonat, model, 
                  type="af",cen=cen, range=c(cen,100), sim=TRUE,nsim=10000)
sim_hot <- as.data.frame(sim_hot)
sim_hot$scenario <- "Counterfactual"
sim_hot$dataset <- "gswp3-w5e5"
set.seed(42)
sim_cold <- attrdl(Qneonat ,cb.temp,data_f2$ve_neonat, model, 
                   type="af",cen=cen, range=c(0,cen), sim=TRUE,nsim=10000)
sim_cold <- as.data.frame(sim_cold)

sim <- cbind(sim_hot,sim_cold) 
write.csv(sim, here("Results", "Attributable_ve_neonats", "Monte_Carlo", "Counterfactual_gswp3-w5e5_all_countries.csv"))


############################################################################
### Calculate ve_neonats attributable to climate change in this scenario  ##
############################################################################
library(reshape2)
library(here)
library(tidyverse)
options(scipen=999)


here::here()
options(digits = 5) 
ve_neonats_f <- read.csv(here("Results", "Attributable_ve_neonats","gswp3-w5e5",
                               "obsclim_gswp3-w5e5.csv")) %>%
  dplyr::select(country, attr_f_hot:sim) %>% 
  dplyr::select(-scen) %>%
  dplyr::rename(
    attr_f_h_f = attr_f_hot,
    attr_f_hl_f = attr_f_hl,
    attr_f_hh_f = attr_f_hh,
    attr_f_c_f = attr_f_cold,
    attr_f_cl_f = attr_f_cl,
    attr_f_ch_f = attr_f_ch)

ve_neonats_cf <- read.csv(here("Results", "Attributable_ve_neonats", "gswp3-w5e5", "counterclim_gswp3-w5e5.csv")) %>%
  dplyr::select(country, attr_f_hot:sim) %>%
  dplyr::select(-scen) %>%
  dplyr::rename(
    attr_f_h_cf = attr_f_hot,
    attr_f_hl_cf = attr_f_hl,
    attr_f_hh_cf = attr_f_hh,
    attr_f_c_cf = attr_f_cold,
    attr_f_cl_cf = attr_f_cl,
    attr_f_ch_cf = attr_f_ch)




ve_neonats <- ve_neonats_f %>%
  left_join(ve_neonats_cf, by=c("country"="country", "sim"="sim")) %>% 
  mutate(
    attr_frac_cold = attr_f_c_f - attr_f_c_cf,
    attr_frac_hot = attr_f_h_f - attr_f_h_cf,
    attr_frac_coldl = attr_f_cl_f - attr_f_cl_cf,
    attr_frac_hotl = attr_f_hl_f - attr_f_hl_cf,
    attr_frac_coldh = attr_f_ch_f - attr_f_ch_cf,
    attr_frac_hoth = attr_f_hh_f - attr_f_hh_cf)


attr_frac_1 <- melt(ve_neonats[,c("country", "attr_frac_hot", "attr_frac_cold")], id=c("country")) 
attr_frac_2 <- melt(ve_neonats[,c("country", "attr_frac_hotl", 
                                  "attr_frac_coldl")], id=c("country")) %>% 
  dplyr::rename("value_l" ="value") %>% 
  mutate(variable=as.character(variable)) %>% 
  mutate(variable= substr(variable,1,nchar(variable)-1)) 

attr_frac_3 <- melt(ve_neonats[,c("country", "attr_frac_hoth", 
                                  "attr_frac_coldh")], id=c("country")) %>% 
  dplyr::rename("value_h" ="value") %>% 
  mutate(variable=as.character(variable)) %>% 
  mutate(variable= substr(variable,1,nchar(variable)-1))


attr_frac_4 <- merge(attr_frac_1, attr_frac_2, 
                     by= c("country", "variable"), all.x = TRUE,  all.y=TRUE)
attr_frac <- merge(attr_frac_4, attr_frac_3, 
                   by= c("country", "variable"), all.x = TRUE,  all.y=TRUE) %>% 
  mutate_if(is.numeric, ~ . * 100) %>%  #multiply by 100 to get the fraction in %
  mutate(variable = recode(variable,
                           "attr_frac_hot"  = "Hot", "attr_frac_cold" = "Cold"))

#Plot attributable fraction by country (geom bar figure with both hot and cold)
jpeg(here("Figures", "ve_neonats" , "Stage II", "gswp3-w5e5", "attr_farc_cc_gswp3-w5e5.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

ggplot(attr_frac, aes(x=country, y=value, fill=variable)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_errorbar( aes(x=country, ymin=value_l, ymax=value_h), width=0.4, colour="black",
                 position=position_dodge(.9), alpha=0.7, size=0.5)+
  labs(title="Fraction of ve_neonats attributable to climate change [gswp3-w5e5]\nLag 0-2",
       x ="Country", y = "Attributable fraction")+
  theme_minimal()+
  scale_fill_manual("",values = c("Hot" = "#eb8060", "Cold" = "#a1e9f0"))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()
############################################################################
############################################################################
library(reshape2)
library(here)
library(tidyverse)
library("readxl")
options(scipen=999)


here::here()
options(digits = 5) 
ve_neonats_f <- read.csv(here("Results", "Attributable_ve_neonats", "gswp3-w5e5", "obsclim_gswp3-w5e5.csv")) %>%
  dplyr::select(country, attr_f_hot:sim) %>% 
  dplyr::select(-scen) %>%
  dplyr::rename(
    attr_f_h_f = attr_f_hot,
    attr_f_hl_f = attr_f_hl,
    attr_f_hh_f = attr_f_hh,
    attr_f_c_f = attr_f_cold,
    attr_f_cl_f = attr_f_cl,
    attr_f_ch_f = attr_f_ch)

ve_neonats_cf <- read.csv(here("Results", "Attributable_ve_neonats", "gswp3-w5e5", "counterclim_gswp3-w5e5.csv")) %>%
  dplyr::select(country, attr_f_hot:sim) %>%
  dplyr::select(-scen) %>%
  dplyr::rename(
    attr_f_h_cf = attr_f_hot,
    attr_f_hl_cf = attr_f_hl,
    attr_f_hh_cf = attr_f_hh,
    attr_f_c_cf = attr_f_cold,
    attr_f_cl_cf = attr_f_cl,
    attr_f_ch_cf = attr_f_ch)




ve_neonats <- ve_neonats_f %>%
  left_join(ve_neonats_cf, by=c("country"="country", "sim"="sim")) %>% 
  mutate(
    attr_frac_c_cc = attr_f_c_f - attr_f_c_cf,
    attr_frac_h_cc = attr_f_h_f - attr_f_h_cf,
    attr_frac_cl_cc = attr_f_cl_f - attr_f_cl_cf,
    attr_frac_hl_cc = attr_f_hl_f - attr_f_hl_cf,
    attr_frac_ch_cc = attr_f_ch_f - attr_f_ch_cf,
    attr_frac_hh_cc = attr_f_hh_f - attr_f_hh_cf)


ve_neonats_hot <- ve_neonats %>% 
  select(country, attr_f_h_cf, attr_f_h_f, attr_frac_h_cc, attr_frac_hl_cc, attr_frac_hh_cc) %>%
  mutate_if(is.numeric, ~ . * 100) %>%  #multiply by 100 to get the fraction in %
  dplyr::rename("Counterfactual" = "attr_f_h_cf",
                "Factual" = "attr_f_h_f") 
ve_neonats_cold <- ve_neonats %>% 
  select(country, attr_f_c_cf, attr_f_c_f, attr_frac_c_cc, attr_frac_cl_cc, attr_frac_ch_cc) %>%
  mutate_if(is.numeric, ~ . * 100) %>%  #multiply by 100 to get the fraction in %
  dplyr::rename("Counterfactual" = "attr_f_c_cf",
                "Factual" = "attr_f_c_f") 

continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 


ve_neonats_cold <- left_join(ve_neonats_cold, continent)

ve_neonats_hot <- left_join(ve_neonats_hot, continent)
attr_frac <- melt(ve_neonats_hot, id.vars= c("country", "continent"), 
                  measure.vars = c("Factual", "Counterfactual", "attr_frac_h_cc", "attr_frac_hl_cc", "attr_frac_hh_cc")) 



jpeg(here("Figures", "ve_neonats" , "Stage II", "gswp3-w5e5", "attr_farc_heat_gswp3-w5e5.jpg"),
     units = 'in', res = 300, width = 10, height = 6)
 
ve_neonats_hot %>%
  ggplot(aes(x = Counterfactual, y = fct_reorder(country, Factual))) +
  geom_col(aes(x = Factual), fill="#d9716c", alpha = 0.5, width = 0.8) +
  geom_col(width = 0.4, show.legend= FALSE, fill="#a3524e") +
  facet_grid(continent ~ ., scales = "free", space = "free", ) +
  labs(title= "gswp3-w5e5",   
       x ="Percentage", y = "")+
  theme(axis.text.x= element_text(face= "bold"))+
  geom_vline(xintercept=0)+
  theme_minimal()
dev.off()





jpeg(here("Figures", "ve_neonats" , "Stage II", "gswp3-w5e5", "attr_farc_cc_gswp3-w5e5.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

ggplot(ve_neonats_hot, aes(x=fct_reorder(country, attr_frac_h_cc), y=attr_frac_h_cc)) + 
  geom_pointrange(aes(ymin=attr_frac_hl_cc, ymax=attr_frac_hh_cc), size=0.6, colour="#a32d12" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+

  labs(title= "gswp3-w5e5", 
       x ="", y = "Percentage")+
  theme_minimal()
dev.off()


attr_frac <- melt(ve_neonats_cold, id.vars=c("country", "continent"), 
                  measure.vars =c("Factual", "Counterfactual", "attr_frac_c_cc", "attr_frac_cl_cc", "attr_frac_ch_cc")) 



jpeg(here("Figures", "ve_neonats" , "Stage II", "gswp3-w5e5", "attr_farc_cold_gswp3-w5e5.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

ve_neonats_cold %>%
  ggplot(aes(x = Counterfactual, y = fct_reorder(country, Factual))) +
  geom_col(aes(x = Factual), fill="#2b50cc", alpha = 0.5, width = 0.8) +
  geom_col(width = 0.4, show.legend= FALSE, fill="#3454bf") +
  geom_vline(xintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  labs(title= "gswp3-w5e5",  
       x ="Percentage", y = "")+
  theme_minimal()
dev.off()




jpeg(here("Figures", "ve_neonats" , "Stage II", "gswp3-w5e5", "attr_farc_cc_gswp3-w5e5_cold.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

ggplot(ve_neonats_cold, aes(x=fct_reorder(country, attr_frac_c_cc), y=attr_frac_c_cc)) + 
  geom_pointrange(aes(ymin=attr_frac_cl_cc, ymax=attr_frac_ch_cc), size=0.6, colour="#2d52cc" )+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  geom_hline(yintercept=0)+
  coord_flip()+
  labs(title= "gswp3-w5e5", 
       x ="", y = "Percentage")+
  theme_minimal()
dev.off()



######################################################################
######           OPTIMAL TEMP VS. Latitude and average temp    ###### 
######################################################################
library(reshape2)
library(here)
library(tidyverse)
library("readxl")
options(scipen=999)
#Upload data with latitude & continent
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 


#Upload data with average annual temperature for each country 
#for the study period
#### Load the daily temperature data

mean_temp <- read.csv(here("Data/Data_1st stage/Linked_data_ISIMIP/gswp3-w5e5/Cluster/tmean_gswp3-w5e5_annual.csv")) [,-1]
country_identifier <- read.csv(here("Data/Data_1st stage/country_identifier.csv")) [,-1]



#Upload data with optimum temperature in absolute terms for each country

MMT <- read.csv(here("Data/Data_1st stage/Linked_data_ISIMIP/gswp3-w5e5/Cluster/MMT_gswp3-w5e5.csv")) [,-1]


str(country_identifier$SurveyId)
all <- country_identifier %>% 
  left_join(MMT) %>% 
  filter(tmean_pct== 41)%>%  # filter only the ref temp pct (for ve_neonats)
  left_join(mean_temp) %>%
  left_join(continent, by=c("CountryName" ="country")) %>% 
  mutate(latitude=as.numeric(latitude),
         latitude=round(latitude, digits=0)) %>% 
  dplyr::rename(Country= CountryName,
                Continent=continent)
library(ggplot2)
library(hrbrthemes)

g1 <- ggplot(all, aes(x=mean_temp, y=mean_pct, shape= Continent, colour =  Country)) + 
  geom_point(size=6) +
  labs(title="A.",
       x ="Average temperature (C°)", y = "Minimum-ve_neonat temperature (C°)")+
  theme_ipsum()

g2 <- ggplot(all, aes(x=latitude, y=mean_pct, shape= Continent, colour =  Country)) + 
  geom_point(size=6) +
  labs(title="B.",
       x ="Average latitude", y = "Minimum-ve_neonat temperature (C°)")+
  geom_vline(xintercept = 0, lty=2)+
  theme_ipsum()


library("patchwork")

jpeg(here("Figures","ve_neonats", "Stage I", "gswp3-w5e5", "MST_gswp3-w5e5.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

g1 + g2 + plot_layout(guides = "collect")
dev.off()