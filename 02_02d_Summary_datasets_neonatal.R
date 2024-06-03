#################################################################################
#################################################################################
#######                                                                    ######
####### PROJECT: Temperature-related neonatal deaths attributable to       ######
#######          climate change in 29 low- and middle-income countries     ######                         
#######                                                                    ######
#######           CODE: SUMMARY RESULTS - Neonatal Deaths                  ######
################################################################################# 
#################################################################################  


######### 1. Combine MC simulations #########

#Calculates mean share of neonatal deaths from hot and cold and CIs across
#in all neonatal deaths all 3 datasets
library(reshape2)
library(here)
library(tidyverse)
options(scipen=999)
here::here()


fac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                       "Factual_gswp3-w5e5.csv"))[,-1]
fac_2 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                       "Factual_20crv3-era5.csv"))[,-1]
fac_3 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                       "Factual_20crv3-w5e5.csv"))[,-1]

fac <- rbind.data.frame(fac_1, fac_2, fac_3) %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)



countfac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                            "Counterfactual_gswp3-w5e5.csv"))[,-1]
countfac_2 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                            "Counterfactual_20crv3-era5.csv"))[,-1]
countfac_3 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                            "Counterfactual_20crv3-w5e5.csv"))[,-1]

countfac <- rbind.data.frame(countfac_1, countfac_2, countfac_3) %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

#Put together the results of all MC simulations, for all reanalysis data sets for both scenarios
summary <- fac %>% 
  cbind(countfac) 
summary <- summary[, !duplicated(colnames(summary))]  
summary <- summary %>%  
  select(-scenario) 


######### 2. Shares and rates temp-attributable mortality by country in factual #########
#Calculate mean and upper and lower bounds for share of heat- and cold-related deaths
#in total neonatal mortality for each country
set.seed(42)
summary_countries <- summary %>% 
  group_by(country) %>% 
  summarise_at(vars(sim_hot_f, sim_hot_cf, sim_cold_f, sim_cold_cf),
               list(mean=mean, upp =~quantile(., probs = 0.95),
                    low =~quantile(., probs = 0.05))) %>% #, 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %


library("readxl")
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 

#Calculate neonatal deaths by country for the whole period (2001-2019)
neonat_data <- read_excel(here("Data/Neonatal_deaths_data_unicef.xlsx")) %>% 
  dplyr::select(country, Year, total_neonat_deaths) %>% 
  dplyr::mutate(total_neonat_deaths= as.numeric(total_neonat_deaths)) %>% 
  group_by(country) %>% 
  summarise(total_neonat_deaths = sum(total_neonat_deaths))
#Calculate births by country for the whole period (2001-2019)
births_data <- read_excel(here("Data/Births_data_unicef.xlsx")) %>% 
  dplyr::select(country, Year, births_total) %>% 
  dplyr::mutate(births_total= as.numeric(births_total)) %>% 
  group_by(country) %>% 
  summarise(births_total = sum(births_total))

births_data_UNICEF <- births_data %>% 
  mutate(births_total= births_total*1000)  %>% #Birth are per 1000, so need to multiply by 1000
  dplyr::filter(country %in% c("Burundi", "Ethiopia", "Liberia", "Sierra Leone", "Rwanda",
                               "Uganda", "Cameroon", "Benin", "Guinea", "Angola", "Nigeria",
                               "Senegal", "Tanzania", "Mali", "Zimbabwe", "South Africa",
                               "Malawi", "Zambia", "Philippines", "Timor-Leste", "Bangladesh",
                               "Pakistan", "Nepal", "India","Jordan", "Tajikistan", "Armenia",
                               "Haiti", "Albania")) %>% 
  write.csv(here("Results", "Tables_SM", "births_data_UNICEF.csv"), row.names=FALSE)

#Calculate rate per 1,000 births
summary_countries2 <- summary_countries %>% 
  select(country, sim_hot_f_mean, sim_cold_f_mean, sim_hot_cf_mean,sim_cold_cf_mean) %>% 
  left_join(neonat_data) %>%
  left_join(births_data) %>% 
  mutate(total_heat_f = (sim_hot_f_mean/100)*total_neonat_deaths,
         total_cold_f = (sim_cold_f_mean/100)*total_neonat_deaths,
         total_heat_cf = (sim_hot_cf_mean/100)*total_neonat_deaths,
         total_cold_cf = (sim_cold_cf_mean/100)*total_neonat_deaths) %>% 
  summarise(total_heat_f = sum(total_heat_f), 
            total_heat_cf = sum(total_heat_cf),
            total_cold_cf = sum(total_cold_cf),
            total_cold_f = sum(total_cold_f)),
         rate_heat_ = (total_heat/(births_total*1000))*100000,
         rate_cold = (total_cold/(births_total*1000))*100000,
         rate_both = (total_heat+total_cold)/(births_total*1000)*100000) %>% 
  dplyr::select(country, rate_heat, rate_cold, rate_both, total_heat, total_cold) 

#% increase in heat-related due to climate change for all countries
heat_cc_increase <- summary_countries2 %>% 
  select(total_heat_f, total_heat_cf) %>% 
  mutate(heat_cc_increase =((total_heat_f-total_heat_cf)/total_heat_cf))
#% increase in heat-related due to climate change by country

summary_countries_seperate_heat <- summary_countries %>% 
  select(country, sim_hot_f_mean, sim_cold_f_mean, sim_hot_cf_mean,sim_cold_cf_mean) %>% 
  left_join(neonat_data) %>%
  left_join(births_data) %>% 
  mutate(total_heat_f = (sim_hot_f_mean/100)*total_neonat_deaths,
         total_cold_f = (sim_cold_f_mean/100)*total_neonat_deaths,
         total_heat_cf = (sim_hot_cf_mean/100)*total_neonat_deaths,
         total_cold_cf = (sim_cold_cf_mean/100)*total_neonat_deaths) %>% 
  mutate(heat_cc_increase =((total_heat_f-total_heat_cf)/total_heat_cf)) %>% 
  write.csv(here("Results", "Tables_SM", "Increase_heat-related_neonat_country.csv"), 
            row.names=FALSE)



cold_cc_decrease <- summary_countries2 %>% 
  select(total_cold_f, total_cold_cf) %>% 
  mutate(cold_cc_increase =((total_cold_f-total_cold_cf)/total_cold_cf))

#% decrease in cold-related due to climate change for by country

summary_countries_seperate_cold <- summary_countries %>% 
  select(country, sim_hot_f_mean, sim_cold_f_mean, sim_hot_cf_mean,sim_cold_cf_mean) %>% 
  left_join(neonat_data) %>%
  left_join(births_data) %>% 
  mutate(total_heat_f = (sim_hot_f_mean/100)*total_neonat_deaths,
         total_cold_f = (sim_cold_f_mean/100)*total_neonat_deaths,
         total_heat_cf = (sim_hot_cf_mean/100)*total_neonat_deaths,
         total_cold_cf = (sim_cold_cf_mean/100)*total_neonat_deaths) %>% 
  mutate(cold_cc_increase =((total_cold_f-total_cold_cf)/total_cold_cf)) %>% 
  write.csv(here("Results", "Tables_SM", "Increase_cold-related_neonat_country.csv"), 
            row.names=FALSE)



######### 3. Share temp-attributable mortality for all countries #########

library("readxl")
#Calculate mean and upper and lower bounds for share of heat- and cold-related deaths
#in total neonatal mortality for all countries altogether
#Upload country-specific neonatal mortality data

summary_world <- summary_countries %>% 
  left_join(neonat_data) %>% 
  select(country, sim_hot_f_mean,sim_cold_f_mean,sim_hot_f_upp,sim_cold_f_upp,
         sim_hot_f_low, sim_cold_f_low, total_neonat_deaths) %>% 
  mutate(total_heat_f = (sim_hot_f_mean/100)*total_neonat_deaths,
         total_heat_f_upp = (sim_hot_f_upp/100)*total_neonat_deaths,
         total_heat_f_low = (sim_hot_f_low/100)*total_neonat_deaths,
         total_cold_f = (sim_cold_f_mean/100)*total_neonat_deaths,
         total_cold_f_upp = (sim_cold_f_upp/100)*total_neonat_deaths,
         total_cold_f_low = (sim_cold_f_low/100)*total_neonat_deaths) %>% 
  dplyr::select(country, total_heat_f,total_heat_f_upp,total_heat_f_low,
                total_cold_f, total_cold_f_upp, total_cold_f_low, total_neonat_deaths) %>% 
  summarise(total_heat_f = sum(total_heat_f), 
            total_heat_f_upp = sum(total_heat_f_upp),
            total_heat_f_low = sum(total_heat_f_low),
            total_cold_f = sum(total_cold_f),
            total_cold_f_upp = sum(total_cold_f_upp),
            total_cold_f_low = sum(total_cold_f_low),
            total_neonat_deaths=sum(total_neonat_deaths)) %>% 
  mutate(pc_heat = (total_heat_f/total_neonat_deaths)*100,
         pc_heat_upp = (total_heat_f_upp/total_neonat_deaths)*100,
         pc_heat_low = (total_heat_f_low/total_neonat_deaths)*100,
         pc_cold = (total_cold_f/total_neonat_deaths)*100,
         pc_cold_upp = (total_cold_f_upp/total_neonat_deaths)*100,
         pc_cold_low = (total_cold_f_low/total_neonat_deaths)*100,
         pc_both = ((total_heat_f+total_cold_f)/total_neonat_deaths)*100,
         pc_both_upp = ((total_heat_f_upp + total_cold_f_upp)/total_neonat_deaths)*100,
         pc_both_low = ((total_heat_f_low + total_cold_f_low)/total_neonat_deaths)*100)


###### Burden of climate change with respect to total heat/cold related neonatal deaths


######### 4. Share temp-related mortality due to CC by country #########
#Different results depending if I use this method (first, taking a mean
#across all observations and then doing the calculations or first doing the
#calculations and then taking a mean.
prop_temp <- summary %>% 
  group_by(country) %>% 
  summarise_at(vars(sim_hot_f, sim_hot_cf, sim_cold_f, sim_cold_cf),
               list(mean=mean)) %>% 
  mutate(prop_hot  = ((sim_hot_f_mean - sim_hot_cf_mean)/sim_hot_f_mean),
         prop_cold = ((sim_cold_cf_mean - sim_cold_f_mean)/sim_cold_cf_mean)) %>% 
  mutate_if(is.numeric, ~ . * 100) %>% 
  select(country, prop_hot, prop_cold) %>% 
  mutate_if(is.numeric, round, digits=0) 


######### 5. Share temp-related mortality due to CC for all countries #########
summary_cc_world <- summary_countries %>% 
  select(country, sim_hot_f_mean, sim_cold_cf_mean) %>% 
  left_join(neonat_data) %>%
  mutate(total_heat = (sim_hot_f_mean/100)*total_neonat_deaths, #Total deaths due to heat in fac
         total_cold = (sim_cold_cf_mean/100)*total_neonat_deaths) %>%  #Total deaths due to cold in count
  dplyr::select(country,total_heat, total_cold) %>% 
  left_join(prop_temp) %>% 
  mutate(total_heat_cc = total_heat*prop_hot/100, #Total heat-related due to climate change in fac (increase)
         total_cold_cc = total_cold*prop_cold/100) #Total cold related due to climate change in counter (reduction)
  
  
#Share of all heat-related neonat deaths due to cc (%)
summary_cc_world <- summary_cc_world %>% 
    summarise(total_heat = sum(total_heat), 
              total_heat_cc = sum(total_heat_cc),
              total_cold = sum(total_cold), 
              total_cold_cc = sum(total_cold_cc)) %>% 
    mutate(p_h_cc= (total_heat_cc/total_heat)*100,
           p_c_cc= (total_cold_cc/total_cold)*100)


#######      Fig. Share temp-attributable mortality due to CC in total mortality by country  #######



library("readxl")
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 
summary_ci <- summary %>%  
  mutate(sim_hot  = sim_hot_f - sim_hot_cf,
         sim_cold = sim_cold_f - sim_cold_cf) %>% 
  select(country, sim_hot,sim_cold) %>% 
  group_by(country) %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),  
                    low =~quantile(., probs = 0.025))) %>%   
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %



df_h <- left_join(summary_ci, continent)
# #Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "af_h_cc_all.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(df_h, aes(x= fct_reorder(country, sim_hot_mean), y= sim_hot_mean)) + 
  
  geom_point(colour="#a32d12", size=1.7) +
  geom_pointrange(aes(ymin = sim_hot_low, ymax = sim_hot_upp), size=0.6, colour="#a32d12" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  labs(x ="", y = "Heat-related neonatal deaths attributable to \nclimate change (percentage)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))
dev.off()

df_c <- left_join(summary_ci, continent)

jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "af_c_cc_all.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(df_c, aes(x= fct_reorder(country, sim_cold_mean), y= sim_cold_mean)) + 
  
  geom_pointrange(aes(ymin = sim_cold_low, ymax = sim_cold_upp), size=0.6, colour="#3454bf" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  
  labs(x ="", y = "Cold-related neonatal deaths attributable to \nclimate change (percentage)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))

dev.off()


#######      Fig. Share temp-attributable mortality in fac-countery by country  #######

 means_hot <- summary_countries %>% 
   select(country, sim_hot_f_mean, sim_hot_cf_mean) %>% 
   gather(scen, mean, sim_hot_f_mean: sim_hot_cf_mean, factor_key=TRUE) %>% 
   mutate(scen = ifelse(scen == "sim_hot_cf_mean", "Counterfactual" , "Factual"))
# 
 library("readxl")
 continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 
 
 means_hot <- left_join(means_hot, continent)
 means_hot <- spread(means_hot, scen, mean)
 library(forcats)
 
# #Order countries based on the results of the attribution part
 means_hot$country <- factor(means_hot$country,                              
                              levels = c("Uganda", "Liberia", "Sierra Leone",
                                         "Ethiopia", "Rwanda", "Cameroon",
                                         "Benin", "Guinea", "Burundi", "Senegal",
                                         "Mali", "Angola", "Nigeria", "Tanzania", 
                                         "South Africa", "Zimbabwe", "Malawi", "Zambia",
                                         "Philippines", "Timor-Leste", "Jordan", "Nepal",
                                         "Bangladesh", "Tajikistan", "Pakistan", "India",
                                         "Armenia","Haiti", "Albania"))
means_hot$country <-fct_rev(means_hot$country)
 
 
 jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "attr_farc_heat_scen.jpg"),
      units = 'in', res = 300, width = 4, height = 6)
 
 means_hot %>%
   ggplot(aes(x = Counterfactual, y = country)) +
   geom_col(aes(x = Factual), fill="#d9716c", alpha = 0.5, width = 0.8) +
   geom_col(width = 0.4, show.legend= FALSE, fill="#a3524e") +
   
   facet_grid(continent ~ ., scales = "free", space = "free", ) +
   labs(x ="Percentage", y = "")+
   theme(axis.text.x= element_text(face= "bold"))+
   geom_vline(xintercept=0)+
   theme_minimal()
 dev.off()
 
#Neonatal deaths ATTRIBUTABLE TO COLD - FACTUAL AND COUTERFACTUAL
means_cold <- summary_countries %>% 
   select(country, sim_cold_f_mean, sim_cold_cf_mean) %>% 
   gather(scen, mean, sim_cold_f_mean: sim_cold_cf_mean, factor_key=TRUE) %>% 
   mutate(scen = ifelse(scen == "sim_cold_cf_mean", "Counterfactual" , "Factual"))

library("readxl")
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 
means_cold <- left_join(means_cold, continent)
means_cold <- spread(means_cold, scen, mean)
library(forcats)
 
#Order countries based on the results of the attribution part
means_cold$country <- factor(means_cold$country,                                   
                              levels = c("Zimbabwe", "Zambia", "Malawi", "South Africa", "Angola",  "Mali", "Nigeria", "Senegal",
                                         "Tanzania",  "Burundi", "Sierra Leone", "Benin",  "Guinea", "Ethiopia", "Cameroon",
                                         "Liberia", "Rwanda", "Uganda", "Tajikistan", "Nepal", "India",  "Bangladesh", "Pakistan",
                                         "Armenia", "Jordan", "Timor-Leste", "Philippines",  "Haiti", "Albania"))
 
means_cold$country <-fct_rev(means_cold$country)

 
 
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "attr_farc_cold_scen.jpg"),
      units = 'in', res = 300, width = 4, height = 6)
 
 means_cold %>%
   ggplot(aes(x = Counterfactual, y = country)) +
   geom_col(aes(x = Factual), fill="#2b50cc", alpha = 0.5, width = 0.8) +
   geom_col(width = 0.4, show.legend= FALSE, fill="#3454bf") +
   
   facet_grid(continent ~ ., scales = "free", space = "free", ) +
   labs(x ="Percentage", y = "")+
   theme(axis.text.x= element_text(face= "bold"))+
   geom_vline(xintercept=0)+
   theme_minimal()
 dev.off()


#######   6. RATES (per 100,000) ####### 

library(reshape2)
library(here)
library(tidyverse)
options(scipen=999)
here::here()
library("readxl")
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 


neonat_data <- read_excel(here("Data/Neonatal_deaths_data_unicef.xlsx")) 

#Calculate births by country for the whole period (2001-2019)
births_data <- read_excel(here("Data/Births_data_unicef.xlsx")) %>% 
  dplyr::select(country, Year, births_total) %>% 
  dplyr::mutate(births_total= as.numeric(births_total)) %>% 
  group_by(country) %>% 
  summarise(births_total = sum(births_total))
#Calculate neonatal deaths by country for the whole period (2001-2019)
neonat_data <- neonat_data %>% 
  dplyr::select(country, Year, total_neonat_deaths) %>% 
  dplyr::mutate(total_neonat_deaths= as.numeric(total_neonat_deaths)) %>% 
  group_by(country) %>% 
  summarise(total_neonat_deaths = sum(total_neonat_deaths))



#########  6.1 Baseline neonatal mortality rate ##########
neonat_rate <- births_data %>% 
  left_join(neonat_data) %>% 
  mutate(rate = (total_neonat_deaths/(births_total*1000))*100000) %>% 
  dplyr::filter(country %in% c("Burundi", "Ethiopia", "Liberia", "Sierra Leone", "Rwanda",
                               "Uganda", "Cameroon", "Benin", "Guinea", "Angola", "Nigeria",
                               "Senegal", "Tanzania", "Mali", "Zimbabwe", "South Africa",
                               "Malawi", "Zambia", "Philippines", "Timor-Leste", "Bangladesh",
                               "Pakistan", "Nepal", "India","Jordan", "Tajikistan", "Armenia",
                               "Haiti", "Albania")) %>% 
  left_join(continent) 

#Order countries based on the results of the heat attribution part
neonat_rate$country <- factor(neonat_rate$country,                                    
             levels = c("Sierra Leone", "Ethiopia", "Liberia", "Mali", "Guinea",
                   "Benin", "Cameroon", "Nigeria", "Angola", "Uganda", "Rwanda", 
               "Burundi", "Senegal", "Tanzania", "Zimbabwe", "Malawi", "Zambia",
         "South Africa", "Timor-Leste", "Philippines", "Pakistan", "Bangladesh",
       "Nepal", "India", "Tajikistan", "Jordan", "Armenia", "Haiti", "Albania"))

neonat_rate$country <-fct_rev(neonat_rate$country)

jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "neonat_rate.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(neonat_rate, aes(x=country, y=rate)) +
  geom_bar(stat="identity",fill="#FF9900", colour="black") +  
  facet_grid(continent ~ ., scales = "free", space = "free") +
  labs(fill = "Category")+
  labs(x ="Country", y = "Neonatal mortality rate (per 100,000)")+  
  theme_minimal()+

  coord_flip()
dev.off()      


#########  6.2 Heat/cold-related neonatal mortality rate due to climate change  #########
library("readxl")
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 
summary_ci <- summary %>%  
  mutate(sim_hot  = sim_hot_f - sim_hot_cf,
         sim_cold = sim_cold_f - sim_cold_cf) %>% 
  select(country, sim_hot,sim_cold) %>% 
  group_by(country) %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975), 
                    low =~quantile(., probs = 0.025))) %>%   
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %



df_h <- left_join(summary_ci, continent)
df_h2 <- df_h %>% 
  dplyr::select(country, sim_hot_mean, sim_hot_upp, sim_hot_low, continent) %>% 
  mutate(sim_hot_mean=sim_hot_mean/100,
         sim_hot_upp= sim_hot_upp/100,
         sim_hot_low=sim_hot_low/100
  )
#Calculate heat-related neonatal mortality rate due to climate change
#Neonatal deaths are in total;births need to be multiplied by 1000
heat_total_cc <- df_h2 %>% 
  left_join(neonat_data) %>%
  left_join(births_data) %>% 
  mutate(total_mean = round(sim_hot_mean*total_neonat_deaths, digits = 0),
         total_upp  = round(sim_hot_upp*total_neonat_deaths, digits = 0),
         total_low  = round(sim_hot_low*total_neonat_deaths, digits = 0),
         total_ci = paste0(total_mean, " (", total_low, ", ", total_upp, ")"),
         rate_mean = round((total_mean/(births_total*1000))*100000, digits = 0),
         rate_upp  = round((total_upp/(births_total*1000))*100000, digits = 0),
         rate_low  = round((total_low /(births_total*1000))*100000, digits = 0),
         rate_ci = paste0(rate_mean, " (", rate_low, ", ", rate_upp, ")"))
write.csv(heat_total_cc, here("Results", "Tables_SM", "Heat_neonatal_CC_revision.csv"), row.names=FALSE)
#Figure - heat-related neonatal deaths rate attributable to climate change

jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "heat_cc_rate.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(heat_total_cc, aes(x= fct_reorder(country, rate_mean), y= rate_mean)) + 
  geom_pointrange(aes(ymin = rate_low, ymax = rate_upp), size=0.6, colour="#a32d12" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  
  labs( x ="", y = "Change in heat-related neonatal mortality rate\nattributable to climate change (per 100,000)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))
dev.off()


df_c <- left_join(summary_ci, continent) %>% 
  dplyr::select(country, sim_cold_mean, sim_cold_upp, sim_cold_low, continent) %>% 
  mutate(sim_cold_mean =sim_cold_mean/100,
         sim_cold_upp=sim_cold_upp/100,
         sim_cold_low=sim_cold_low/100 )
#Calculate heat-related neonatal mortality rate due to climate change
#Neonatal deaths are in total;births need to be multiplied by 1000
cold_total_cc <- df_c %>% 
  left_join(neonat_data) %>%
  left_join(births_data) %>% 
  mutate(total_mean = round(sim_cold_mean*total_neonat_deaths, digits = 0),
         total_upp  = round(sim_cold_upp*total_neonat_deaths, digits = 0),
         total_low  = round(sim_cold_low*total_neonat_deaths, digits = 0),
         total_ci = paste0(total_mean, " (", total_low, ", ", total_upp, ")"),
         rate_mean = round((total_mean/(births_total*1000))*100000, digits = 0),
         rate_upp  = round((total_upp/(births_total*1000))*100000, digits = 0),
         rate_low  = round((total_low /(births_total*1000))*100000, digits = 0),
         rate_ci = paste0(rate_mean, " (", rate_low, ", ", rate_upp, ")"))  
write.csv(cold_total_cc, here("Results", "Tables_SM", "Cold_neonat_CC_revision.csv"), row.names=FALSE)
#Order countries based on the results of the heat attribution part
cold_total_cc$country <- factor(cold_total_cc$country,                                   
           levels = c("Sierra Leone", "Ethiopia", "Liberia", "Mali", "Guinea",
                 "Benin", "Cameroon", "Nigeria", "Angola", "Uganda", "Rwanda", 
             "Burundi", "Senegal", "Tanzania", "Zimbabwe", "Malawi", "Zambia",
       "South Africa", "Timor-Leste", "Philippines", "Pakistan", "Bangladesh",
      "Nepal", "India", "Tajikistan", "Jordan", "Armenia", "Haiti", "Albania"))


cold_total_cc$country <-fct_rev(cold_total_cc$country)
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "cold_cc_rate.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(cold_total_cc, aes(x= country, y= rate_mean)) + 
  geom_pointrange(aes(ymin = rate_low, ymax = rate_upp), size=0.6, colour="#3454bf" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  
  labs(x ="", y = "Change in cold-related neonatal mortality rate\n attributable to climate change (per 100,000)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))
dev.off()

#########  6.3 Temp-attributable mortality by 6 categories #########
library(reshape2)
library(here)
library(tidyverse)
options(scipen=999)

here::here()

extr_h1 <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                         "obsclim_gswp3-w5e5_extr_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))

extr_h2 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-era5", 
                         "obsclim_20crv3-era5_extr_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))

extr_h3 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-w5e5", 
                         "obsclim_20crv3-w5e5_extr_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))



mod_h1 <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                        "obsclim_gswp3-w5e5_mod_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))

mod_h2 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-era5", 
                        "obsclim_20crv3-era5_mod_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))

mod_h3 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-w5e5", 
                        "obsclim_20crv3-w5e5_mod_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))


mild_h1 <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", 
                         "obsclim_gswp3-w5e5_mild_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))

mild_h2 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-era5", 
                         "obsclim_20crv3-era5_mild_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))

mild_h3 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-w5e5", 
                         "obsclim_20crv3-w5e5_mild_hot.csv")) %>% 
  dplyr::select(c("country", "attr_f_hot","range", "scen", "sim"))

hot_attr <- bind_rows(extr_h1, extr_h2, extr_h3,
                      mild_h1, mild_h2, mild_h3,
                      mod_h1, mod_h2, mod_h3) %>% 
  dplyr::rename("attr_f_total" = "attr_f_hot")

hot_attr <- hot_attr %>% 
  group_by(country, range) %>% 
  summarise(attr_f_total=mean(attr_f_total))

extr_c1 <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", "obsclim_gswp3-w5e5_extr_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))

extr_c2 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-era5", "obsclim_20crv3-era5_extr_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))

extr_c3 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-w5e5", "obsclim_20crv3-w5e5_extr_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))

mod_c1 <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", "obsclim_gswp3-w5e5_mod_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))
mod_c2 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-era5", "obsclim_20crv3-era5_mod_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))
mod_c3 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-w5e5", "obsclim_20crv3-w5e5_mod_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))


mild_c1 <- read.csv(here("Results", "Attributable_neonatal", "gswp3-w5e5", "obsclim_gswp3-w5e5_mild_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))
mild_c2 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-era5", "obsclim_20crv3-era5_mild_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))
mild_c3 <- read.csv(here("Results", "Attributable_neonatal", "20crv3-w5e5", "obsclim_20crv3-w5e5_mild_cold.csv")) %>% 
  dplyr::select(c("country", "attr_f_cold","range", "scen", "sim"))

cold_attr <- bind_rows(extr_c1, extr_c2, extr_c3,
                       mod_c1, mod_c2, mod_c3,
                       mild_c1, mild_c2, mild_c3) %>% 
  dplyr::rename("attr_f_total" = "attr_f_cold")

cold_attr <- cold_attr %>% 
  group_by(country, range) %>% 
  summarise(attr_f_total=mean(attr_f_total))


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
#Calculate neonatal deaths by country for the whole period (2001-2019)
neonat_data <- neonat_data %>% 
  dplyr::select(country, Year, total_neonat_deaths) %>% 
  dplyr::mutate(total_neonat_deaths= as.numeric(total_neonat_deaths)) %>% 
  group_by(country) %>% 
  summarise(total_neonat_deaths = sum(total_neonat_deaths))
#Neonatal deaths per 100,000 live births

neonat_data_relative <- neonat_data %>% 
  left_join(births_data) %>%
  mutate(neonat_relative = (total_neonat_deaths/births_total*1000)*100000)

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
         range=factor(range, c("extremely hot", "moderately hot", "mildly hot", "mildly cold",
                               "moderately cold", "extremely cold"))) 
str(total_temp)
total_temp <- attr_frac_1 %>% 
  dplyr::select(country, range, value) %>% 
  spread(range, value) 

attr_frac_1 <- left_join(attr_frac_1, continent) 

  

attr_frac_1$country <- factor(attr_frac_1$country,                                   
             levels = c("Sierra Leone", "Ethiopia", "Liberia", "Mali", "Guinea",
                   "Benin", "Cameroon", "Nigeria", "Angola", "Uganda", "Rwanda", 
               "Burundi", "Senegal", "Tanzania", "Zimbabwe", "Malawi", "Zambia",
         "South Africa", "Timor-Leste", "Philippines", "Pakistan", "Bangladesh",
       "Nepal", "India", "Tajikistan", "Jordan", "Armenia", "Haiti", "Albania"))


attr_frac_1$country <-fct_rev(attr_frac_1$country)
#######      Fig. Share temp-attributable mortality by category ####### 
#lot attributable fraction by country

jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "heat_cold_rate_6categories_nonames.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(attr_frac_1, aes(x=country, y=value, fill=range)) +
  geom_bar(position = "stack", stat="identity", color="black")+  
  scale_fill_manual(values=c('#e85b54', '#eb7c60', '#eba160','#9bd2f2',  '#4193f0', '#415ef0'))+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  labs(fill = "Temperature range")+
  labs(x ="Country", y = "Heat- and cold-related neonatal mortality rate (per 100,000)")+  #41.2
  theme_minimal()+
  theme(axis.title  = element_text(size = 10),)+
  coord_flip()

dev.off()

######### 7. Share temp-related mortality due to CC by country for each dataset #########

library(reshape2)
library(here)
library(tidyverse)
options(scipen=999)


here::here()
###  1) gswp3-w5e5

fac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                       "Factual_gswp3-w5e5.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)

countfac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                            "Counterfactual_gswp3-w5e5.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

summary <- fac_1 %>% 
  cbind(countfac_1) 
summary <- summary[, !duplicated(colnames(summary))] 
summary <- summary %>%  
  select(-scenario) %>% 
  mutate(sim_hot  = sim_hot_f- sim_hot_cf,
         sim_cold = sim_cold_f- sim_cold_cf) %>% 
  select(country, sim_hot,sim_cold)

summary_ci <- summary %>% 
  group_by(country) %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),
                    low =~quantile(., probs = 0.025))) %>% 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %

library("readxl")
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 


df_h <- left_join(summary_ci, continent)
df_h_1 <- df_h %>% 
  select(country, sim_hot_mean, continent) %>% 
  dplyr::rename("gswp3-w5e5" = sim_hot_mean)

# #Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "gswp3-w5e5", "af_h_cc_gswp3-w5e5.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(df_h, aes(x= fct_reorder(country, sim_hot_mean), y= sim_hot_mean)) + 
  geom_pointrange(aes(ymin = sim_hot_low, ymax = sim_hot_upp), size=0.6, colour="#a32d12" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  labs(title= "gswp3-w5e5", 
       x ="", y = "Heat-related neonatal deaths attributable to \nclimate change (percentage)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))

theme_minimal()
dev.off()


df_c <- left_join(summary_ci, continent)
df_c_1 <- df_c %>% 
  select(country, sim_cold_mean, continent) %>% 
  dplyr::rename("gswp3-w5e5" = sim_cold_mean)


#Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "gswp3-w5e5", "af_c_cc_gswp3-w5e5.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(df_c, aes(x= fct_reorder(country, sim_cold_mean), y= sim_cold_mean)) + 
  
  geom_pointrange(aes(ymin = sim_cold_low, ymax = sim_cold_upp), size=0.6, colour="#3454bf" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  
  labs(x ="", y = "Cold-related neonatal deaths attributable to \nclimate change (percentage)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))
dev.off()

###  2) 20crv3-era5

fac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", "Factual_20crv3-era5.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)

countfac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", "Counterfactual_20crv3-era5.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

summary <- fac_1 %>% 
  cbind(countfac_1) 
summary <- summary[, !duplicated(colnames(summary))] 
summary <- summary %>%  
  select(-scenario) %>% 
  mutate(sim_hot  = sim_hot_f- sim_hot_cf,
         sim_cold = sim_cold_f- sim_cold_cf) %>% 
  select(country, sim_hot,sim_cold)

summary_ci <- summary %>% 
  group_by(country) %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),
                    low =~quantile(., probs = 0.025))) %>% 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %

library("readxl")
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 


df_h <- left_join(summary_ci, continent)
df_h_2 <- df_h %>% 
  select(country, sim_hot_mean, continent) %>%
  dplyr::rename("20crv3-era5" = sim_hot_mean)

# #Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "20crv3-era5", "af_h_cc_20crv3-era5.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(df_h, aes(x= fct_reorder(country, sim_hot_mean), y= sim_hot_mean)) + 
  geom_pointrange(aes(ymin = sim_hot_low, ymax = sim_hot_upp), size=0.6, colour="#a32d12" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  labs(title= "20crv3-era5", x ="", y = "Heat-related neonatal deaths attributable
       to \nclimate change (percentage)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))
dev.off()





df_c <- left_join(summary_ci, continent)
df_c_2 <- df_c %>% 
  select(country, sim_cold_mean, continent) %>% 
  dplyr::rename("20crv3-era5" = sim_cold_mean)

#Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "20crv3-era5", 
    "af_c_cc_20crv3-era5.jpg"), units = 'in', res = 300, width = 4, height = 6)

ggplot(df_c, aes(x= fct_reorder(country, sim_cold_mean), y= sim_cold_mean)) + 
  
  geom_pointrange(aes(ymin = sim_cold_low, ymax = sim_cold_upp), size=0.6, colour="#3454bf" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  labs(title= "20crv3-era5", 
       x ="", y = "Cold-related neonatal deaths attributable to \nclimate change (percentage)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))
dev.off()

###  3) 20crv3-w5e5

fac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                       "Factual_20crv3-w5e5.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)

countfac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                            "Counterfactual_20crv3-w5e5.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

summary <- fac_1 %>% 
  cbind(countfac_1) 
summary <- summary[, !duplicated(colnames(summary))] 
summary <- summary %>%  
  select(-scenario) %>% 
  mutate(sim_hot  = sim_hot_f- sim_hot_cf,
         sim_cold = sim_cold_f- sim_cold_cf) %>% 
  select(country, sim_hot,sim_cold)

summary_ci <- summary %>% 
  group_by(country) %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),
                    low =~quantile(., probs = 0.025))) %>% 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %

library("readxl")
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) 


df_h <- left_join(summary_ci, continent)
df_h_3 <- df_h %>% 
  select(country, sim_hot_mean, continent) %>% 
  dplyr::rename("20crv3-w5e5" = sim_hot_mean)


#Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "20crv3-w5e5", "attr_farc_cc_20crv3-w5e5_mod_extr.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(df_h, aes(x= fct_reorder(country, sim_hot_mean), y= sim_hot_mean)) + 

  geom_pointrange(aes(ymin = sim_hot_low, ymax = sim_hot_upp), size=0.6, colour="#a32d12" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  labs(title= "20crv3-w5e5",
       x ="", y = "Heat-related neonatal deaths attributable to \nclimate change (percentage)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))

dev.off()





df_c <- left_join(summary_ci, continent)
df_c_3 <- df_c %>% 
  select(country, sim_cold_mean, continent) %>% 
  dplyr::rename("20crv3-w5e5" = sim_cold_mean)

#Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "20crv3-w5e5", "af_c_cc_20crv3-w5e5.jpg"),
     units = 'in', res = 300, width = 4, height = 6)

ggplot(df_c, aes(x= fct_reorder(country, sim_cold_mean), y= sim_cold_mean)) + 
  
  geom_pointrange(aes(ymin = sim_cold_low, ymax = sim_cold_upp), size=0.6, colour="#3454bf" )+
  geom_hline(yintercept=0)+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  coord_flip()+
  labs(title= "20crv3-w5e5", x ="", y = "Cold-related neonatal deaths attributable 
       to \nclimate change (percentage)")+
  theme_minimal()+
  theme(axis.title  = element_text(size = 10))
dev.off()


######### 8. Dotplot - dataset comparison in share temp-related mortality due to CC by #########


hot_datasets <- df_h_1 %>% 
  left_join(df_h_2) %>% 
  left_join(df_h_3) %>% 
  gather(key = "Dataset", value = "Percentage", "gswp3-w5e5", "20crv3-era5", 
         "20crv3-w5e5")


cold_datasets <- df_c_1 %>% 
  left_join(df_c_2) %>% 
  left_join(df_c_3) %>% 
  gather(key = "Dataset", value = "Percentage", "gswp3-w5e5", "20crv3-era5", 
         "20crv3-w5e5")


theme_dotplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(.75)),
        axis.ticks.y = element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_line(size = 0.5),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())

jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", 
   "heat_datasets_updated.jpg"), units = 'in', res = 300, width = 5, height = 6)

ggplot(hot_datasets, aes(Percentage, country, color = Dataset)) +
  geom_point(size=4) +  
  labs(x= "Heat-related neonatal mortality attributed to climate change (%)")+
  geom_vline(xintercept=0, colour="dark grey")+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  theme_minimal() +
  scale_color_manual(values = c("#bb0032", "#E7B800", "#FC4E07"))
dev.off() 


jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "cold_datasets_updated.jpg"),
     units = 'in', res = 300, width = 5, height = 6)

ggplot(cold_datasets, aes(Percentage, country, color = Dataset)) +
  geom_point(size=4) +  
  labs(x= "Cold-related neonatal mortality attributed to climate change (%)")+
  geom_vline(xintercept=0, colour="dark grey")+
  facet_grid(continent ~ ., scales = "free", space = "free") +
  theme_minimal() +
  scale_color_manual(values = c("#4b00bb", "#00bbbb","#005abb"))
dev.off() 


######### 9. Comparison across datasets ######################
library(reshape2)
library(here)
library(tidyverse)
options(scipen=999)
here::here()

###  1) gswp3-w5e5

fac_1 <- read.csv(here("Results", "Attributable_neonat", "Monte_Carlo", 
                       "Factual_gswp3-w5e5_all.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)

countfac_1 <- read.csv(here("Results", "Attributable_neonat", "Monte_Carlo",
                            "Counterfactual_gswp3-w5e5_all.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

summary1 <- fac_1 %>% 
  cbind(countfac_1) 
summary1 <- summary1[, !duplicated(colnames(summary1))] 
summary1 <- summary1 %>%  
  select(-scenario) %>% 
  mutate(sim_hot  = sim_hot_f- sim_hot_cf,
         sim_cold = sim_cold_f- sim_cold_cf) %>% 
  select(sim_hot,sim_cold)

summary_ci1 <- summary1 %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),
                    low =~quantile(., probs = 0.025))) %>% 
  mutate(dataset= "gswp3-w5e5") %>% 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %


###  1) 20crv3-era5

fac_2 <- read.csv(here("Results", "Attributable_neonat", "Monte_Carlo", 
                       "Factual_20crv3-era5_all.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)

countfac_2 <- read.csv(here("Results", "Attributable_neonat", "Monte_Carlo", 
                            "Counterfactual_20crv3-era5_all.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

summary2 <- fac_2 %>% 
  cbind(countfac_2) 
summary2 <- summary2[, !duplicated(colnames(summary2))] 
summary2 <- summary2 %>%  
  select(-scenario) %>% 
  mutate(sim_hot  = sim_hot_f- sim_hot_cf,
         sim_cold = sim_cold_f- sim_cold_cf) %>% 
  select(sim_hot,sim_cold)

summary_ci2 <- summary2 %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),
                    low =~quantile(., probs = 0.025))) %>% 
  mutate(dataset= "20crv3-era5") %>% 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %


###  1) 20crv3-w5e5

fac_3 <- read.csv(here("Results", "Attributable_neonat", "Monte_Carlo", 
                       "Factual_20crv3-w5e5_all.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)

countfac_3 <- read.csv(here("Results", "Attributable_neonat", "Monte_Carlo", 
                            "Counterfactual_20crv3-w5e5_all.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

summary3 <- fac_3 %>% 
  cbind(countfac_3) 
summary3 <- summary3[, !duplicated(colnames(summary3))] 
summary3 <- summary3 %>%  
  select(-scenario) %>% 
  mutate(sim_hot  = sim_hot_f- sim_hot_cf,
         sim_cold = sim_cold_f- sim_cold_cf) %>% 
  select(sim_hot,sim_cold)

summary_ci3 <- summary3 %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),
                    low =~quantile(., probs = 0.025))) %>% 
  mutate(dataset= "20crv3-w5e5") %>% 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %


all_df <- summary_ci1 %>% 
  rbind(summary_ci2) %>% 
  rbind(summary_ci3)
# #Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "heat_allcountr.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

ggplot(all_df, aes(x= fct_reorder(dataset, sim_hot_mean), y= sim_hot_mean)) + 
  geom_pointrange(aes(ymin = sim_hot_low, ymax = sim_hot_upp), size=0.6, colour="#a32d12" )+
  geom_hline(yintercept=0)+
  labs(title= "Heat impacts across countries by dataset", 
       x ="", y = "Percentage")+
  theme_minimal()

dev.off()


#Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "cold_allcountr.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

ggplot(all_df, aes(x=  fct_reorder(dataset, sim_cold_mean), y= sim_cold_mean)) + 
  
  geom_pointrange(aes(ymin = sim_cold_low, ymax = sim_cold_upp), size=0.6, colour="#3454bf" )+
  geom_hline(yintercept=0)+
  labs(title= "Cold impacts across countries by dataset",   
       x ="", y = "Percentage")+
  theme_minimal()
dev.off()


#########  10. Comparison across datasets - world ######################

library(reshape2)
library(here)
library(tidyverse)
options(scipen=999)
here::here()

###  1) gswp3-w5e5

fac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                      "Factual_gswp3-w5e5_all_countries.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)

countfac_1 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                            "Counterfactual_gswp3-w5e5_all_countries.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

summary1 <- fac_1 %>% 
  cbind(countfac_1) 
summary1 <- summary1[, !duplicated(colnames(summary1))] 
summary1 <- summary1 %>%  
  select(-scenario) %>% 
  mutate(sim_hot  = sim_hot_f- sim_hot_cf,
         sim_cold = sim_cold_f- sim_cold_cf) %>% 
  select(sim_hot,sim_cold)

summary_ci1 <- summary1 %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),
                    low =~quantile(., probs = 0.025))) %>% 
  mutate(dataset= "gswp3-w5e5") %>% 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %



###  1) 20crv3-era5

fac_2 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                       "Factual_20crv3-era5_all_countries.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)

countfac_2 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                            "Counterfactual_20crv3-era5_all_countries.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

summary2 <- fac_2 %>% 
  cbind(countfac_2) 
summary2 <- summary2[, !duplicated(colnames(summary2))] 
summary2 <- summary2 %>%  
  select(-scenario) %>% 
  mutate(sim_hot  = sim_hot_f- sim_hot_cf,
         sim_cold = sim_cold_f- sim_cold_cf) %>% 
  select(sim_hot,sim_cold)

summary_ci2 <- summary2 %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),
                    low =~quantile(., probs = 0.025))) %>% 
  mutate(dataset= "20crv3-era5") %>% 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %


###  1) 20crv3-w5e5

fac_3 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                       "Factual_20crv3-w5e5_all_countries.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_f = sim_hot,
                sim_cold_f = sim_cold)

countfac_3 <- read.csv(here("Results", "Attributable_neonatal", "Monte_Carlo", 
                            "Counterfactual_20crv3-w5e5_all_countries.csv"))[,-1] %>% 
  dplyr::rename(sim_hot_cf = sim_hot,
                sim_cold_cf = sim_cold)

summary3 <- fac_3 %>% 
  cbind(countfac_3) 
summary3 <- summary3[, !duplicated(colnames(summary3))] 
summary3 <- summary3 %>%  
  select(-scenario) %>% 
  mutate(sim_hot  = sim_hot_f- sim_hot_cf,
         sim_cold = sim_cold_f- sim_cold_cf) %>% 
  select(sim_hot,sim_cold)

summary_ci3 <- summary3 %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),
                    low =~quantile(., probs = 0.025))) %>% 
  mutate(dataset= "20crv3-w5e5") %>% 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %


all_df <- summary_ci1 %>% 
  rbind(summary_ci2) %>% 
  rbind(summary_ci3)
#Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "heat_allcountries.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

ggplot(all_df, aes(x= fct_reorder(dataset, sim_hot_mean), y= sim_hot_mean)) + 
  geom_pointrange(aes(ymin = sim_hot_low, ymax = sim_hot_upp), size=0.6, colour="#a32d12" )+
  geom_hline(yintercept=0)+
  labs(title= "Heat impacts across all countries by dataset", 
       x ="", y = "Percentage")+
  theme_minimal()

dev.off()


# #Plot attributable fraction by country
jpeg(here("Figures", "Neonatal mortality" , "Stage II", "all_datasets", "cold_allcountr.jpg"),
     units = 'in', res = 300, width = 10, height = 6)

ggplot(all_df, aes(x=  fct_reorder(dataset, sim_cold_mean), y= sim_cold_mean)) + 
  
  geom_pointrange(aes(ymin = sim_cold_low, ymax = sim_cold_upp), size=0.6, colour="#3454bf" )+
  geom_hline(yintercept=0)+
  labs(title= "Cold impacts across all countries by dataset",  
       x ="", y = "Percentage")+
  theme_minimal()
dev.off()
###  All 3 datasets combined for all counties

fac <- rbind.data.frame(fac_1, fac_2, fac_3) 
countfac <- rbind.data.frame(countfac_1, countfac_2, countfac_3) 
summary <- fac %>% 
  cbind(countfac) 
summary <- summary[, !duplicated(colnames(summary))]  
summary <- summary %>%  
  select(-scenario) 

#Calclulate as a proportion of all heat/cold-related neonatal deaths
prop_temp <- summary %>% 
  summarise_at(vars(sim_hot_f, sim_hot_cf, sim_cold_f, sim_cold_cf),
               list(mean=mean)) %>% 
  mutate(prop_hot  = ((sim_hot_f_mean - sim_hot_cf_mean)/sim_hot_f_mean),
         prop_cold = ((sim_cold_f_mean - sim_cold_cf_mean)/sim_cold_cf_mean)) %>% 
  mutate_if(is.numeric, ~ . * 100) %>% 
  select(prop_hot, prop_cold)
#Share of heat- and cold-related neonatal deaths in factual scenario
summary_countries <- summary %>% 
  summarise_at(vars(sim_hot_f, sim_hot_cf, sim_cold_f, sim_cold_cf),
               list(mean=mean, upp =~quantile(., probs = 0.95),
                    low =~quantile(., probs = 0.05))) %>% #, 
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %



#Share of total neonatal deaths associated to non-optimal temperatures
summary_countries2 <- summary %>% 
  dplyr::select(sim_hot_f, sim_hot_cf, sim_cold_f, sim_cold_cf) %>% 
  mutate(factual= sim_hot_f+ sim_cold_f,
         counterfactual = sim_hot_cf+ sim_cold_cf) %>% 
  summarise_at(vars(factual, counterfactual),
               list(mean=mean, upp =~quantile(., probs = 0.95),
                    low =~quantile(., probs = 0.05))) %>% #, 
  mutate_if(is.numeric, ~ . * 100)

summary <- summary %>%  
  mutate(sim_hot  = sim_hot_f - sim_hot_cf,
         sim_cold = sim_cold_f - sim_cold_cf) %>% 
  select(sim_hot,sim_cold)


#Share of heat- and cold-related neonatal deaths attributable to climate change
summary_ci <- summary %>% 
  summarise_at(vars(sim_hot, sim_cold),
               list(mean=mean, upp =~quantile(., probs = 0.975),  
                    low =~quantile(., probs = 0.025))) %>%   #
  mutate_if(is.numeric, ~ . * 100) #multiply by 100 to get the fraction in %

######### 11. Optimal temperature vs. Latitude and average temp    ###### 
library(reshape2)
library(here)
library(tidyverse)
library("readxl")
options(scipen=999)

#Upload data with latitude & continent
continent <- read_excel(here("Data/Data_1st stage/country_continent.xlsx")) %>% 
dplyr::rename(Country = country)
 
#Upload data with average annual temperature for each country 
#for the study period
#Upoad the daily temperature data


mean_temp <- read.csv(here("Data/Data_2nd stage/ISIMIP data/Annual_TMP_by_country.csv")) [,-1]
mean_temp <- mean_temp %>%
   dplyr::select(Country, year, tmp_obsclim) %>% 
   dplyr::filter(as.numeric(year) >=2001 & as.numeric(year) <=2019)

mean_temp_wide <- spread(mean_temp, year, tmp_obsclim) %>% 
   mutate(mean = rowMeans(across('2001': '2019'))) %>% 
   select(Country, mean) %>% 
   mutate(mean=round(mean,0))

 country_identifier <- read.csv(here("Data/Data_1st stage/country_identified.csv")) [,-1]

 
 MMT_1 <- read.csv(here("Data/Data_1st stage/MMT_gswp3-w5e5.csv")) [,-1] %>% 
   dplyr::rename(mean_pct_1 = mean_pct)
 MMT_2 <- read.csv(here("Data/Data_1st stage/MMT_20crv3-era5.csv")) [,-1] %>% 
   dplyr::rename(mean_pct_2 = mean_pct)

 MMT <- MMT_1 %>% 
   left_join(MMT_2) %>% 
   rowwise() %>%
   mutate(mean_pct = mean(c_across(c('mean_pct_1', 'mean_pct_2')), na.rm=TRUE)) %>% 
   ungroup() %>% 
   dplyr::select(SurveyId, tmean_pct,mean_pct)
 

all <- country_identifier %>% 
   left_join(MMT) %>% 
   filter(tmean_pct== 53) %>% 
   mutate(mean_pct=round(mean_pct,0)) %>% 
   dplyr::rename(Country = CountryName,
                 MMT = mean_pct) %>% #Average of 51.9 and 53.6 filter only the ref temp pct 
   left_join(mean_temp_wide) %>%
   left_join(continent) %>% 
   mutate(latitude=as.numeric(latitude),
          latitude=round(latitude, digits=0)) %>% 
   dplyr::rename(Continent=continent)
 library(ggplot2)
 library(hrbrthemes)
 

