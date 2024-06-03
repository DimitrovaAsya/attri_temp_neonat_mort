
#################################################################################
#################################################################################
#######                                                                    ######
####### PROJECT: Temperature-related neonatal deaths attributable to       ######
#######          climate change in 29 low- and middle-income countries     ######                         
#######                                                                    ######
#######                                                                    ###### 
#######             CODE: BLUP method - Neonatal mortality                 ######
#################################################################################
#################################################################################
#################################################################################
# FIRST STAGE
# - DEFINE THE CROSS-BASIS MATRICES 
# - BUILD OBJECTS TO STORE THE RESULTS
# - RUN THE CONDITIONAL LOGISTIC MODELS
# - REDUCE THE FITTED MAIN MODEL TO SUMMARIES
# - STORE THE RESULTS
####################################################################

library(here)
library(survival)
library(dlnm)
library(tidyverse)

here::here()

load(here("Data/Data_1st stage/Linked_data_ISIMIP/gswp3-w5e5/obsclim_gswp3-w5e5_2.RData"))

print(max(neonat_mort$tmean_lag0, na.rm=TRUE)-min(neonat_mort$tmean_lag0, na.rm=TRUE))
min(neonat_mort$tmean_lag0, na.rm=TRUE)
max(neonat_mort$tmean_lag0, na.rm=TRUE)


data <- neonat_mort  %>%
  mutate(year  = format(case_date, format = "%Y"),
         month = format(case_date, format = "%m")) %>% 
  dplyr::mutate(SurveyYear = substr(SurveyId, 3, 6)) %>% 
  mutate(SurveyYear = as.numeric(SurveyYear)) %>% 
  mutate(year = as.numeric(year)) %>% 
  filter(year_int - year <= 14) %>% 
  drop_na(c(tmean_lag0, tmean_lag1,tmean_lag2, tmean_lag3, 
            tmean_lag4, tmean_lag5, tmean_lag6, tmean_lag7))   

data_neonat <- data %>%   # Number of neonatal deaths observations
  summarise(neonat_mort = sum(neonat_mort))

data_neonat_country <- data %>%   
  group_by(CountryName) %>% 
  summarise(neonat_mort = sum(neonat_mort))

 data <- data %>% 
   dplyr::select(c(CountryName, neonat_mort, strata_id,tmean_lag0, tmean_lag1,
                   tmean_lag2)) 

# LOAD THE PACKAGE
library(dlnm) ; library(splines) ; library(xtable)

# CHECK VERSION OF THE PACKAGE
if(packageVersion("dlnm")<"2.2.0")
  stop("update dlnm package to version >= 2.2.0")


#COUNTRIES
countries <- as.character(unique(data$CountryName))

####################################################################
# CREATE A LIST OF DATAFRAMES FOR EACH OF THE 29 COUNTRIES' SERIES

datalist <- lapply(countries, function(country) data[data$CountryName==country,])
names(datalist) <- countries
m <- length(datalist)

# TEMPERATURE RANGES
ranges <- t(sapply(datalist, function(x) range(x$tmean_lag0,na.rm=T)))

lag <- c(0,2)
bound <- colMeans(ranges)
yall <- matrix(NA,length(datalist),2,dimnames=list(countries,paste("b",seq(2),sep="")))

Sall <- vector("list",length(datalist))
names(Sall) <- countries
Shot <- Scold <-Sall

####################################################################
arglag <- list(fun="ns",knots=logknots(2,nk=1))  #For lag 2 days

knots <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 98)


length(knots)
#Create a list that contains the lists that will store the models for each country
df_1 <- vector("list",length(datalist))
names(df_1) <- countries
#Create a list to store the models for each country 
df_2 <- vector("list",length(knots))
names(df_2) <- knots 

for(i in seq(datalist)) {
  # PRINT
  cat(i,"")
 
  # LOAD
  sub <- datalist[[i]]
  Qstillbirth <- sub  %>% 
    dplyr::select(tmean_lag0:tmean_lag2) 
  
  Qstillbirth <- as.matrix(Qstillbirth) 
  Qstillbirth <- unname(Qstillbirth)

  for (k in knots) {
    knotperc <- k
    varknots <- mean(sapply(datalist,function(x) quantile(x$tmean_lag0, knotperc/100))) 
    argvar=list(fun="ns", knots = varknots, bound = bound)
    arglag <- list(fun="ns",knots=logknots(2,nk=1))
    # DEFINE THE CROSS-BASES
    cb <- crossbasis(Qstillbirth,lag=lag,argvar=argvar,arglag=arglag)
    # RUN THE FIRST-STAGE MODELS
    model <- clogit(neonat_mort ~ cb + strata(strata_id), data=sub)
    
    df_2[[paste0(k)]] <- model

  }
  df_1[[i]] <- df_2
}
  
#Apply the aictab to each country list
library(AICcmodavg)
List <- lapply(df_1, aictab)
str(tables)   
#Calculates AICc for all models for each country and makes a table with ranking      
aictab(cand.set = models_ns, modnames = model.names_ns)
#Select a model that has the minimum value of the sum of the AIC in all the 25 countries
tables <- bind_rows(List, .id = "Country") %>% 
  dplyr::select(Country, Modnames, AICc)

tables_wide <- spread(tables, Country, AICc) %>% 
  mutate(Total = rowSums(tables_wide[ , 2:25]))
######--> 98 percentile shown as the best model

lag <- c(0,2)
bound <- colMeans(ranges)

knotperc <- 98
varknots <- mean(sapply(datalist,function(x) quantile(x$tmean_lag0, knotperc/100))) 
argvar=list(fun="ns", knots = varknots, bound = bound)
arglag <- list(fun="ns",knots=logknots(2,nk=1))

yall <- matrix(NA,length(datalist),2,dimnames=list(countries,paste("b",seq(2),sep="")))


Sall <- vector("list",length(datalist))
names(Sall) <- countries
Shot <- Scold <-Sall

####################################################################
arglag <- list(fun="ns",knots=logknots(2,nk=1))  #For lag 2 days


# RUN THE MODEL FOR EACH COUNTRY

# WARNING FOR PREDICTION BEYOND BOUNDARIES SUPPRESSED
options(warn=-1)

# LOOP FOR COUNTREIES
for(i in seq(datalist)) {
  
  # PRINT
  cat(i,"")
  str(sub)
  # LOAD
  sub <- datalist[[i]]
  Qstillbirth <- sub  %>% 
    dplyr::select(tmean_lag0:tmean_lag2) 
  
  Qstillbirth <- as.matrix(Qstillbirth) 
  Qstillbirth <- unname(Qstillbirth)
  
  # DEFINE THE CROSS-BASES
  cb <- crossbasis(Qstillbirth,lag=lag,argvar=argvar,arglag=arglag)
 
  # RUN THE FIRST-STAGE MODELS
  model <- clogit(neonat_mort ~ cb + strata(strata_id), data=sub)

  crall <- crossreduce(cb,model, cen= 26)
 
  # STORE THE RESULTS
  # EXTRACT AND SAVE THE RELATED COEF AND VCOV
  # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
  yall[i,] <- coef(crall)
  Sall[[i]] <- vcov(crall)

}



 ####################################################################
 # SECOND STAGE
 # - RUN THE MULTIVARIATE META-ANALYTICAL MODELS WITH mvmeta
 # - CREATE BASIS VARIABLES USING onebasis, TO BE USED FOR PREDICTION
 # - OBTAIN PREDICTIONS THROUGH crosspred (dlnm)
 ####################################################################
 
 ####################################################################
 # PERFORM MULTIVARIATE META-ANALYSIS
 
 # LOAD THE PACKAGES (mvmeta PACKAGE IS ASSUMED TO BE INSTALLED)
 library(mvmeta)
 
 # SELECT THE ESTIMATION METHOD
 method <- "reml"
 # IN THE CURRENT VERSION, SET control=list(showiter=T) TO 
 #   INSPECT THE OPTIMIZATION SEARCH
 
 # OVERALL CUMULATIVE SUMMARY FOR THE MAIN MODEL
 mvall <- mvmeta(yall~1,Sall,method=method)
 summary(mvall)
 

 ####################################################################
 # CREATE BASES FOR PREDICTION
 
 # BASES OF TEMPERATURE AND LAG USED TO PREDICT, EQUAL TO THAT USED FOR ESTIMATION
 # COMPUTED USING THE ATTRIBUTES OF THE CROSS-BASIS USED IN ESTIMATION
 xvar <- seq(bound[1],bound[2],by=0.1)
 bvar <- do.call("onebasis",c(list(x=xvar),attr(cb,"argvar")))
 xlag <- 0:20/10
 blag <- do.call("onebasis",c(list(x=xlag),attr(cb,"arglag")))
 
 ####################################################################
 # REGION-SPECIFIC FIRST-STAGE SUMMARIES
 
 regall <- lapply(seq(nrow(yall)),function(i) crosspred(bvar,coef=yall[i,],
                    vcov=Sall[[i]],model.link="log",cen=26))
 
 ####################################################################
 # PREDICTION FOR A GRID OF TEMPERATURE AND LAG VALUES
 
 # OVERALL CUMULATIVE SUMMARY ASSOCIATION FOR MAIN MODEL
 cpall <- crosspred(bvar,coef=coef(mvall),vcov=vcov(mvall),
                    model.link="log",by=0.1,from=bound[1],to=bound[2],cen=26)
 

 ###################################################################
 # OVERALL CUMULATIVE SUMMARY ASSOCIATION
 
 # PLOT

 jpeg(here("Figures","Neonatal mortality", "Stage I", "gswp3-w5e5", "pt_Lag0-2_gswp3-w5e5_BLUP.jpg"),
      units = 'in', res = 300, width = 10, height = 6)
 
 par(mar=c(5,4,1,1)+0.1,cex.axis=0.9,mgp=c(2.5,1,0))
 layout(matrix(1:1,ncol=1))
 
 plot(cpall,type="n",ylab="RR",ylim=c(.8,2),xlab="Temperature (C)")
 for(i in seq(regall)) lines(regall[[i]],ptype="overall",col=grey(0.5),lty=2)
 abline(h=1)
 lines(cpall,col=2,lwd=2)
 mtext("Main model: first-stage and pooled estimates",cex=1)
 legend ("top",c("Pooled (with 95%CI)","First-stage region-specific"),
         lty=c(1,2),lwd=1.5,col=c(2,grey(0.7)),bty="n",inset=0.1,cex=0.8)
 dev.off()
 
 # POINT OF MINIMUM MORTALITY
 cpall$predvar[which.min(cpall$allRRfit)]
 round(sum(regEngWales$tmean<17.1)/nrow(regEngWales)*100,1)
 
 # Q TEST AND I-SQUARE
 (qall <- qtest(mvall))
 round(((qall$Q-qall$df)/qall$Q)[1]*100,1)
 (qall2 <- qtest(mvall2))
 round(((qall2$Q-qall2$df)/qall2$Q)[1]*100,1)
 (qall3 <- qtest(mvall3))
 round(((qall3$Q-qall3$df)/qall3$Q)[1]*100,1)
 








