
# Load data ----


## Read catch at age, change units
catage = as.matrix(read.csv("catch_at_age_matrix.csv", row.names=1)) #icatage is in thousands of indv

## Read in weight at age over years
W = as.matrix(read.csv("W_at_age_matrix.csv", row.names=1)) #in kg

## Read in proportion mature and reshape it
PM = as.matrix(read.csv("prop_mature.csv", row.names=1))
PM_mat = matrix(data=PM, nrow=56, ncol=10, byrow=TRUE) # change the rows and columns to be the same as the W matrix

## Read in natural mortality 
M = read.csv("nat_mort.csv", row.names=1)$M

## Read in survey CPUE for the CPUE tuning
CPUE = as.matrix(read.csv("IBTS_Q1.csv",  row.names=1))


## Read in some functions we need
source("vpa_function.R")

## We also need to define the length of the survey CPUE time-series
first_year_CPUE = 1983
last_year_CPUE = 2019

#We'll assume that the terminal year 2018 had similar 
#catchability and recruitment as in the most recent years  
years = as.character(first_year_CPUE:(last_year_CPUE - 1)) # without the terminal year
years_to_T = as.character(first_year_CPUE:(last_year_CPUE))

#### ------- THIS PART I AM CONFUSED WHAT SHOULD THEY BE?????????????
## We need to make some rough initial guess 
Fterm_guess = c(0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

model = vpa(C=catage, M=M, Fterm=Fterm_guess, Fages=5)
str(model) #see the complicated structure of the results
str(model$N)

N_est = model$N
Z_est = model$Z


#Forecast N in 2019: recruits, survivors
#Assume recruitment will be similar to recent years
Nage_1 = mean( N_est[years, "age1"])

#Use the decay function to forecast older ages: Nt+1 = Nt*exp(-Z)
survivors = N_est["2018", ]*exp(- Z_est["2018", ])

#new recruits, survivors, plus group
N_est2019 = c(Nage_1, survivors[1:8], sum(survivors[9:10]))

#Bind the forecast to the bottom of N_est
N_est_new = rbind( N_est, "2019"=N_est2019)

# Use CPUE from IBTS_Q1 and N from the cohort model to calculate log_q = log(CPUE)-log(N_est_new) by age and year. Note that the survey started later than catch data. 
# Only use years and ages that we have both catch and survey data.
dimnames(N_est_new)
dimnames(CPUE)

log_q = log(CPUE)-log(N_est_new[21:57,1:6]) 

#Step 2: Then average log(q) across recent years to get average log(q) by age (Avg_log_q)
Avg_log_q = colMeans(log_q[years, ])
Avg_log_q_mat = matrix(rep(Avg_log_q, length(years_to_T)), nrow= length(years_to_T), byrow=TRUE)

objective = function(Fterm1to10) {
  Fterm = Fterm1to10
  Fterm["age10"] = mean(Fterm[6:9]) #F at age 7 = mean of F at ages 3-6
  model = vpa(C=catage, M=M, Fterm= Fterm, Fages=3)
  N_est = model[['N']]
  Z_est = model[['Z']] 
  
  
  #Forecast N in 2019: recruits, survivors
  #Assume recruitment will be similar to recent years
  Nage_1 = mean( N_est[years, "age1"])
  
  #Use the decay function to forecast older ages: Nt+1 = Nt*exp(-Z)
  survivors = N_est["2018", ]*exp(- Z_est["2018", ])
  
  #new recruits, survivors, plus group
  N_est2019 = c(Nage_1, survivors[1:8], sum(survivors[9:10]))
  
  #Bind the forecast to the bottom of N_est
  N_est_new = rbind( N_est, "2019"=N_est2019)
  
  #Calculate log residuals from deterministic CPUE=q*N, i.e. log(CPUE)=log(q)+log(N)
  resids = log(CPUE[years_to_T,]) - Avg_log_q_mat- log(N_est_new[years_to_T,1:6])
  
  SSR = sum(resids^2)# sum of squared residuals
  
  return(SSR)
}
#Set starting values for the optimizer
Fterm_guess = colMeans(model$F)

#Optimize / tune the parameter values according to the objective function
opt = optim(fn = objective, par = Fterm_guess[1:10], 
            lower = 0.1, upper = 2, method =  "L-BFGS-B") #bound F values

opt

#Store the optimized Fterm and the VPA that comes from it
Fterm= opt$par
Fterm = c(Fterm, mean(Fterm[6:9])) 
Fmodel = vpa(catage, M=M, Fterm= Fterm, Fages=5) #THIS DOESNT RUN!


###################################################
########### PLOT RESULTS ######################
###################################################

#Make a data frame (results) with one column for year and other columns for each of the following: SSB, Recruitment, Fbar, Catch (in biomass units)
results = data.frame(Year = as.integer(rownames(catage)))

N_est = model$N

## 1 Calculate SSB for all years
results$SSB = rowSums(N_est[,1:10]*W*PM_mat)

## 2 Recruitment for all years
results$Recruitment = N_est[,1] #leave blank so it is all rows ,1 (columns 1)

F_est = model$F

## 3 Fbar for ages 2 to 4 by year
results["F"] = rowMeans(F_est[,2:4]) 

## 4 Catch for all years #rowsum, 
# Hint: use input data catage and W
results$Catch = rowSums(catage[,1:10]*W)

str(results)

#Reshape data to make a multiplanel plot
results2 = melt(results, id.vars="Year")

#Multipanel plot
ggplot(data = results2 , aes(Year, value))+
  geom_line(data=subset(results2, variable %in% c("SSB","F")))+
  geom_col(data=subset(results2, variable %in% c("Recruitment","Catch")))+
  facet_wrap(~variable, scale="free", labeller = label_parsed)+
  ylab(NULL)

# Nicer ggplots
library("patchwork")
library("tidyverse")


Catch_plot <- results %>%
              ggplot(aes(x=Year, y=Catch)) +
                geom_bar(stat="identity", fill="steelblue",alpha=0.75)+
              theme_minimal()+ 
            ggtitle("Catches") +
            xlab(" ") + ylab("Catches ( )")+
            theme(plot.title = element_text(hjust = 0.5))

Catch_plot

Rec_plot <- results %>%
  ggplot(aes(x=Year, y=Recruitment)) +
  geom_bar(stat="identity", fill="grey")+
  theme_minimal()+ 
  ggtitle("Recruitment") +
  xlab(" ") + ylab("Recruitment ( )")+
  theme(plot.title = element_text(hjust = 0.5))

Rec_plot

F_plot <- results %>%
  ggplot(aes(x=Year, y=F)) +
  geom_line()+
  theme_minimal()+ 
  ggtitle("Fishing pressure") +
  xlab(" ") + ylab("F (ages 2-4)")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 1.3))+
  geom_hline(yintercept=.54, linetype="dashed", color = "black")+
  geom_hline(yintercept=.39, linetype="dotted", color = "black")+
  geom_hline(yintercept=.31, color = "orange")

F_plot

SSB_plot <- results %>%
  ggplot(aes(x=Year, y=SSB)) +
  geom_line()+
  theme_minimal()+ 
  ggtitle("Spawning Stock Biomass") +
  xlab(" ") + ylab("SSB ( )")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 3.5e05))+
  geom_hline(yintercept=1.07e05, linetype="dashed", color = "black")+
  geom_hline(yintercept=1.5e05, color = "orange")+
  geom_hline(yintercept=1.5e05, linetype="dotted", color = "black")



Plot_full <- (Catch_plot + Rec_plot) / (F_plot + SSB_plot)
  
Plot_full


#Plot Recruitment vs SSB
ggplot(results, aes(SSB, Recruitment))+
  geom_text(aes(label=Year))

# callmeans
## Calculate how exploitation varies by age
## hint: take F by age average across years, 
## then rescale by dividing the vector by its maximum

Exp <- colMeans(F_est[, ])/max(Exp)
Exp

#    age1      age2      age3      age4      age5      age6      age7 
# 0.4257367 1.0000000 0.9642775 0.8130401 0.6817649 0.5765303 0.8199515 

#Saving results
#option 1
save(model, results, file="NScod_VPA_results.Rdata")
#then use load("NScod_VPA_results.Rdata") to get them back

#option 2
write.csv(model$N, file="NScod_VPA_N.csv")
write.csv(model$F, file="NScod_VPA_F.csv")
write.csv(model$Z, file="NScod_VPA_Z.csv")
write.csv(results, file="NScod_VPA_results.csv")
#then use read.csv to get them back


