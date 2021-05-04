
# ---- Load needed packages 

library("ggplot2")
library("tidyverse")
library("patchwork")

#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
#  ---- ---- ---- ---- ----  Load data ---- ---- ---- ---- ---- 
#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

## Read catch at age, change units
catage = as.matrix(read.csv("catch_at_age_matrix.csv", row.names=1)) # catage is in thousands of indv

## Read in weight at age over years
W = as.matrix(read.csv("W_at_age_matrix.csv", row.names=1)) #in kg... so SSB and landings will be in tonnes after N*W with W in kg

## Read in proportion mature and reshape it
PM = as.matrix(read.csv("prop_mature.csv", row.names=1))
PM_mat = matrix(data=PM, nrow=56, ncol=10, byrow=TRUE) # change the rows and columns to be the same as the W matrix

## Read in natural mortality 
M = read.csv("nat_mort.csv", row.names=1)$M

## Read in survey CPUE for the CPUE tuning
CPUE = as.matrix(read.csv("IBTS_Q1.csv",  row.names=1))

#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
#  ---- ---- ---- ---- ----  CPUE tuning of terminal F ---- ---- ---- ---- ---- 
#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

## Read in some functions we need
source("vpa_function.R")

## Some rough initial guess 
Fterm_guess = c(0.3, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

## We also need to define the length of the survey CPUE time-series
first_year_CPUE = as.numeric(rownames(CPUE)[1])
last_year_CPUE = as.numeric(rownames(CPUE[nrow(CPUE), , drop = FALSE]))

# Code based off of CPU_tuning_script.R

#We'll assume that the terminal year 2018 had similar 
# catchability and recruitment as in the most recent years  
years = as.character(first_year_CPUE:(last_year_CPUE - 1)) # without the terminal year
years_to_T = as.character(first_year_CPUE:(last_year_CPUE)) # last year

model = vpa(C=catage, M=M, Fterm=Fterm_guess, Fages=3) 
str(model) # see the complicated structure of the results
str(model$N)

N_est = model$N # N estimates
Z_est = model$Z # Z -> total mortality

#Forecast N in 2019: recruits, survivors

#Assume recruitment will be similar to the ten recent years
Nage_1 = mean( N_est[years[27:36], "age1"]) 

#Use the decay function to forecast older ages: Nt+1 = Nt*exp(-Z)

survivors = tail(N_est,1)*exp(- tail(Z_est,1))  # tail is pulling year 2018

#new recruits, survivors, plus group
N_est2019 = c(Nage_1, survivors[1:(length(survivors)-2)], sum(survivors[,(ncol(survivors)-1):ncol(survivors)]))

#Bind the forecast to the bottom of N_est
N_est_new = rbind( N_est, "2019"=N_est2019)

# Use CPUE from IBTS_Q1 and N from the cohort model to calculate log_q = log(CPUE)-log(N_est_new) by age and year.
# Note that the survey started later than catch data. so...Only use years and ages that we have both catch and survey data.
dimnames(N_est_new)
dimnames(CPUE)
# only ages 1-6 overlap between the two (columns 1-6)
# only years 1983 - 2019 overlap between the two (rows 21-57 in N_est)
log_q = log(CPUE)-log(N_est_new[21:57,1:6]) 

#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
#  ---- ---- ---- ---- ----  Step 2 ---- ---- ---- ---- ---- 
#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

# Then average log(q) across recent years to get average log(q) by age (Avg_log_q)
Avg_log_q = colMeans(log_q[years, ])
Avg_log_q_mat = matrix(rep(Avg_log_q, length(years_to_T)), nrow= length(years_to_T), byrow=TRUE)

#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
#  ---- ---- ----   now the objective function ---- ---- 
#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 

objective = function(Fterm1to10) {
  Fterm = Fterm1to10
  Fterm["age10"] = mean(Fterm[7:9]) #F at age 10 = mean of F at ages 7-9
  model = vpa(C=catage, M=M, Fterm= Fterm, Fages=3)
  N_est = model[['N']]
  Z_est = model[['Z']] 
  
  
  #Forecast N in 2019: recruits, survivors
  #Assume recruitment will be similar to recent ten years
  Nage_1 = mean( N_est[years[27:36], "age1"])
  
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
opt = optim(fn = objective, par = Fterm_guess[1:9], 
            lower = 0.1, upper = 2, method =  "L-BFGS-B") #bound F values

opt

#Store the optimized Fterm and the VPA that comes from it
Fterm= opt$par
Fterm = c(Fterm, mean(Fterm[7:9])) 

## see Fterm that you can then use in the vpa-function below
Fterm

model = vpa(catage, M=M, Fterm= Fterm, Fages=3)


#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 
#  ---- ---- ---- ---- ----  PLOT RESULTS ---- ---- ---- ---- ---- 
#  ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- 


#Make a data frame (results) with one column for year and other columns for each of the following: SSB, Recruitment, Fbar, Catch (in biomass units)
results = data.frame(Year = as.integer(rownames(catage)))

N_est = model$N #note that this is N on January 1 in this example
Z_est = model$Z #this is total mortality integrated over the entire year 
F_est = model$F #this is fishing mortality integrated over the entire year 


## 1) Calculate SSB for all years
results$SSB = rowSums(model$N*PM_mat*W)

## 2) Recruitment for all years
results$Recruitment = model$N[,1] #leave blank so it is all rows ,1 (columns 1)

## 3) Fbar for ages 2 to 4 by year
results$F_2to4 = rowMeans(model$F[,2:4])

## 4) Catch for all years #rowsum, 
# Hint: use input data catage and W
results$Catch = rowSums(catage * W)

str(results)

Catch_plot <- results %>%
  ggplot(aes(x=Year, y=Catch)) +
  geom_bar(stat="identity", fill="steelblue",alpha=0.75)+
  theme_classic()+ 
  ggtitle("Catches") +
  xlab(" ") + ylab("Catches (weight in kg)")+
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_x_continuous(breaks=seq(1963, 2018, 10))+
  scale_y_continuous(breaks=c(0,1.0e05, 2.0e05, 3.0e05,4e05,5e05,6e05),limits = c(0,6e05))


Catch_plot

Rec_plot <- results %>%
  ggplot(aes(x=Year, y=Recruitment)) +
  geom_bar(stat="identity", fill="grey")+
  theme_classic()+ 
  ggtitle("Recruitment") +
  xlab(" ") + ylab("Recruitment (No. Indv.)")+
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_x_continuous(breaks=seq(1963, 2018, 10))+
  scale_y_continuous(breaks=c(0,0.5e06,1.0e06, 1.5e06, 2.0e06, 2.5e06, 3.0e06, 3.5e06),limits = c(0,3.5e06))


Rec_plot

F_plot <- results %>%
  ggplot(aes(x=Year, y=F_2to4)) +
  geom_line()+
  theme_classic()+
  ggtitle("Fishing pressure") +
  xlab(" ") + ylab("F (ages 2-4)")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks=seq(1963, 2018, 10))+ 
  scale_y_continuous(breaks=c(0,.2, .4, .6, .8, 1.0, 1.2),limits = c(0,1.4))+ 
  geom_hline(yintercept=0.554, color = "orange")+ 
  geom_label(
    label="Fmsy = 0.554", 
    x=1983,
    y=.55,
    label.size = 0.1,
    color = "orange"
  )


F_plot

SSB_plot <- results %>%
  ggplot(aes(x=Year, y=SSB)) +
  geom_line()+
  theme_classic()+
  ggtitle("Spawning Stock Biomass") +
  xlab(" ") + ylab("SSB (weight in kg)")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  scale_x_continuous(breaks=seq(1963, 2018, 10))+
  scale_y_continuous(breaks=c(0,0.5e05,1.0e05, 1.5e05, 2.0e05, 2.5e05, 3.0e05, 3.5e05),limits = c(0,3.5e05))+
  geom_hline(yintercept=107814.8, color = "orange")+
  geom_label(
    label="Bpa = 107814.8", 
    x=2013,
    y=107814.8,
    label.size = 0.1,
    color = "orange"
  )+
  geom_hline(yintercept=65819.41, color = "lightblue")+
  geom_label(
    label="Blim = 65819.41", 
    x=1973,
    y=65819.41,
    label.size = 0.1,
    color = "lightblue"
  )

SSB_plot

Plot_full <- (Catch_plot + Rec_plot) / (F_plot + SSB_plot)

Plot_full+
  plot_annotation(
    title = 'Stock Development Overtime',
    theme = theme(plot.title = element_text(hjust = 0.5)),
    caption = "Cod in Subarea 4, Dividion 7.d, and Subdiveision 20")


#############################################################################################################################
########### Make a stock-recruitment plot, so we can say something about Blim (see ICES guidelines). ######################## 
########### Note shift the time-series so that I match up SSB and Recruitment at age-1 (i.e. SSB produce recruits that are ##               ##
########### 0-year old, but we don't see them until they are 1-year old) ####################################################                                                ##
#############################################################################################################################

Years = 1963:2018
Recruitment = model$N[,1]
SSB = rowSums(model$N*PM_mat*W)

S = data.frame(SSB[1:55])
names(S)[names(S) == "SSB.1.55."] <- "SSB"
rownames(S) <- NULL

R = data.frame(Recruitment[2:56])
names(R)[names(R) == "Recruitment.2.56."] <- "Recruitment"
rownames(R) <- NULL

Y = data.frame(Years[1:55])
names(Y)[names(Y) == "Years.1.55."] <- "Year"

# make a dataframe of the results: 
rickercurve <- cbind(S,R,Y)

# do a linear model of SSB (predictor) vs log(R/SSB) (response)
mod = lm(log(Recruitment/SSB)~SSB, rickercurve)
mod 
# use the linear model to calculate a and b for the Ricker curve
# a = exp(intercept), b = -slope
a = exp(1.298e+00)
b = -1.059e-06  

# get predictions of R either using the coefficients,
# or use the linear model to make predictions of log(R/SSB) and transform them to R
rickercurve$pred = predict(mod)
rickercurve = transform(rickercurve, predR = exp(pred)*SSB)

# use the linear model to calculate a and b for the Ricker curve
# a = exp(intercept), b = -slope

#Plot Recruitment vs SSB
ggplot(rickercurve, aes(SSB, Recruitment))+
  geom_text(aes(label=Year))+ 
  scale_x_continuous(breaks=c(0,0.5e05,1.0e05, 1.5e05, 2.0e05, 2.5e05, 3.0e05, 3.5e05),limits = c(0,3.5e05))+
  scale_y_continuous(breaks=c(0,0.5e06,1.0e06, 1.5e06, 2.0e06, 2.5e06, 3.0e06),limits = c(0,3e06))+ 
  geom_line(aes(y=predR))+
  annotate(xmin = 50000, xmax = 80000, 
           ymin = 840000, ymax = 1000000, 
           geom = "rect", alpha = 0.2,
           fill = "red")+
  theme_classic()


#By changing Fmult we can change the fishing pressure, but the exploitation pattern (relative difference between age groups) remains the same.
# If Fmult=1 we continue with the same F (and exploitation pattern). 

F_mult=1  

dat_2019 = data.frame(
  W_5year_ave = colMeans(W[52:56,]),
  M = M,
  prop_mature = PM,
  N = N_est2019,
  F = Fterm*F_mult,
  Z = (Fterm*F_mult)+M
)


mean(dat_2019$F[2:4]) #Often when we present F, we present an average accross selected age groups (in this case I chose age 2-4)
# 0.7919829 

dat_2019 = transform(dat_2019,  C=F/(Z)*N*(1-exp(-(Z)))) #Here I calculate Catch numbers at age (using the catch equation) and place it in dat_2019

SSB_2019 = sum(dat_2019$N * dat_2019$W_5year_ave * dat_2019$Proportion_Mature) #Here I calculate SSB in 2008 (January 1)
# 71350.21

landings_2019 = sum(dat_2019$C * dat_2019$W_5year_ave) #Here I calculate catch (or landings or yield; has many names) in tons
# 63040.84


######################################################################
########### Calculate Fmsy ###########################################
######################################################################


# do a linear model of SSB (predictor) vs log(R/SSB) (response)
mod = lm(log(Recruitment/SSB)~SSB, rickercurve)

mod 

#this is the parameters for the Ricker-model
a = exp(mod$coefficients[1]) 
b = -(mod$coefficients[2]) 

#The fucntion is called calc_assymptotic_SSB_landings(). 
# Here we use the Ricker stock-recruitment model instead of a Beverton & Holt (see Fmsy_function.R)

calc_assymptotic_SSB_landings = function(Fmult, N=dat_2019$N, W=dat_2019$W, M=dat_2019$M, PM=dat_2019$prop_mature, E=dat_2019$F, nyears=1) {
  F = E*Fmult
  Z = F + M
  for(y in 1:nyears){
    SSB = sum(N * W * PM)
    R = a*SSB*exp(-b*SSB)  #Ricker assumption
    C = F/(Z)*N*(1-exp(-Z))
    survivors = N * exp(-Z)
    Nnew = c(R, survivors[1:8], sum(survivors[9:10]))
    SSB = sum(Nnew * W * PM)
    landings = sum(C * W)
    N = Nnew  #update the N at age for next year of the for-loop
  }
  return(list(SSB, landings)) #after nyears what are the values
}


F_mult_range = seq(0,2, by=.1) #We need to define a range of Fmult that we loop through to figure out which Fmult gives us Fmsy

long_term_results = data.frame(Fmult= F_mult_range, SSB = NA, yield = NA) #We produce an empty data frame to be filled in using the following code


#Here we loop through all the different values of Fmult specified by F_mult_range
for(i in 1:length(F_mult_range)){
  long_term_results[i,2:3] = 
    calc_assymptotic_SSB_landings(Fmult= F_mult_range[i], nyears=1000, N=dat_2019$N, W=dat_2019$W, M=dat_2019$M, PM=dat_2019$Proportion_Mature, E=dat_2019$F)
}

#Here we identify the Fmult that gives us maximum yield and then we calculate the corresponding F values (per age group) and lastly we average accross age 2-4 and this gives us Fmsy
XXX = which.max(long_term_results$yield)
F_MSY_all_ages = dat_2019$F * long_term_results[XXX,1]
Fmsy = mean(F_MSY_all_ages[2:4])
Fmsy
# 0.554388

long_term_plot <- long_term_results %>%
  ggplot(aes(x=Fmult, y=yield)) +
  geom_line()

long_term_plot <- long_term_results %>%
  ggplot(aes(x=Fmult, y=SSB)) +
  geom_line()

ggplot() + 
  geom_line(data = long_term_results, aes(x = Fmult, y = yield, colour = "Yield")) +
  geom_line(data = long_term_results, aes(x = Fmult, y = SSB ,colour = "SSB")) +
  labs(y = "1000s of tonnes",
       x = "Fmult",
       colour = "Parameter")+ 
  theme(legend.position = c(0.8, 0.8))+
  theme_classic()

# ============================================================= #
#### now lets project diff. outcomes for yrs with diff fmult.
# ============================================================= #

# Now lets make F effort at the FMSY 

### Fmult is .7 --> @ Fsmy
F_mult= .7 

dat_2019 = data.frame(
  W_5year_ave = colMeans(W[52:56,]),
  M = M,
  prop_mature = PM,
  N = N_est2019,
  F = Fterm*F_mult,
  Z = (Fterm*F_mult)+M
)


mean(dat_2019$F[2:4]) #Often when we present F, we present an average accross selected age groups (in this case I chose age 2-4)

dat_2019 = transform(dat_2019,  C=F/(Z)*N*(1-exp(-(Z)))) #Here I calculate Catch numbers at age (using the catch equation) and place it in dat_2019

SSB_2019 = sum(dat_2019$N * dat_2019$W_5year_ave * dat_2019$Proportion_Mature) #Here I calculate SSB in 2008 (January 1)

landings_2019 = sum(dat_2019$C * dat_2019$W_5year_ave)

#  --------- now 2020 ---------- #

dat_2020 = data.frame(
  W_5year_ave = dat_2019$W_5year_ave,
  M = dat_2019$M,
  prop_mature = dat_2019$Proportion_Mature,
  N = c(dat_2019$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2020$F[2:4])

# note that we assume that recruitment in 2020 (age-1) is the same is in 2019 
dat_2020[1,4]

#Here we calculate the survivors from 2019 and insert them into dat_2020
survivors = dat_2019$N * exp(-dat_2019$Z)
dat_2020$N[2:9]=survivors[1:8]
dat_2020$N[10]=sum(survivors[9:10])

dat_2020 = transform(dat_2020,  C=F/(Z)*N*(1-exp(-(Z))))
SSB_2020 = sum(dat_2020$N * dat_2020$W_5year_ave * dat_2020$prop_mature)

landings_2020 = sum(dat_2020$C * dat_2020$W_5year_ave)

#  --------- now 2021 ---------- #

dat_2021 = data.frame(
  W_5year_ave = dat_2020$W_5year_ave,
  M = dat_2020$M,
  prop_mature = dat_2020$prop_mature,
  N = c(dat_2020$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2021$F[2:4])

# note that we assume that recruitment in is the same is in 2019 
dat_2021[1,4]

#Here we calculate the survivors from 2019 and insert them into dat_2020
survivors = dat_2020$N * exp(-dat_2020$Z)
dat_2021$N[2:9]=survivors[1:8]
dat_2021$N[10]=sum(survivors[9:10])

dat_2021 = transform(dat_2021,  C=F/(Z)*N*(1-exp(-(Z))))
SSB_2021 = sum(dat_2021$N * dat_2021$W_5year_ave * dat_2021$prop_mature)

landings_2021 = sum(dat_2021$C * dat_2021$W_5year_ave)

#  --------- now 2022 ---------- #

dat_2022 = data.frame(
  W_5year_ave = dat_2021$W_5year_ave,
  M = dat_2021$M,
  prop_mature = dat_2021$prop_mature,
  N = c(dat_2021$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2022$F[2:4])

# note that we assume that recruitment in 2020 (age-1) is the same is in 2019 
dat_2022[1,4]

#Here we calculate the survivors from 2019 and insert them into dat_2021
survivors = dat_2021$N * exp(-dat_2021$Z)
dat_2022$N[2:9]=survivors[1:8]
dat_2022$N[10]=sum(survivors[9:10])

dat_2022 = transform(dat_2022,  C=F/(Z)*N*(1-exp(-(Z))))
SSB_2022 = sum(dat_2022$N * dat_2022$W_5year_ave * dat_2022$prop_mature)

landings_2022 = sum(dat_2022$C * dat_2022$W_5year_ave)

#  --------- now 2023 ---------- #

dat_2023 = data.frame(
  W_5year_ave = dat_2022$W_5year_ave,
  M = dat_2022$M,
  prop_mature = dat_2022$prop_mature,
  N = c(dat_2022$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2023$F[2:4])

dat_2023[1,4]

survivors = dat_2022$N * exp(-dat_2022$Z)
dat_2023$N[2:9]=survivors[1:8]
dat_2023$N[10]=sum(survivors[9:10])

dat_2023 = transform(dat_2023,  C=F/(Z)*N*(1-exp(-(Z))))
SSB_2023 = sum(dat_2023$N * dat_2023$W_5year_ave * dat_2023$prop_mature)

landings_2023 = sum(dat_2023$C * dat_2023$W_5year_ave)

#  --------- now 2024 ---------- #

dat_2024 = data.frame(
  W_5year_ave = dat_2023$W_5year_ave,
  M = dat_2023$M,
  prop_mature = dat_2023$prop_mature,
  N = c(dat_2023$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2024$F[2:4])

dat_2024[1,4]

survivors = dat_2023$N * exp(-dat_2023$Z)
dat_2024$N[2:9]=survivors[1:8]
dat_2024$N[10]=sum(survivors[9:10])

dat_2024 = transform(dat_2024,  C=F/(Z)*N*(1-exp(-(Z))))
SSB_2024 = sum(dat_2024$N * dat_2024$W_5year_ave * dat_2024$prop_mature)

landings_2024 = sum(dat_2024$C * dat_2024$W_5year_ave)


#  ------ ------ ------ ------ ------ ------ ------ ------ ------  ------ ------ ------  ------ ------
#  ------ ------ Different lower Fmult values - recover is still very slow @ Fmsy ------  ------ ------
# ------ ------ ------ ------ ------ ------ ------ ------ ------  ------ ------ ------  ------ ------

F_mult= .45

dat_2019_lower = data.frame(
  W_5year_ave = colMeans(W[52:56,]),
  M = M,
  prop_mature = PM,
  N = N_est2019,
  F = Fterm*F_mult,
  Z = (Fterm*F_mult)+M
)


mean(dat_2019_lower$F[2:4]) #Often when we present F, we present an average accross selected age groups (in this case I chose age 2-4)

dat_2019_lower = transform(dat_2019_lower,  C=F/(Z)*N*(1-exp(-(Z)))) #Here I calculate Catch numbers at age (using the catch equation) and place it in dat_2019

SSB_2019_lower = sum(dat_2019_lower$N * dat_2019_lower$W_5year_ave * dat_2019_lower$Proportion_Mature) #Here I calculate SSB in 2008 (January 1)

landings_2019_lower = sum(dat_2019_lower$C * dat_2019_lower$W_5year_ave)


#  --------- now 2020 ---------- #

dat_2020_lower = data.frame(
  W_5year_ave = dat_2019_lower$W_5year_ave,
  M = dat_2019_lower$M,
  prop_mature = dat_2019_lower$Proportion_Mature,
  N = c(dat_2019_lower$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2020_lower$F[2:4])

# note that we assume that recruitment in 2020 (age-1) is the same is in 2019 
dat_2020_lower[1,4]

#Here we calculate the survivors from 2019 and insert them into dat_2020
survivors = dat_2019_lower$N * exp(-dat_2019_lower$Z)
dat_2020_lower$N[2:9]=survivors[1:8]
dat_2020_lower$N[10]=sum(survivors[9:10])

dat_2020_lower = transform(dat_2020_lower,  C=F/(Z)*N*(1-exp(-(Z))))

SSB_2020_lower = sum(dat_2020_lower$N * dat_2020_lower$W_5year_ave * dat_2020_lower$prop_mature)

landings_2020_lower = sum(dat_2020_lower$C * dat_2020_lower$W_5year_ave)


#  --------- now 2021 ---------- #

dat_2021_lower = data.frame(
  W_5year_ave = dat_2020_lower$W_5year_ave,
  M = dat_2020_lower$M,
  prop_mature = dat_2020_lower$prop_mature,
  N = c(dat_2020_lower$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2021$F[2:4])

# note that we assume that recruitment in is the same is in 2019 
dat_2021_lower[1,4]

survivors = dat_2020_lower$N * exp(-dat_2020_lower$Z)
dat_2021_lower$N[2:9]=survivors[1:8]
dat_2021_lower$N[10]=sum(survivors[9:10])

dat_2021_lower = transform(dat_2021_lower,  C=F/(Z)*N*(1-exp(-(Z))))
SSB_2021_lower = sum(dat_2021_lower$N * dat_2021_lower$W_5year_ave * dat_2021_lower$prop_mature)

landings_2021_lower = sum(dat_2021_lower$C * dat_2021_lower$W_5year_ave)

#  --------- now 2022 ---------- #

dat_2022_lower = data.frame(
  W_5year_ave = dat_2021_lower$W_5year_ave,
  M = dat_2021_lower$M,
  prop_mature = dat_2021_lower$prop_mature,
  N = c(dat_2021_lower$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2022_lower$F[2:4])

# note that we assume that recruitment in 2020 (age-1) is the same is in 2019 
dat_2022_lower[1,4]

survivors = dat_2021_lower$N * exp(-dat_2021_lower$Z)
dat_2022_lower$N[2:9]=survivors[1:8]
dat_2022_lower$N[10]=sum(survivors[9:10])

dat_2022_lower = transform(dat_2022_lower,  C=F/(Z)*N*(1-exp(-(Z))))
SSB_2022_lower = sum(dat_2022_lower$N * dat_2022_lower$W_5year_ave * dat_2022_lower$prop_mature)

landings_2022_lower = sum(dat_2022_lower$C * dat_2022_lower$W_5year_ave)

#  --------- now 2023 ---------- #

dat_2023_lower = data.frame(
  W_5year_ave = dat_2022_lower$W_5year_ave,
  M = dat_2022_lower$M,
  prop_mature = dat_2022_lower$prop_mature,
  N = c(dat_2022_lower$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2023_lower$F[2:4])

dat_2023_lower[1,4]

survivors = dat_2022_lower$N * exp(-dat_2022_lower$Z)
dat_2023_lower$N[2:9]=survivors[1:8]
dat_2023_lower$N[10]=sum(survivors[9:10])

dat_2023_lower = transform(dat_2023_lower,  C=F/(Z)*N*(1-exp(-(Z))))
SSB_2023_lower = sum(dat_2023_lower$N * dat_2023_lower$W_5year_ave * dat_2023_lower$prop_mature)

landings_2023_lower = sum(dat_2023_lower$C * dat_2023_lower$W_5year_ave)

#  --------- now 2024 ---------- #

dat_2024_lower = data.frame(
  W_5year_ave = dat_2023_lower$W_5year_ave,
  M = dat_2023_lower$M,
  prop_mature = dat_2023_lower$prop_mature,
  N = c(dat_2023_lower$N[1], NA, NA, NA, NA, NA, NA, NA, NA, NA),
  F = (Fterm*F_mult),
  Z = (Fterm*F_mult)+M
)

mean(dat_2024_lower$F[2:4])

dat_2024_lower[1,4]

survivors = dat_2023_lower$N * exp(-dat_2023_lower$Z)
dat_2024_lower$N[2:9]=survivors[1:8]
dat_2024_lower$N[10]=sum(survivors[9:10])

dat_2024_lower = transform(dat_2024_lower,  C=F/(Z)*N*(1-exp(-(Z))))
SSB_2024_lower = sum(dat_2024_lower$N * dat_2024_lower$W_5year_ave * dat_2024_lower$prop_mature)

landings_2024_lower = sum(dat_2024_lower$C * dat_2024_lower$W_5year_ave)



# ===== ICES plots with the projected data as well ======= #

Catch_plot <- results %>%
  ggplot(aes(x=Year, y=Catch)) +
  geom_bar(stat="identity", fill="steelblue",alpha=0.75)+
  theme_classic()+
  ggtitle("Catches") +
  xlab(" ") + ylab("Catches (weight in tonnes)")+
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_x_continuous(breaks=seq(1963, 2018, 10))+
  scale_y_continuous(breaks=c(0,1.0e05, 2.0e05, 3.0e05,4e05,5e05,6e05),limits = c(0,6e05))+
  geom_point(aes(x=2019, y=landings_2019), colour="green")+
  geom_point(aes(x=2020, y=landings_2020), colour="green")+
  geom_point(aes(x=2021, y=landings_2021), colour="green")+
  geom_point(aes(x=2022, y=landings_2022), colour="green")+
  geom_point(aes(x=2023, y=landings_2023), colour="green")+
  geom_point(aes(x=2024, y=landings_2024 ), colour="green")+
  geom_point(aes(x=2019, y=landings_2019_lower), colour="pink")+
  geom_point(aes(x=2020, y=landings_2020_lower), colour="pink")+
  geom_point(aes(x=2021, y=landings_2021_lower), colour="pink")+
  geom_point(aes(x=2022, y=landings_2022_lower), colour="pink")+
  geom_point(aes(x=2023, y=landings_2023_lower), colour="pink")+
  geom_point(aes(x=2024, y=landings_2024_lower ), colour="pink")

Catch_plot

Rec_plot <- results %>%
  ggplot(aes(x=Year, y=Recruitment)) +
  geom_bar(stat="identity", fill="grey")+
  theme_classic()+
  ggtitle("Recruitment") +
  xlab(" ") + ylab("Recruitment (No. Indv.)")+
  theme(plot.title = element_text(hjust = 0.5))+ 
  scale_x_continuous(breaks=seq(1963, 2018, 10))+
  scale_y_continuous(breaks=c(0,0.5e06,1.0e06, 1.5e06, 2.0e06, 2.5e06, 3.0e06, 3.5e06),limits = c(0,3.5e06))+
  geom_point(aes(x=2019, y=dat_2019[1,4]), colour="green",size =2.5)+
  geom_point(aes(x=2020, y=dat_2020[1,4]), colour="green",size =2.5)+
  geom_point(aes(x=2021, y=dat_2021[1,4]), colour="green",size =2.5)+
  geom_point(aes(x=2022, y=dat_2022[1,4]), colour="green",size =2.5)+
  geom_point(aes(x=2023, y=dat_2023[1,4]), colour="green",size =2.5)+
  geom_point(aes(x=2024, y=dat_2024[1,4]), colour="green",size =2.5)+ 
  geom_point(aes(x=2019, y=dat_2019_lower[1,4]), colour="pink")+
  geom_point(aes(x=2020, y=dat_2020_lower[1,4]), colour="pink")+
  geom_point(aes(x=2021, y=dat_2021_lower[1,4]), colour="pink")+
  geom_point(aes(x=2022, y=dat_2022_lower[1,4]), colour="pink")+
  geom_point(aes(x=2023, y=dat_2023_lower[1,4]), colour="pink")+
  geom_point(aes(x=2024, y=dat_2024_lower[1,4]), colour="pink")


Rec_plot


F_plot <- results %>%
  ggplot(aes(x=Year, y=F_2to4)) +
  geom_line()+
  theme_classic()+
  ggtitle("Fishing pressure") +
  xlab(" ") + ylab("F (ages 2-4)")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept=Fmsy, color = "orange")+ 
  geom_label(
    label="Fmsy = 0.554", 
    x=1983,
    y=.55,
    label.size = 0.1,
    color = "orange"
  )+
  scale_x_continuous(breaks=seq(1963, 2018, 10))+ 
  scale_y_continuous(breaks=c(0,.2, .4, .6, .8, 1.0, 1.2),limits = c(0,1.4)) +
  geom_point(aes(x=2019, y=mean(dat_2019$F[2:4])), colour="green" )+
  geom_point(aes(x=2020, y=mean(dat_2020$F[2:4])), colour="green")+
  geom_point(aes(x=2021, y=mean(dat_2021$F[2:4])), colour="green")+
  geom_point(aes(x=2022, y=mean(dat_2022$F[2:4])), colour="green")+
  geom_point(aes(x=2023, y=mean(dat_2023$F[2:4])), colour="green")+
  geom_point(aes(x=2024, y=mean(dat_2024$F[2:4])), colour="green")+
  geom_point(aes(x=2019, y=mean(dat_2019_lower$F[2:4])), colour="pink")+
  geom_point(aes(x=2020, y=mean(dat_2020_lower$F[2:4])), colour="pink")+
  geom_point(aes(x=2021, y=mean(dat_2021_lower$F[2:4])), colour="pink")+
  geom_point(aes(x=2022, y=mean(dat_2022_lower$F[2:4])), colour="pink")+
  geom_point(aes(x=2023, y=mean(dat_2023_lower$F[2:4])), colour="pink")+
  geom_point(aes(x=2024, y=mean(dat_2024_lower$F[2:4])), colour="pink")

F_plot


########### Inserting Blim and Bpa reference point ###################

Blim = 	65819.41 #I decided to use the 1993 SSB value based on type-1 in the ICES guidelines

Bpa = Blim * exp(0.3 * 1.645)

SSB_plot <- results %>%
  ggplot(aes(x=Year, y=SSB)) +
  geom_line()+
  theme_classic()+
  ggtitle("Spawning Stock Biomass") +
  xlab(" ") + ylab("SSB (weight in tonnes)")+
  theme(plot.title = element_text(hjust = 0.5)) + 
  #scale_y_continuous(limits = c(0, 3.5e05))+
  geom_hline(yintercept=Bpa, color = "orange")+
  geom_label(
    label="Bpa = 107814.8", 
    x=2013,
    y=Bpa,
    label.size = 0.1,
    color = "orange"
  )+
  geom_hline(yintercept=Blim, color = "lightblue")+
  geom_label(
    label="Blim = 65819.41", 
    x=1973,
    y=Blim,
    label.size = 0.1,
    color = "lightblue"
  )+
  scale_x_continuous(breaks=seq(1963, 2018, 10))+
  scale_y_continuous(breaks=c(0,0.5e05,1.0e05, 1.5e05, 2.0e05, 2.5e05, 3.0e05, 3.5e05),limits = c(0,3.5e05))+
  geom_point(aes(x=2019, y=SSB_2019), colour="green")+
  geom_point(aes(x=2020, y=SSB_2020), colour="green")+
  geom_point(aes(x=2021, y=SSB_2021 ), colour="green")+
  geom_point(aes(x=2022, y=SSB_2022 ), colour="green")+
  geom_point(aes(x=2023, y=SSB_2023), colour="green")+
  geom_point(aes(x=2024, y=SSB_2024), colour="green")+
  geom_point(aes(x=2019, y=SSB_2019_lower), colour="pink")+
  geom_point(aes(x=2020, y=SSB_2020_lower), colour="pink")+
  geom_point(aes(x=2021, y=SSB_2021_lower ), colour="pink")+
  geom_point(aes(x=2022, y=SSB_2022_lower ), colour="pink")+
  geom_point(aes(x=2023, y=SSB_2023_lower), colour="pink")+
  geom_point(aes(x=2024, y=SSB_2024_lower), colour="pink")


SSB_plot

Plot_full <- (Catch_plot + Rec_plot) / (F_plot + SSB_plot)

Plot_full+
  plot_annotation(
    title = 'Stock Development Overtime',
    theme = theme(plot.title = element_text(hjust = 0.5)),
    caption = "Cod in Subarea 4, Dividion 7.d, and Subdiveision 20")




