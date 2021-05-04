
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

