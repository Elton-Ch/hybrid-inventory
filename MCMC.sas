PROC MCMC DATA=Input_data OUTPOST=sfposterior DIAG=ALL DIC
		  NBI=10000 NTU=10000 NMC=100000 THIN=10 STATS=ALL
		   MCHISTORY=BRIEF PLOTS(SMOOTH)=ALL seed=12345 ;
	
	PARMS sf1-sf&num_regions 1;
	PARMS sigma2 %SYSEVALF((&SE_percent_error/100)**2);	
	PRIOR sf1-sf&num_regions ~ normal(mean=1, var=%SYSEVALF((&SP_percent_error/100)**2));	
	PRIOR sigma2 ~ igamma(shape=&shape_a,scale=&scale_b);	

	BEGINNODATA;
		std = sqrt(sigma2);
	ENDNODATA;

	mu = sf&num_regions*r&num_regions;
	
	MODEL obs ~ normal(mu, sd=std);
	
	preddist outpred=pred stats=NONE;
RUN;
