"# Pharmacogenetic-Prediction-Models-Code" 

The analysis done in this paper was performed using the R (v3.4.3) programming language executed on a (Linux) CentOs operating system.  As specified below, some scripts were executed as single (stand alone) scripts within R and some were submitted as jobs using Slurm on the University of Alabama at Birmingham's Cheaha computer.  Although the code may be modified to run as scripts within R, the simulations and models take considerable time to run, necessitating the use of a job system.

The following are instructions for the user:

Install necessary R packages: SimCorrMix, BhGLM, earth, ggplot2, gridExtra, grid, plyr, dplyr, xtable, Rmisc, psych

Create folder "Main" and Copy "bootstrap.bh2.R" and "glmNet2.R" to this folder.


Normal(0, 1) error terms:

1) Create folder "Main/Normal_Errors" to hold datasets.  Copy SIM_VARS_NORMAL, "sim_vars_normal.R", SIM_Y, and "sim_y.R" to this folder.

2) Simulate data sets: Within "Main/Normal_Errors", submit array job SIM_VARS_NORMAL to execute "sim_vars_normal.R". # 36 jobs run in array

3) Create 3 folders "Main/Normal_Errors/Beta1", "Main/Normal_Errors/Beta2", "Main/Normal_Errors/Beta3" to hold results from small, medium, and large treatment effects.

4) Create data sets with outcomes for three different treatment effect sizes: Within "Main/Normal_Errors", submit array job SIM_Y to execute "sim_y.R". # 36 jobs run in array

5) Run the models: Within each of the folders "Main/Normal_Errors/Beta1", "Main/Normal_Errors/Beta2", "Main/Normal_Errors/Beta3", copy and submit the following jobs and execute the following scripts:

	a) LASSO to execute Lasso.R # 36 jobs run in array
	
	b) ELASTICNET25 to execute EN25.R # 36 jobs run in array
	
	c) ELASTICNET50 to execute EN50.R # 36 jobs run in array
	
	d) ELASTICNET75 to execute EN75.R # 36 jobs run in array
	
	e) BLASSO to execute BLasso.R # 36 jobs run in array
	
	f) MARS_MODELS to execute MARS.R # 180 jobs run in array

6) Summarize results for each treatment effect size: Copy "summaries_beta1_normal.R" to "Main/Normal_Errors/Beta1" folder; 
   Copy "summaries_beta2_normal.R" to "Main/Normal_Errors/Beta2" folder; 
   Copy "summaries_beta3_normal.R" to "Main/Normal_Errors/Beta3" folder.
   Execute each of scripts within those folders.

7) Summarize results: Copy "summaries_normal.R" to "Main/Normal_Errors/" folder.  Execute the script.

Logistic(0, 1) error terms:

1) Create folder "Main/Logistic_Errors" to hold datasets.  Copy SIM_VARS_LOGISTIC, "sim_vars_logistic.R", SIM_Y, and "sim_y.R" to this folder.

2) Simulate data sets: Within "Main/Logistic_Errors", submit array job SIM_VARS_LOGISTIC to execute "sim_vars_logistic.R". # 36 jobs run in array

3) Create 3 folders "Main/Logistic_Errors/Beta1", "Main/Logistic_Errors/Beta2", "Main/Logistic_Errors/Beta3" to hold results from small, medium, and large treatment effects.

4) Create data sets with outcomes for three different treatment effect sizes: Within "Main/Logistic_Errors", submit array job SIM_Y to execute "sim_y.R". # 36 jobs run in array

5) Run the models: Within each of the folders "Main/Logistic_Errors/Beta1", "Main/Logistic_Errors/Beta2", "Main/Logistic_Errors/Beta3", copy and submit the following jobs and execute the following scripts:

	a) LASSO to execute Lasso.R # 36 jobs run in array
	
	b) ELASTICNET25 to execute EN25.R # 36 jobs run in array
	
	c) ELASTICNET50 to execute EN50.R # 36 jobs run in array
	
	d) ELASTICNET75 to execute EN75.R # 36 jobs run in array
	
	e) BLASSO to execute BLasso.R # 36 jobs run in array
	
	f) MARS_MODELS to execute MARS.R # 180 jobs run in array

6) Summarize results for each treatment effect size: Copy "summaries_beta1_logistic.R" to "Main/Logistic_Errors/Beta1" folder; 
   Copy "summaries_beta2_logistic.R" to "Main/Logistic_Errors/Beta2" folder; 
   Copy "summaries_beta3_logistic.R" to "Main/Logistic_Errors/Beta3" folder.
   Execute each of scripts within those folders.

7) Summarize results: Copy "summaries_logistic.R" to "Main/Logistic_Errors/" folder.  Execute the script.
