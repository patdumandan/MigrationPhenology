data {
                 
                 vector<lower=0,upper=0>[4] zeros4; // vector for random effects distribution
                 int N; //no.of observations (length of hms$year)
                 int Nsp; //no.of populations
                 real years[N];//explanatory variable (centered around mean)
                 int count[N];//no.of individuals per year
                 real obs_days[N]; // obs.effort
                 int spcode[N];// id for each population
                 
                 }
                 
                 parameters {
                 
                 vector[2] beta;            // fixed effects, intercept and slopes
                 //real<lower=-0.86,upper=0.83> bkpoint;// fixed effects, knotpoint (bounding specified to help convergence)
                 vector<lower=0>[2] sigma_sp;   // level 2 error sd per parameter
                 real<lower=0> phi;        // level 1 error sd
                 vector[2] u[Nsp];         //  level 2 error sd per species per parameter
                 cholesky_factor_corr[2] L_u_Corr; // cholesky factor for the random effects correlation matrix
                 
                 }
                 
                 transformed parameters {  
                 vector[2] alpha[Nsp];     // random effects
                 real<lower=0> y_mu[N];              // mean parameter based on regression equation
                 
                 for (i in 1:Nsp) {
                 alpha[i,1] = beta[1] + u[i,1];
                 alpha[i,2] = beta[2] + u[i,2];
                 
                 }
                 
                 
                 //=====================
                 // regression equation
                 //=====================
                 
                 for (j in 1:N){
                 
                 
                 y_mu[j]= exp(alpha[spcode[j],1]+alpha[spcode[j],2]* (years[j]+obs_days[j]));
                 
                 }
                 
                 } 
                 
                 model {
                 
                 //========
                 // priors
                 //========
                 
                 beta[1] ~ normal(0, 10); // prior: fixed effect, intercept
                 beta[2] ~ normal(0, 10);   // prior: fixed effect, slope before knot
                 
                 // prior: fixed effect, knot point
                 
                 sigma_sp[1] ~ cauchy(0,2.5);    // prior: random effect sd, intercept
                 sigma_sp[2] ~ cauchy(0,2.5);    // prior: random effect sd, slope before bkpt
                 
                 phi ~ cauchy(0,2.5);       // prior: level 1 error sd
                 
                 L_u_Corr ~ lkj_corr_cholesky(1);
                 // prior: cholesky factor for random effects correlation matrix
                 // NB. this prior is the lkj correlation distribution with shape parameter 1 
                 // which is equivalent to a uniform distribution over the possible correlation 
                 // matrices (where a shape parameter > 1 would have resulted in an upside down
                 // U-shaped distribution with the mode being located at the identity matrix)
                 
                 //=============================
                 // random effects distribution
                 //=============================
                 
                 for (i in 1:Nsp) u[i] ~ multi_normal_cholesky(zeros4, diag_pre_multiply(sigma_sp, L_u_Corr));
                 // NB. the second parameter here is the cholesky factor L 
                 // (for the correlation matrix). It only uses the sd rather 
                 // than the variances since Sigma = L*L'
                 
                 //==================
                 // model likelihood
                 //==================
                 
                 count ~ neg_binomial_2(y_mu, phi); // likelihood for the observed data
                 
                 }
                 
                 generated quantities {
                 
                 int <lower=0> y_mu_pred[N];   // predicted mean
                 corr_matrix[2] u_Corr;   // random effects correlation matrix
                 matrix[2,2] u_Sigma;     // random effects covariance matrix
                 vector[N] log_lik;
                 vector[2] alpha_tosave[Nsp];
                 // monitor random effects for a subset of patients only
                 // (for plotting predictions) and do not monitor 'alpha' 
                 // in the model above (since it consumes too much memory!)
                 
                 //==================================================
                 // predicted mean outcome using regression equation
                 //==================================================
                 
                 for (i in 1:Nsp) {  
                 alpha_tosave[i] = alpha[i];
                 }
                 
               //  y_mu_pred = neg_binomial_2_rng(y_mu, phi);
                 
                 for (n in 1:N) {
                 log_lik[n] = neg_binomial_2_lpmf(count [n] |y_mu[n],phi);
                 y_mu_pred[n] = neg_binomial_2_rng(y_mu[n], phi);
                 }
                 //=====================================================
                 // recover the correlation and covariance matrices
                 // using the cholesky factor of the correlation matrix
                 //=====================================================
                 
                 u_Corr = multiply_lower_tri_self_transpose(L_u_Corr);    
                 // correlation matrix: u_Corr = L_u_Corr * L_u_Corr'
                 
                 u_Sigma = quad_form_diag(u_Corr, sigma_sp);
                 // covariance matrix: u_Sigma = diag(sigma_sp) * u_Corr * diag(sigma_sp)
                 