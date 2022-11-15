data{
  
  int N; //length of observations
  int Nsp; //number of species
  //vector[N] div_peak;
  int div_peak[N]; // shannon index per year
  vector[N] years; // standardized years
  
}

parameters {
  
  real alpha; //intercept
  real <lower=-0.81, upper=0.85> bkpt; // bkpt year est
  real slope1; //pre-bkpoint slope
  real slope2; //post-bkpoint slope
  real phi; //process error

}

transformed parameters{
  
  real <lower=0> y_mu[N];
  
  for (i in 1:N) {
    
    if (years[i]<bkpt){
      
      y_mu[i]=alpha+slope1*years[i]-bkpt;
    }
    
    else {
      y_mu[i]=alpha+(slope2*years[i])-bkpt;
    }
  }
}

model {
  
  alpha ~normal(0,10);
  slope1~ normal(0,10);
  slope2~ normal(0,10);
  bkpt~ normal(0,10);
  phi ~cauchy(0,2.5);
  
  div_peak~ poisson(y_mu);
  //if peak diversity date
  //for (j in 1:N) {
  //div_peak ~ normal(y_mu, phi);
  //}
}
