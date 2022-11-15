data{
  
  int N; //length of observations
  real div_index[N]; // shannon index per year
  vector[N] years; // standardized years
  vector[N] bird_time; //peak migration dates per species per year
  vector[N] bird_count;// peak numbers per species per year 
  
}

parameters {
  
  real alpha; //intercept
  real b0; //year effect
  real time_eff; // slope; timing effect vary per species
  real count_eff; //slope2; intensity effect vary per species
  real phi; //process error

}

model {
  
  alpha ~normal(0,10);
  b0~ normal(0,10);
  time_eff~ normal(0,10);
  count_eff~ normal(0,10);
  phi ~normal(0,10);
  
  div_index ~ normal(alpha+ years*b0+bird_time*time_eff+bird_count*count_eff, phi);
  
}