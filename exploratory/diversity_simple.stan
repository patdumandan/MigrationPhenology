data{
  
  int N; //length of observations
  int Nsp; //number of species
  real div_index[N]; // shannon index per year
  vector[N] years; // standardized years
 // real sp_time[Nsp]; //peak migration dates per species per year
  //real sp_count[Nsp];// peak numbers per species per year 
  //int spcode[N]; //code for each population
  
}

parameters {
  
  real alpha; //intercept
  real b0; //year effect
 // real time_eff[Nsp]; // slope; timing effect vary per species
//  real count_eff[Nsp]; //slope2; intensity effect vary per species
  real phi; //process error
  //real sigma_sp;// species error

}

model {
  
  alpha ~normal(0,10);
  b0~ normal(0,10);
  //time_eff~ normal(0,10);
  //count_eff~ normal(0,10);
  phi ~normal(0,10);
  //sigma_sp~normal(0,10);
  
  div_index ~ normal(alpha+ years*b0, phi);
  
}