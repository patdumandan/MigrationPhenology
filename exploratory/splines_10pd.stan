data { 

  int N; //rows of observations 
  int num_basis1; //no. of basis (order-1) 
  real pd_10[N]; // 10% passage date (log-transformed)
  vector[N] years; //years
  matrix[num_basis1, N] B1; //matrix of coefficients of splines(rows), length of years
 
} 
 
parameters { 
  
  row_vector[num_basis1] b_raw; // smooth term for year
  real a0; //intercept
  real<lower=0> sigma; //process error
  real <lower=0> tau; // error term for noncentered parameterization of spline coefficients (year)

} 
 
transformed parameters { 
 
  row_vector[num_basis1] b; //noncentered parameters of splines for year coef
  vector[N] Y_hat; 
  
  b = b_raw*tau;
  
 Y_hat=a0*years + to_vector(b*B1); //predictions estimated as parameters
 }

model { 
  
  a0~normal(0,1);
  b_raw ~ normal(0,1); 
  tau ~ cauchy(0, 1);  
  sigma ~ cauchy(0, 1); 
  
  //likelihood
  pd_10~ normal(Y_hat, sigma);
  
  }
  
generated quantities {
  real y_pred[N];
  
  for (n in 1:N)
  {
  y_pred[n]= normal_rng(Y_hat[n], sigma);
  }
  
}
  