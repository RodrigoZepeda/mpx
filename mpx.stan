//Model for Monkepox data
//----------------------------------------
//Author: Rodrigo Zepeda-Tello rzepeda17@gmail.com
//Adapted from: https://www.medrxiv.org/content/10.1101/2022.08.17.22278897v1.full.pdf
functions {
  vector sir_pair(real t, vector y, vector theta, int N){
      //State      y = (S, I, R, PSS, PSI, PII, PSR, PIR, PRR, T); length(y) = 10
      //Parameter theta = (B, mu, rho, sigma, delta, h, phi); length(theta = 7)
      //Integer N > 0 represents the total population (normalizing constant)
      
      //State
      //--------------------------------------------------
      vector[10] dy;
      real S   = y[1];  //Susceptible
      real I   = y[2];  //Infected
      real R   = y[3];  //Recovered
      real PSS = y[4];  //Susceptible-Susceptible pairing
      real PSI = y[5];  //Susceptible-Infected pairing
      real PII = y[6];  //Infected-Infected pairing
      real PSR = y[7];  //Susceptible-Recovered pairing
      real PIR = y[8];  //Infected-Recovered pairing
      real PRR = y[9];  //Recovered-Recovered pairing
      real T   = y[10]; //Total infected
      
      //Parameters
      //--------------------------------------------------
      real B     = theta[1]; //Birth-rate 
      real mu    = theta[2]; //Death-rate 
      real rho   = theta[3]; //Pair formation rate
      real sigma = theta[4]; //Pair dissolution rate 
      real delta = theta[5]; //Infection recovery rate
      real h     = theta[6]; //Probability of transmission
      real phi   = theta[7]; //Contact rate within partnership
      
      //Model
      //--------------------------------------------------
      dy[1]  = B - (mu + rho)*S + (sigma + mu)*(2.0*PSS + PSI + PSR);                      //dS/dt
      dy[2]  = -(mu + rho + delta)*I + (sigma + mu)*(2.0*PII + PSI + PIR);                 //dI/dt
      dy[3]  = -(mu + rho)*R + delta*I + (sigma + mu)*(2.0*PRR + PSR + PIR);               //dR/dt
      dy[4]  = rho/2.0*S^2/N - (sigma + 2.0*mu)*PSS;                                       //dPSS/dt
      dy[5]  = rho*(1.0 - h)*S*I/N - (sigma + phi*h + 2.0*mu + delta)*PSI;                 //dPSI/dt
      dy[6]  = rho/2.0*I^2/N + rho*h*S*I/N + phi*h*PSI - (sigma + 2.0*mu + 2.0*delta)*PII; //dPII/dt
      dy[7]  = delta*PSI + rho*S*R/N - (sigma + 2.0*mu)*PSR;                               //dPSR/dt
      dy[8]  = rho*I*R/N + delta*PII - (sigma + 2.0*mu + delta)*PIR;                       //dPIR/dt
      dy[9]  = delta*PIR + rho/2.0*R^2.0/N - (sigma + 2.0*mu)*PRR;                         //dPRR/dt
      dy[10] = rho*h*S*I/N + phi*h*PSI - (mu + delta)*T - delta*PII;                        //dT/dt
      
      return dy;
  }
}

data {
  int<lower=1> Nweeks;       //Number of weeks measured by model
  int<lower=0> N;            //Total size of the population
  int cases[Nweeks];         //Total mpx cases registered
  vector[10] initial_values; //Initial values
}

transformed data {
  array[Nweeks] real<lower = 0> weeks; //Vector of sequential weeks
  real initial_time  = 0;              //Starting time for model
  int model_dim = 10;                  //Total number of variables in ODE
  int infected_index = 10;             //Index for the infectious
  for (i in 1:Nweeks){
    weeks[i] = i;
  }
}

parameters {
  real<lower=0> B;          //Birth-rate 
  real<lower=0> mu;         //Death-rate 
  real<lower=0> rho_inv;    //Pair formation inverse rate (days)
  real<lower=0> sigma_inv;  //Pair dissolution inverse rate (days)
  real<lower=0> delta;      //Infection recovery rate
  real<lower=0> h;          //Probability of transmission
  real<lower=0> phi;        //Contact rate within partnership
  real<lower=0> varphi_inv; //Noise associated to measurement (negative binomial second parameter)
}

transformed parameters{
  vector[model_dim] mpx[Nweeks];
  real varphi = 1.0 / varphi_inv;
  real rho    = 1.0 / rho_inv;
  real sigma  = 1.0 / sigma_inv;
  {
    vector[7] theta;
    theta[1] = B;
    theta[2] = mu;
    theta[3] = rho;
    theta[4] = sigma;
    theta[5] = delta;
    theta[6] = h;
    theta[7] = phi;
    
    mpx = ode_rk45(sir_pair, initial_values, initial_time, weeks, theta, N);
  }
}

model {
  //priors
  B          ~ normal(0.0, 0.001);
  mu         ~ normal(1.0 / 18250.0, 0.001);
  rho_inv    ~ exponential(15.0 / 7.0);        //Pair formation rate
  sigma_inv  ~ exponential(1.0 / 7.0);         //Pair dissolution rate
  delta      ~ normal(7.0 / 30.0, 0.001);      //Fixed parameter from the literature
  h          ~ normal(0.9, 0.001);
  phi        ~ normal(1.0, 0.001);
  varphi_inv ~ exponential(5.0);
  
  //sampling distribution
  for (week in 1:Nweeks) 
    cases[week] ~ neg_binomial_2(mpx[week][infected_index], phi);
}

