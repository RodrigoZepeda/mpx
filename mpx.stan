//Model for Monkepox data
//----------------------------------------
//Author: Rodrigo Zepeda-Tello rzepeda17@gmail.com
//Adapted from: https://www.medrxiv.org/content/10.1101/2022.08.17.22278897v1.full.pdf
functions {
  vector sir_pair(real t, vector y, vector theta, int N) {
    //State      y = (S, I, R, PSS, PSI, PII, PSR, PIR, PRR, T); length(y) = 10
    //Parameter theta = (B, mu, rho, sigma, delta, h, phi); length(theta = 7)
    //Integer N > 0 represents the total population (normalizing constant)
    
    //State
    //--------------------------------------------------
    vector[10] dy;
    real S = y[1];    //Susceptible
    real I = y[2];    //Infected
    real R = y[3];    //Recovered
    real PSS = y[4];  //Susceptible-Susceptible pairing
    real PSI = y[5];  //Susceptible-Infected pairing
    real PII = y[6];  //Infected-Infected pairing
    real PSR = y[7];  //Susceptible-Recovered pairing
    real PIR = y[8];  //Infected-Recovered pairing
    real PRR = y[9];  //Recovered-Recovered pairing
    real Inc = y[10]; //Incidence 
    
    //Parameters
    //--------------------------------------------------
    real rho   = theta[1]; //Pair formation rate
    real sigma = theta[2]; //Pair dissolution rate 
    real delta = theta[3]; //Infection recovery rate
    real h     = theta[4]; //Probability of transmission
    real phi   = theta[5]; //Contact rate within partnership
    
    //Model
    //--------------------------------------------------
    dy[1]  = -rho*S + sigma*(2.0*PSS + PSI + PSR);        //dS/dt
    dy[2]  = -(rho + delta)*I + sigma*(2.0*PII + PSI + PIR);   //dI/dt
    dy[3]  = -rho*R + delta*I + sigma*(2.0*PRR + PSR + PIR); //dR/dt
    dy[4]  = rho/2.0*S^2/N - sigma*PSS;                         //dPSS/dt
    dy[5]  = rho*(1.0 - h)*S*I/N - (sigma + phi*h + delta)*PSI;   //dPSI/dt
    dy[6]  = rho/2.0*I^2/N + rho*h*S*I/N +phi*h*PSI - (sigma + 2.0*delta)*PII; //dPII/dt
    dy[7]  = delta*PSI + rho*S*R/N - sigma*PSR;           //dPSR/dt
    dy[8]  = rho*I*R/N + delta*PII - (sigma + delta)*PIR; //dPIR/dt
    dy[9]  = delta*PIR + rho/2.0*R^2.0/N - sigma*PRR;     //dPRR/dt
    dy[10] = rho*h*S*I/N + phi*h*PSI;   //dT
    
    return dy;
  }
}

data {
  int<lower=1> Nweeks;           //Number of weeks measured by model
  int<lower=0> N;                //Total size of the population
  array[Nweeks] int cases;       //Total mpx cases registered
  vector[10] initial_values;
}

transformed data {
  array[Nweeks] real<lower=0> weeks; //Vector of sequential weeks
  real initial_time = 0;             //Starting time for model
  int model_dim = 10;                //Total number of variables in ODE
  for (i in 1 : Nweeks) {
    weeks[i] = i;
  }
}

parameters {
  real<lower=0> rho;                //Pair formation inverse rate (days)
  real<lower=0> sigma;              //Pair dissolution inverse rate (days)
  real<lower=0> delta;              //Infection recovery rate
  real<lower=0,upper=1> h;          //Probability of transmission
  real<lower=0> phi;                //Contact rate within partnership
  real<lower=0> varphi;             //Noise associated to measurement (negative binomial second parameter)
  real<lower=0,upper=1> p_detection;//Probability of detection
}

transformed parameters {
  array[Nweeks] vector[model_dim] mpx; //ODE model cases
  vector[Nweeks] incidence;
  {
    vector[5] theta;
    theta[1] = rho;
    theta[2] = sigma;
    theta[3] = delta;
    theta[4] = h;
    theta[5] = phi;
    
    mpx = ode_rk45_tol(sir_pair, initial_values, initial_time, weeks, 1.0e-8, 1.0e-8, 10000,  theta, N);
  }
  
  incidence[1] = mpx[1][10]*p_detection;
  for (week in 2:Nweeks)
    incidence[week] = (mpx[week][10] - mpx[week - 1][10])*p_detection + 1.0;
  
}

model {
  rho         ~ normal(0.0, 0.1);   //Pair formation rate (average 15.0 days = 15/7 weeks)
  sigma       ~ normal(0.0, 0.1); //Pair dissolution rate = 1 day
  delta       ~ normal(0.0, 0.01);      //Fixed parameter from the literature
  h           ~ normal(0.9, 0.01);
  phi         ~ normal(1.0, 0.01);
  p_detection ~ beta(1,2);
  phi         ~ normal(0.0, 0.1);
  
  //sampling distribution
  cases ~ neg_binomial_2(incidence, varphi);
}


