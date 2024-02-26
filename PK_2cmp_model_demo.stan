// Andrew P. Woodward, 2023.
// This Stan program defines a population pharmacokinetic model.
// The setup allows for either IV or PO dosing, for single-dose administration only, with exact solutions for the structural model.
// In this example there is a categorical covariate applied to all the pharmacokinetic parameters.
// The likelihood function for this model allows for left-censoring.

functions {
  // 'PK_model_IV()' contains the structural model, and returns predicted concentrations for a single set of parameters.
  // Model definitions are adapted from Bertrand et. al. 2008.
  real[] PK_model_IV(real [] time, real dose, real Cl, real ClQ, real ka1, real Vc, real Vp, real F){
      // Initialize internal variables.
      int nt;
      real B;
      real A;
      real alpha;
      real beta;
      real result [size(time)];
      nt = size(time);
      beta =  0.5*(((ClQ/Vc)+(ClQ/Vp)+(Cl/Vc))-sqrt((((ClQ/Vc)+(ClQ/Vp)+(Cl/Vc))^2)-(4*((ClQ/Vp)*(Cl/Vc)))));
      alpha = ((ClQ/Vp)*(Cl/Vc))/beta;
      A = (1/Vc)*((alpha-(ClQ/Vp))/(alpha-beta));
      B = (1/Vc)*((beta-(ClQ/Vp))/(beta-alpha));
       // Generate the solution (predicted concentration) of the PK model for each observation.
        for(i in 1:nt){
          result[i] = (((dose)*((A*exp(-alpha*time[i]))+(B*exp(-beta*time[i])))));
        }
    return result;
  }
  // 'PK_model_PO()' contains the structural model, and returns predicted concentrations for a single set of parameters.
  // Model definitions are adapted from Bertrand et. al. 2008.
  real[] PK_model_PO(real [] time, real dose, real Cl, real ClQ, real ka1, real Vc, real Vp, real F){
      // Initialize internal variables.
      int nt;
      real B;
      real A;
      real alpha;
      real beta;
      real result [size(time)];
      nt = size(time);
      beta =  0.5*(((ClQ/Vc)+(ClQ/Vp)+(Cl/Vc))-sqrt((((ClQ/Vc)+(ClQ/Vp)+(Cl/Vc))^2)-(4*((ClQ/Vp)*(Cl/Vc)))));
      alpha = ((ClQ/Vp)*(Cl/Vc))/beta;
      A = (ka1/Vc)*(((ClQ/Vp)-alpha)/((ka1-alpha)*(beta-alpha)));
      B = (ka1/Vc)*(((ClQ/Vp)-beta)/((ka1-beta)*(alpha-beta)));
       // Generate the solution (predicted concentration) of the PK model for each observation.
        for(i in 1:nt){
          result[i] = ((dose*F)*((A*exp(-alpha*time[i]))+(B*exp(-beta*time[i]))-((A+B)*exp(-ka1*time[i]))));
        }
    return result;
  }
  real beta_slope_fun(real Cl, real ClQ, real Vc, real Vp){
    real beta;
      beta =  0.5*(((ClQ/Vc)+(ClQ/Vp)+(Cl/Vc))-sqrt((((ClQ/Vc)+(ClQ/Vp)+(Cl/Vc))^2)-(4*((ClQ/Vp)*(Cl/Vc)))));
    return beta;
  }
  // 'PK_caller()' is a high-level function that is called in the transformed parameters block.
  // This calls the total structural model on each subject-occasion. 
  // It accepts from the data a count of observations J, count of subjects Nsub, and subject integer labels IDvec.
  // The return is an (J) vector containing the model predictions to match the data.
  real [] PK_caller(real [] time, real [] dose, int [] route, real [] Cl, real [] ClQ, real [] ka1, real [] Vc, real [] Vp, real [] F, int J, int Nocc, int [] IDvec, int [] Nobs){
      real raw_predictions[J];  
      int row_start;   
      int row_end;
      row_start = 1;
       // for each subject-occasion;
        for (S in 1:Nocc){
         row_end = (row_start-1) + (Nobs[S]);
         if (route[S] == 0){
         raw_predictions[row_start:row_end] = PK_model_IV(time[row_start:row_end], dose[S], Cl[S], ClQ[S], ka1[S], Vc[S], Vp[S], F[S]);
         }else{
         raw_predictions[row_start:row_end] = PK_model_PO(time[row_start:row_end], dose[S], Cl[S], ClQ[S], ka1[S], Vc[S], Vp[S], F[S]);  
         }
         row_start = row_end+1;
        }
        for (L in 1:J){
          if (time[L] == 0){
            raw_predictions[L] = raw_predictions[L] + 0.00001;
          }
        }
    return raw_predictions;
  }
}

data{
  // no data are actually defined here; we bring the data via RStan.
  int<lower = 1>  J; // the number of time points.
  int<lower = 1>  K; // the total number of non-censored observations.
  int<lower = 1>  P; // the total number of censored observations.
  int<lower = 1>  Nsub; // the number of subjects.
  int<lower = 1>  Nocc; // the number of subject-occasions.
  real<lower = 0> time[J]; // the times of sampling.
  real<lower = 0> conc[K]; // the observed concentrations.
  int<lower = 1>  IDvec[J]; // a vector of subject labels for the times.
  int<lower = 1>  Nobs[Nocc]; // the count of observations for each subject-occasion.
  real<lower = 0> dose[Nocc]; // a vector of doses.
  int<lower = 1>  sub_ind[Nocc]; // a vector of subject labels for the subject-occasions.
  int<lower = 0>  route[Nocc]; // the dose type {0: IV, 1: PO}.
  int<lower = 0>  species[Nocc]; // the species covariate {0: goat, 1: sheep}.
  int<lower = 1>  noncens_vec[K]; // vector index for the non-censored observations.
  int<lower = 1>  cens_vec[P]; // vector index for the censored observations.
  real<lower = 0>  cens_conc[P]; // a vector of upper limits for the censored concentration observations.
}

parameters{
  // the error parameters.
  real<lower = 0> sigma_b;
  // the covariates.
  real species_Cl;
  real species_ClQ;
  real species_ka1;
  real species_Vc;
  real species_Vp;
  real species_F;
  // the population PK parameters.
  // these are on linear-predictor scale, so are unconstrained.
  real Cl_mu;
  real ClQ_mu;
  real ka1_mu;
  real Vc_mu;
  real Vp_mu;
  real F_mu;
  // the between-subject variance parameters.
  // note that their lower zero bound is applied here (so their priors are truncated).
  real <lower = 0> Cl_sd;
  real <lower = 0> ClQ_sd;
  real <lower = 0> ka1_sd;
  real <lower = 0> Vc_sd;
  real <lower = 0> Vp_sd;
  real <lower = 0> F_sd;
  // vectors of subject z-scores (the individual parameters).
  real Cl_ZS  [Nsub];
  real ClQ_ZS [Nsub];
  real ka1_ZS [Nsub];
  real Vc_ZS  [Nsub];
  real Vp_ZS  [Nsub];
  real F_ZS   [Nsub];
}

transformed parameters{
  // assemble vectors of the subject-occasion level parameters.
  real<lower = 0> output_predictions[J];  
  real<lower = 0> Cl_VEC  [Nocc];
  real<lower = 0> ClQ_VEC [Nocc];
  real<lower = 0> ka1_VEC [Nocc];
  real<lower = 0> Vc_VEC  [Nocc];
  real<lower = 0> Vp_VEC  [Nocc];
  real<lower = 0> F_VEC   [Nocc];
  // populate the subject-occasion level vectors using the parameters and covariates.
  for (C in 1:Nocc){
    Cl_VEC[C]  = exp(      (Cl_mu  + (species_Cl* species[C]) + (Cl_sd *Cl_ZS[sub_ind[C]])));
    ClQ_VEC[C] = exp(      (ClQ_mu + (species_ClQ*species[C]) + (ClQ_sd*ClQ_ZS[sub_ind[C]])));
    ka1_VEC[C] = exp(      (ka1_mu + (species_ka1*species[C]) + (ka1_sd*ka1_ZS[sub_ind[C]])));
    Vc_VEC[C]  = exp(      (Vc_mu  + (species_Vc* species[C]) + (Vc_sd *Vc_ZS[sub_ind[C]])));
    Vp_VEC[C]  = exp(      (Vp_mu  + (species_Vp* species[C]) + (Vp_sd *Vp_ZS[sub_ind[C]])));
    F_VEC[C]   = inv_logit((F_mu   + (species_F*  species[C]) + (F_sd  *F_ZS[sub_ind[C]])));
  }
  // generate the model-predicted concentrations by calling to the structural model.
  output_predictions = PK_caller(time, dose, route, Cl_VEC, ClQ_VEC, ka1_VEC, Vc_VEC, Vp_VEC, F_VEC, J, Nocc, IDvec, Nobs);
}

model{
  // define the priors.
  // the error variance (proportional scale).
  sigma_b ~ normal(0,0.5);
  // the species covariates on the pharmacokinetic parameters.
  species_Cl  ~ normal(0,1);
  species_ClQ ~ normal(0,1);
  species_ka1 ~ normal(0,1);
  species_Vc  ~ normal(0,1);
  species_Vp  ~ normal(0,1);
  species_F   ~ normal(0,1);
  // location of the pharmacokinetic parameters (linear predictor scale).
  Cl_mu   ~ normal(1,3);
  ClQ_mu  ~ normal(1,3);
  ka1_mu  ~ normal(-3,2);
  Vc_mu   ~ normal(7,2);
  Vp_mu   ~ normal(7,2);
  F_mu    ~ normal(0,2);
  // scale of the pharmacokinetic parameters (linear predictor scale).
  Cl_sd   ~ normal(0.5,0.5);
  ClQ_sd  ~ normal(0.5,0.5);
  ka1_sd  ~ normal(0.5,0.5);
  Vc_sd   ~ normal(0.5,0.5);
  Vp_sd   ~ normal(0.5,0.5);
  F_sd    ~ normal(1,1);
  // the z-scores for the individual parameters (for the non-centered parameterization).
  Cl_ZS   ~ normal(0,1);
  ClQ_ZS  ~ normal(0,1);
  ka1_ZS  ~ normal(0,1);
  Vc_ZS   ~ normal(0,1);
  Vp_ZS   ~ normal(0,1);
  F_ZS    ~ normal(0,1);
  // define the likelihood.
  // the likelihood for the non-censored observations, generated from the predicted outputs and the error model (proportional-normal in this case, for simplicity).  
  for (M in 1:K){
    conc[M] ~ normal(output_predictions[noncens_vec[M]], (sigma_b*output_predictions[noncens_vec[M]]));  
  }
  // the likelihood for the censored observations, from the normal cumulative distribution function. The censoring limit may be different for each observation.  
  for (M in 1:P){
    target += normal_lcdf(cens_conc[M] | output_predictions[cens_vec[M]], (sigma_b*output_predictions[cens_vec[M]]));
  }
}

generated quantities{
  // add the terminal half-life to the reported parameters; alternately, this could be done from the samples during postprocessing.
  real beta_VEC [Nocc];
  real beta_mu_goat;
  real beta_mu_sheep;
  for (C in 1:Nocc){
    beta_VEC[C] = beta_slope_fun(Cl_VEC[C],ClQ_VEC[C],Vc_VEC[C],Vp_VEC[C]);
  }
  beta_mu_goat   = beta_slope_fun(exp(Cl_mu),exp(ClQ_mu),exp(Vc_mu),exp(Vp_mu));
  beta_mu_sheep  = beta_slope_fun(exp(Cl_mu+species_Cl),exp(ClQ_mu+species_ClQ),exp(Vc_mu+species_Vc),exp(Vp_mu+species_Vp));
}
