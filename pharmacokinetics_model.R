# AP. Woodward, University of Canberra, 2023.
# This R script includes simulation of a pharmacokinetic experiment and analysis of the resulting data, using a population pharmacokinetic model implemented in Stan (https://doi.org/10.18637/jss.v076.i01).
# The objective of this implementation is to illustrate principles of Bayesian statistics relevant to veterinary pharmacology, especially the specification of priors and post-processing of the model to generate expressive graphics.
#     This R specification and the accompanying Stan code could be easily adapted to studies of common designs.
# This implementation is intended to function comparably to off-the-shelf population pharmacokinetics software popular in veterinary pharmacology practice.
#     Potential benefits highlighted here are stability of the model even with many parameters, and easy generation of uncertainty statements for secondary parameters and functions of primary secondary parameters.
#     Key features of the model implemented here are simultaneous analysis of different routes of administration, analysis of a categorical covariate, and censored responses.

# Load the required packages.
library(deSolve)
library(ggplot2)
library(PKconverter)
library(viridis)
library(rlist)
library(brms)
library(ggpubr)
library(rstan)
library(ggdist)
library(matrixStats)
library(distributional)
library(ggh4x)

# Define the exact solution for a two-compartment model with intravenous bolus administration (specification is adapted from: https://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf).
pk_2cmp_mod_IV <- function(time,dose,Cl,Q,V1,V2){
  k12 = Q/V1
  k21 = Q/V2
  k10 = Cl/V1
  beta  = 0.5*(k12+k21+k10-sqrt((k12+k21+k10)^2-4*k21*k10)) 
  alpha = (k21*k10)/beta 
  A = (1/V1)*((alpha-k21)/(alpha-beta)) 
  B = (1/V1)*((beta-k21)/(beta-alpha))
  Ct = dose*((A*exp(-alpha*time))+(B*exp(-beta*time)))
  return(list(cbind(Ct, time), c(A = A, B = B, alpha = alpha, beta = beta, k12 = k12, k21 = k21, k10 = k10)))
}

# Define the exact solution for a two-compartment model with first-order absorption with bioavailability (specification is adapted from: https://lixoft.com/wp-content/uploads/2016/03/PKPDlibrary.pdf).
pk_2cmp_mod_PO <- function(time,dose,Cl,Q,V1,V2,ka,F){
  k12 = Q/V1
  k21 = Q/V2
  k10 = Cl/V1
  beta  = 0.5*(k12+k21+k10-sqrt((k12+k21+k10)^2-4*k21*k10)) 
  alpha = (k21*k10)/beta 
  A = (ka/V1)*((k21-alpha)/((ka-alpha)*(beta-alpha))) 
  B = (ka/V1)*((k21-beta)/ ((ka-beta)*(alpha-beta)))
  Ct = (dose*F)*(A*exp(-alpha*time)+B*exp(-beta*time)-(A+B)*exp(-ka*time))
  return(list(cbind(Ct, time), c(A = A, B = B, alpha = alpha, beta = beta, k12 = k12, k21 = k21, k10 = k10)))
}

# Simulate a 16-subject crossover PK experiment, where subjects of two species {sheep,goat} receive the same dose by each of IV bolus and first-order PO with bioavailability.
# The study design and specification of the parameters is inspired by (https://doi.org/10.3389/fvets.2020.00586). 
# This simulation is intended to be fairly representative of descriptive pharmacokinetic studies in veterinary applications.
#     The population parameters are log-normally distributed (logit-normal for F).
#     In this case the population parameters are uncorrelated, but their correlation can in principle be modelled (it is often reasonable to do so).
#     Units for the variables are ng/kg (dose), ng/mL (concentration), and min (time).
#     Units for the parameters are therefore mL/kg/min (clearance) and mL/kg (volume).
#     The dose administered is 10000ng/kg (0.01mg/kg).
n_subjects <- 16
subject_PK_IV_list <- vector(length = n_subjects, mode = 'list')
subject_PK_PO_list <- vector(length = n_subjects, mode = 'list')
pk_ind_theta <- data.frame(subject = factor(seq(1,16,1)), species = character(n_subjects), clearance = numeric(n_subjects), Q = numeric(n_subjects), VC = numeric(n_subjects), VP = numeric(n_subjects), ka = numeric(n_subjects), F = numeric(n_subjects))
pk_ind_theta$species <- rep(c('sheep','goat'), each = 8)
pk_ind_theta$clearance[pk_ind_theta$species == 'sheep'] <- exp(rnorm(n_subjects/2, 1.2, 0.3))
pk_ind_theta$clearance[pk_ind_theta$species == 'goat']  <- exp(rnorm(n_subjects/2, 2, 0.3))
pk_ind_theta$Q  <- exp(rnorm(n_subjects, 1.4, 0.2))
pk_ind_theta$VC <- exp(rnorm(n_subjects, 7.2, 0.4))
pk_ind_theta$VP <- exp(rnorm(n_subjects, 7, 0.3))
pk_ind_theta$ka <- exp(rnorm(n_subjects, -4, 0.4))
pk_ind_theta$F  <- plogis(rnorm(n_subjects, 3, 1))
pk_time_seq <- c(0,15,30,45,60,90,120,180,360,720,1440,2160,2880,3600,4320,5760,7200)
  for(i in 1:n_subjects){
    pred_IV_2cmp <- as.data.frame((pk_2cmp_mod_IV(pk_time_seq,1000000,pk_ind_theta$clearance[i],pk_ind_theta$Q[i],pk_ind_theta$VC[i],pk_ind_theta$VP[i]))[[1]])
    pred_IV_2cmp$subject <- i
    pred_IV_2cmp$route   <- 'IV'
    pred_IV_2cmp$species <- pk_ind_theta$species[i]
    subject_PK_IV_list[[i]] <- pred_IV_2cmp
    pred_PO_2cmp <- as.data.frame((pk_2cmp_mod_PO(pk_time_seq,1000000,pk_ind_theta$clearance[i],pk_ind_theta$Q[i],pk_ind_theta$VC[i],pk_ind_theta$VP[i],pk_ind_theta$ka[i],pk_ind_theta$F[i]))[[1]])
    pred_PO_2cmp$subject <- i
    pred_PO_2cmp$route   <- 'PO'
    pred_PO_2cmp$species <- pk_ind_theta$species[i]
    subject_PK_PO_list[[i]] <- pred_PO_2cmp
  }
subject_PK_data <- rbind(list.rbind(subject_PK_IV_list),list.rbind(subject_PK_PO_list))
subject_PK_data$Ct <- (exp(log(subject_PK_data$Ct)+rnorm(length(subject_PK_data$Ct),0,0.25)))
subject_PK_data$CENS <- 0
subject_PK_data$CENS[subject_PK_data$Ct<5] <- -1
subject_PK_data$Ct[subject_PK_data$CENS == -1] <- 5
subject_PK_data$subject <- factor(subject_PK_data$subject)
subject_PK_data$route <- factor(subject_PK_data$route)
subject_PK_data$species <- factor(subject_PK_data$species)
subject_PK_data$time_log1 <- log(subject_PK_data$time+1)
subject_PK_data$subject_route <- subject_PK_data$route:subject_PK_data$subject

# Generate a HGAM-based NCA to visualize the subject-level and population-level trajectories (https://doi.org/10.1101/2023.07.13.548803).
#     For this case we'll forego any parameter estimation from the HGAM and instead implement that via the parametric model; the HGAM will be used only for visualization.
NCA_HGAM_mod <- brm(log(Ct)|cens(CENS) ~ 1 + route*species + s(time_log1, by = interaction(route,species)) + s(time_log1, subject, bs = 'fs') + s(time_log1, subject_route, bs = 'fs'), data = subject_PK_data, chains = 4, cores = 4, iter = 4000, control = list(adapt_delta = 0.95, max_treedepth = 12))
summary(NCA_HGAM_mod)
NCA_HGAM_new <- expand.grid(time_log1 = log(seq(0,1500,1)+1), route = c('IV','PO'), species = c('sheep','goat'), subject = NA, subject_route = NA)
NCA_HGAM_new$time <- seq(0,1500,1)
NCA_HGAM_pop <- cbind(NCA_HGAM_new, exp(fitted(NCA_HGAM_mod, newdata = NCA_HGAM_new)))
NCA_HGAM_linear <- ggplot(data = NCA_HGAM_pop, aes(x = time, y = Estimate)) + geom_line(aes(color = species)) + geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = species), alpha = 0.1) + facet_wrap(~route) + scale_y_log10() + geom_point(data = subject_PK_data, aes(x = time, y = Ct, color = species, shape = as.factor(CENS)), inherit.aes = FALSE, fill = 'white', alpha = 0.5) + scale_shape_manual(values = c(21,19), guide = 'none') + theme_bw() + coord_cartesian(xlim = c(0,1500), ylim = c(0.01,100)) + theme(legend.position = 'top', legend.title = element_blank()) + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_fill_viridis(discrete = TRUE, end = 0.8) + geom_line(data = subject_PK_data, aes(x = time, y = Ct, group = subject, color = species), alpha = 0.1) + xlab('Time (min)') + ylab('Concentration (mg/L)')
NCA_HGAM_log    <- ggplot(data = NCA_HGAM_pop, aes(x = time, y = Estimate)) + geom_line(aes(color = species)) + geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5, fill = species), alpha = 0.1) + facet_wrap(~route) + scale_y_log10() + geom_point(data = subject_PK_data, aes(x = time, y = Ct, color = species, shape = as.factor(CENS)), inherit.aes = FALSE, fill = 'white', alpha = 0.5) + scale_shape_manual(values = c(21,19), guide = 'none') + theme_bw() + coord_cartesian(xlim = c(1,1500), ylim = c(0.01,100)) + theme(legend.position = 'top', legend.title = element_blank()) + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_fill_viridis(discrete = TRUE, end = 0.8) + geom_line(data = subject_PK_data, aes(x = time, y = Ct, group = subject, color = species), alpha = 0.1) + xlab('Time (min)') + ylab('Concentration (mg/L)') + scale_x_log10()    
NCA_HGAM_plots  <- ggarrange(NCA_HGAM_linear, NCA_HGAM_log, nrow = 2, common.legend = TRUE)
ggsave(NCA_HGAM_plots, file = 'NCA_HGAM_plots.svg', units = 'mm', width = 200, height = 220)

# Define a population pharmacokinetic model.
#     For this case we'll presume that a two-compartment model is a reasonable selection (in the real case this wouldn't be known).
#     Take the prior distributions for the PK parameters as log-normal (logit-normal for F); this is a typical choice, but is a key assumption.
#     The model is defined in Stan (https://doi.org/10.18637/jss.v076.i01); then the posterior draws are returned here for post-processing.
subject_PK_cov <- subject_PK_data[which(!duplicated(subject_PK_data[c(3,4)])),]
subject_PK_cov$route_num <- 0
subject_PK_cov$route_num[subject_PK_cov$route == 'PO'] <- 1
subject_PK_cov$species_num <- 0
subject_PK_cov$species_num[subject_PK_cov$species == 'sheep'] <- 1

# Set the initial values for the MCMC chains (this seems generally necessary for nonlinear models).
#     In this case the initial values are fairly plausible under the priors, and are randomly perturbed to add some diversity to the chains.
PK_stan_inits <- function(){
    list(sigma_b = runif(1,0.1,1), 
         species_Cl = runif(1,-1,1), 
         species_ClQ = runif(1,-1,1), 
         species_ka1 = runif(1,-1,1), 
         species_Vc = runif(1,-1,1), 
         species_Vp = runif(1,-1,1), 
         species_F = runif(1,-1,1), 
         Cl_mu = runif(1,-1,2), 
         ClQ_mu = runif(1,-1,2), 
         ka1_mu = runif(1,-5,-1), 
         Vc_mu = runif(1,6,9), 
         Vp_mu = runif(1,6,9), 
         F_mu = runif(1,-2,2), 
         Cl_sd = runif(1,0.2,1), 
         ClQ_sd = runif(1,0.2,1), 
         ka1_sd = runif(1,0.2,1), 
         Vc_sd = runif(1,0.2,1), 
         Vp_sd = runif(1,0.2,1), 
         F_sd = runif(1,0.2,1), 
         Cl_ZS = runif(16,-2,2), 
         ClQ_ZS = runif(16,-2,2), 
         ka1_ZS = runif(16,-2,2), 
         Vc_ZS = runif(16,-2,2),
         Vp_ZS = runif(16,-2,2), 
         F_ZS = runif(16,-2,2))
}

# Define and fit the model.
#     The simulated data includes time-zero observations for IV phases, so those are removed.
subject_PK_data_nozero <- subject_PK_data[!((subject_PK_data$route == 'IV') & (subject_PK_data$time == 0)),]
PK_stan_list <- list(J = length(subject_PK_data_nozero$time), Q = length(subject_PK_data_nozero$Ct), K = length(which(!(subject_PK_data_nozero$CENS == -1))), P = length(which((subject_PK_data_nozero$CENS == -1))), Nsub = length(unique(subject_PK_data_nozero$subject)), Nocc = length(unique(subject_PK_data$subject_route)), time = subject_PK_data_nozero$time, conc = subject_PK_data_nozero$Ct[subject_PK_data_nozero$CENS == 0], IDvec = as.integer(subject_PK_data_nozero$subject), Nobs = as.integer(table(subject_PK_data_nozero$subject_route)), dose = rep(1000000,dim(subject_PK_cov)[1]), sub_ind = as.integer(subject_PK_cov$subject), route = subject_PK_cov$route_num, species = subject_PK_cov$species_num, noncens_vec = (which(!(subject_PK_data_nozero$CENS == -1))), cens_vec = (which((subject_PK_data_nozero$CENS == -1))), cens_conc = subject_PK_data_nozero$Ct[(subject_PK_data_nozero$CENS == -1)])
PK_2cmp_test_mod2 <- stan(file = 'PK_2cmp_model_demo.stan', data = PK_stan_list, chains = 4, cores = 4, iter = 4000, init = PK_stan_inits, control = list(adapt_delta = 0.95))
summary(PK_2cmp_test_mod2)

# Generate visual descriptions of the priors.
PK_prior_parms <- data.frame(parameter = character(12), distribution = character(12))
PK_prior_parms$parameter <- c('Cl_mu','ClQ_mu','Vc_mu','Vp_mu','ka1_mu','F_mu','Cl_sd','ClQ_sd','Vc_sd','Vp_sd','ka1_sd','F_sd')
PK_prior_parms$distribution <- c(dist_normal(1,3),
                                 dist = dist_normal(1,3),
                                 dist = dist_normal(7,2),
                                 dist = dist_normal(7,2),
                                 dist = dist_normal(-3,2),
                                 dist = dist_normal(0,2),
                                 dist_truncated(dist = dist_normal(0.5,0.5), lower = 0),
                                 dist_truncated(dist = dist_normal(0.5,0.5), lower = 0),
                                 dist_truncated(dist = dist_normal(0.5,0.5), lower = 0),
                                 dist_truncated(dist = dist_normal(0.5,0.5), lower = 0),
                                 dist_truncated(dist = dist_normal(0.5,0.5), lower = 0),
                                 dist_truncated(dist = dist_normal(1,1),     lower = 0))
PK_prior_parms$parameter <- factor(PK_prior_parms$parameter, levels = c('Cl_mu','Cl_sd','ClQ_mu','ClQ_sd','Vc_mu','Vc_sd','Vp_mu','Vp_sd','ka1_mu','ka1_sd','F_mu','F_sd'))
PK_prior_labels = c(Cl_mu = 'Population median plasma clearance (mL/kg/min)',ClQ_mu = 'Population median intercompartmental clearance (mL/kg/min)', Vc_mu = 'Population median central volume (mL/kg)', Vp_mu = 'Population median peripheral volume (mL/kg)', ka1_mu = 'Population median absorption rate (1/min)', F_mu = 'Population median bioavailability (%)',
                    Cl_sd = 'Population SD log(plasma clearance)',ClQ_sd = 'Population SD log(intercompartmental clearance)', Vc_sd = 'Population SD log(central volume)', Vp_sd = 'Population SD log(peripheral volume)', ka1_sd = 'Population SD log(absorption rate)', F_sd = 'Population SD logit(bioavailability)')
PK_prior_plots <- ggplot(data = PK_prior_parms, aes(xdist = distribution)) + stat_halfeye(alpha = 0.7, .width = c(.5,0.9)) + facet_wrap(~parameter, scales = 'free', ncol = 2, labeller = labeller(parameter = PK_prior_labels)) + theme_bw() + scale_y_continuous(breaks = NULL) + xlab('') + ylab('') + theme(panel.grid.minor = element_blank()) +
  facetted_pos_scales(y = list(parameter == 'Cl_mu' ~ scale_y_continuous(limits = c(0,0.15), breaks = NULL), parameter == 'ClQ_mu' ~ scale_y_continuous(limits = c(0,0.15), breaks = NULL), parameter == 'F_mu' ~ scale_y_continuous(limits = c(0,0.2), breaks = NULL), parameter == 'ka1_mu' ~ scale_y_continuous(limits = c(0,0.2), breaks = NULL), parameter == 'Vc_mu' ~ scale_y_continuous(limits = c(0,0.2), breaks = NULL), parameter == 'Vp_mu' ~ scale_y_continuous(limits = c(0,0.2), breaks = NULL)), 
                      x = list(parameter == 'Cl_mu' ~ scale_x_continuous(breaks = log(c(0.1,1,10,100,1000)), labels = c(0.1,1,10,100,1000)), parameter == 'ClQ_mu' ~ scale_x_continuous(breaks = log(c(0.1,1,10,100,1000)), labels = c(0.1,1,10,100,1000)), parameter == 'F_mu' ~ scale_x_continuous(breaks = qlogis(c(0.01,0.1,0.25,0.50,0.75,0.9,0.99)), labels = c(0.01,0.1,0.25,0.50,0.75,0.9,0.99)), parameter == 'ka1_mu' ~ scale_x_continuous(breaks = log(c(0.001,0.01,0.1,1)), labels = c(0.001,0.01,0.1,1)), parameter == 'Vc_mu' ~ scale_x_continuous(breaks = log(c(100,1000,10000)), labels = c(100,1000,10000)), parameter == 'Vp_mu' ~ scale_x_continuous(breaks = log(c(100,1000,10000)), labels = c(100,1000,10000))))
ggsave(PK_prior_plots, file = 'PK_prior_plots.svg', units = 'mm', height = 300, width = 240)

# Obtain the posterior samples of the population parameters including the species effects.
PK_2cmp_pop_draws <- extract(PK_2cmp_test_mod2, pars = c('Cl_mu','ClQ_mu','ka1_mu','Vc_mu','Vp_mu','F_mu','beta_mu_goat','beta_mu_sheep','species_Cl','species_ClQ','species_ka1','species_Vc','species_Vp','species_F','Cl_sd','ClQ_sd','ka1_sd','Vc_sd','Vp_sd','F_sd'))
PK_2cmp_pop_theta <- rbind(data.frame(parameter = 'Cl',    sample = PK_2cmp_pop_draws$Cl_mu,                                  species = 'goat'),
                           data.frame(parameter = 'ClQ',   sample = PK_2cmp_pop_draws$ClQ_mu,                                 species = 'goat'),
                           data.frame(parameter = 'ka1',   sample = PK_2cmp_pop_draws$ka1_mu,                                 species = 'goat'),
                           data.frame(parameter = 'Vc',    sample = PK_2cmp_pop_draws$Vc_mu,                                  species = 'goat'),
                           data.frame(parameter = 'Vp',    sample = PK_2cmp_pop_draws$Vp_mu,                                  species = 'goat'),
                           data.frame(parameter = 'F',     sample = PK_2cmp_pop_draws$F_mu,                                   species = 'goat'),
                           data.frame(parameter = 'Thalf', sample = (log(2))/PK_2cmp_pop_draws$beta_mu_goat,                  species = 'goat'),
                           data.frame(parameter = 'AUC',   sample = rep(NA,length(PK_2cmp_pop_draws$Cl_mu)),                  species = 'goat'))
PK_2cmp_pop_theta <- rbind(data.frame(parameter = 'Cl',    sample = PK_2cmp_pop_draws$Cl_mu  + PK_2cmp_pop_draws$species_Cl,  species = 'sheep'),
                           data.frame(parameter = 'ClQ',   sample = PK_2cmp_pop_draws$ClQ_mu + PK_2cmp_pop_draws$species_ClQ, species = 'sheep'),
                           data.frame(parameter = 'ka1',   sample = PK_2cmp_pop_draws$ka1_mu + PK_2cmp_pop_draws$species_ka1, species = 'sheep'),
                           data.frame(parameter = 'Vc',    sample = PK_2cmp_pop_draws$Vc_mu  + PK_2cmp_pop_draws$species_Vc,  species = 'sheep'),
                           data.frame(parameter = 'Vp',    sample = PK_2cmp_pop_draws$Vp_mu  + PK_2cmp_pop_draws$species_Vp,  species = 'sheep'),
                           data.frame(parameter = 'F',     sample = PK_2cmp_pop_draws$F_mu   + PK_2cmp_pop_draws$species_F,   species = 'sheep'),
                           data.frame(parameter = 'Thalf', sample = (log(2))/PK_2cmp_pop_draws$beta_mu_sheep,                 species = 'sheep'),
                           data.frame(parameter = 'AUC',   sample = rep(NA,length(PK_2cmp_pop_draws$Cl_mu)),                  species = 'sheep'),
                           PK_2cmp_pop_theta)
PK_2cmp_pop_theta <- rbind(data.frame(parameter = 'Cl_sd',    sample = PK_2cmp_pop_draws$Cl_sd,                               species = 'SD'),
                           data.frame(parameter = 'ClQ_sd',   sample = PK_2cmp_pop_draws$ClQ_sd,                              species = 'SD'),
                           data.frame(parameter = 'ka1_sd',   sample = PK_2cmp_pop_draws$ka1_sd,                              species = 'SD'),
                           data.frame(parameter = 'Vc_sd',    sample = PK_2cmp_pop_draws$Vc_sd,                               species = 'SD'),
                           data.frame(parameter = 'Vp_sd',    sample = PK_2cmp_pop_draws$Vp_sd,                               species = 'SD'),
                           data.frame(parameter = 'F_sd',     sample = PK_2cmp_pop_draws$F_sd,                                species = 'SD'),
                           PK_2cmp_pop_theta)
PK_2cmp_pop_theta$sample[!((PK_2cmp_pop_theta$parameter == 'F')|(PK_2cmp_pop_theta$parameter == 'Thalf')|(PK_2cmp_pop_theta$species == 'SD'))] <- exp(PK_2cmp_pop_theta$sample[!((PK_2cmp_pop_theta$parameter == 'F')|(PK_2cmp_pop_theta$parameter == 'Thalf')|(PK_2cmp_pop_theta$species == 'SD'))]) 
PK_2cmp_pop_theta$sample[(PK_2cmp_pop_theta$parameter == 'F')]  <- plogis(PK_2cmp_pop_theta$sample[(PK_2cmp_pop_theta$parameter == 'F')]) 
PK_2cmp_pop_theta$sample[((PK_2cmp_pop_theta$parameter == 'AUC')&(PK_2cmp_pop_theta$species == 'goat'))]  <- 1000000/(exp(PK_2cmp_pop_draws$Cl_mu))
PK_2cmp_pop_theta$sample[((PK_2cmp_pop_theta$parameter == 'AUC')&(PK_2cmp_pop_theta$species == 'sheep'))] <- 1000000/(exp(PK_2cmp_pop_draws$Cl_mu + PK_2cmp_pop_draws$species_Cl))
  
# Visualize the posterior distributions of the location parameters.
PK_2cmp_pop_plot_Cl    <- ggplot(data = PK_2cmp_pop_theta[PK_2cmp_pop_theta$parameter == 'Cl',],    aes(x = sample, fill = species))     + stat_slab(alpha = 0.5) + theme_bw() + ylab('') + theme(legend.title = element_blank()) + scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(0,15))   + scale_fill_viridis(discrete = TRUE, end = 0.8) + xlab('Plasma clearance (mL/min/kg)')
PK_2cmp_pop_plot_ClQ   <- ggplot(data = PK_2cmp_pop_theta[PK_2cmp_pop_theta$parameter == 'ClQ',],   aes(x = sample, fill = species))     + stat_slab(alpha = 0.5) + theme_bw() + ylab('') + theme(legend.title = element_blank()) + scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(0,15))   + scale_fill_viridis(discrete = TRUE, end = 0.8) + xlab('Intercompartmental clearance (mL/min/kg)')
PK_2cmp_pop_plot_ka1   <- ggplot(data = PK_2cmp_pop_theta[PK_2cmp_pop_theta$parameter == 'ka1',],   aes(x = sample, fill = species))     + stat_slab(alpha = 0.5) + theme_bw() + ylab('') + theme(legend.title = element_blank()) + scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(0,0.05)) + scale_fill_viridis(discrete = TRUE, end = 0.8) + xlab('Absorption rate constant (1/min)')
PK_2cmp_pop_plot_Vc    <- ggplot(data = PK_2cmp_pop_theta[PK_2cmp_pop_theta$parameter == 'Vc',],    aes(x = sample, fill = species))     + stat_slab(alpha = 0.5) + theme_bw() + ylab('') + theme(legend.title = element_blank()) + scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(0,3000)) + scale_fill_viridis(discrete = TRUE, end = 0.8) + xlab('Central volume of distribution (mL/kg)')
PK_2cmp_pop_plot_Vp    <- ggplot(data = PK_2cmp_pop_theta[PK_2cmp_pop_theta$parameter == 'Vp',],    aes(x = sample, fill = species))     + stat_slab(alpha = 0.5) + theme_bw() + ylab('') + theme(legend.title = element_blank()) + scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(0,3000)) + scale_fill_viridis(discrete = TRUE, end = 0.8) + xlab('Peripheral volume of distribution (mL/kg)')
PK_2cmp_pop_plot_F     <- ggplot(data = PK_2cmp_pop_theta[PK_2cmp_pop_theta$parameter == 'F',],     aes(x = sample*100, fill = species)) + stat_slab(alpha = 0.5) + theme_bw() + ylab('') + theme(legend.title = element_blank()) + scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(75,100)) + scale_fill_viridis(discrete = TRUE, end = 0.8) + xlab('Bioavailability (%)')
PK_2cmp_pop_plot_Thalf <- ggplot(data = PK_2cmp_pop_theta[PK_2cmp_pop_theta$parameter == 'Thalf',], aes(x = sample, fill = species))     + stat_slab(alpha = 0.5) + theme_bw() + ylab('') + theme(legend.title = element_blank()) + scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(0,1000)) + scale_fill_viridis(discrete = TRUE, end = 0.8) + xlab('Terminal half-life (min)')
PK_2cmp_pop_plot_AUC   <- ggplot(data = PK_2cmp_pop_theta[PK_2cmp_pop_theta$parameter == 'AUC',],   aes(x = sample, fill = species))     + stat_slab(alpha = 0.5) + theme_bw() + ylab('') + theme(legend.title = element_blank()) + scale_y_continuous(breaks = NULL) + coord_cartesian(xlim = c(0,500000)) + scale_fill_viridis(discrete = TRUE, end = 0.8) + xlab('Intravenous AUC (ng\U00B7min/mL)') + scale_x_continuous(breaks = c(0,100000,200000,300000,400000,500000), labels = c('0','100000','200000','300000','400000','500000'))
PK_2cmp_pop_plots <- ggarrange(PK_2cmp_pop_plot_Cl,PK_2cmp_pop_plot_ClQ,PK_2cmp_pop_plot_ka1,PK_2cmp_pop_plot_F,PK_2cmp_pop_plot_Vc,PK_2cmp_pop_plot_Vp,PK_2cmp_pop_plot_Thalf,PK_2cmp_pop_plot_AUC, nrow = 4, ncol = 2, common.legend = TRUE, legend = 'top')
ggsave(PK_2cmp_pop_plots, file = 'PK_2cmp_pop_plots.svg', units = 'mm', height = 220, width = 180)

# Collate the posterior estimates of the location parameters.
PK_2cmp_theta_table <- data.frame(parameter = character(length = 22), species = character(length = 22), q5 = numeric(length = 22), q50 = numeric(length = 22), q95 = numeric(length = 22))
PK_2cmp_theta_table[,c(1,2)] <- (unique(PK_2cmp_pop_theta[c('parameter', 'species')]))
for(i in 1:(dim(PK_2cmp_theta_table)[1])){
  PK_2cmp_theta_table[i,3:5] <- round(quantile(PK_2cmp_pop_theta$sample[((PK_2cmp_pop_theta$parameter == PK_2cmp_theta_table$parameter[i])&(PK_2cmp_pop_theta$species == PK_2cmp_theta_table$species[i]))], c(0.05,0.50,0.95)),3)
}
write.csv(PK_2cmp_theta_table, file = 'PK_2cmp_theta_table.csv', row.names = FALSE)

# Extract the population parameter estimates (all samples) and generate prediction graphics.
PK_goat_IV_mat  <- matrix(NA, nrow = length(seq(1,7200,1)), ncol = length(PK_2cmp_pop_draws$Cl_mu))
PK_sheep_IV_mat <- matrix(NA, nrow = length(seq(1,7200,1)), ncol = length(PK_2cmp_pop_draws$Cl_mu))
PK_goat_PO_mat  <- matrix(NA, nrow = length(seq(1,7200,1)), ncol = length(PK_2cmp_pop_draws$Cl_mu))
PK_sheep_PO_mat <- matrix(NA, nrow = length(seq(1,7200,1)), ncol = length(PK_2cmp_pop_draws$Cl_mu))
for(i in 1:length(PK_2cmp_pop_draws$Cl_mu)){
  PK_goat_IV_mat[,i]  <-  ((pk_2cmp_mod_IV(time = seq(1,7200,1), 1000000, exp(PK_2cmp_pop_draws$Cl_mu[i]), exp(PK_2cmp_pop_draws$ClQ_mu[i]), exp(PK_2cmp_pop_draws$Vc_mu[i]), exp(PK_2cmp_pop_draws$Vp_mu[i])))[[1]])[,1]
  PK_sheep_IV_mat[,i] <-  ((pk_2cmp_mod_IV(time = seq(1,7200,1), 1000000, exp(PK_2cmp_pop_draws$Cl_mu[i] + PK_2cmp_pop_draws$species_Cl[i]), exp(PK_2cmp_pop_draws$ClQ_mu[i] + PK_2cmp_pop_draws$species_ClQ[i]), exp(PK_2cmp_pop_draws$Vc_mu[i] + PK_2cmp_pop_draws$species_Vc[i]), exp(PK_2cmp_pop_draws$Vp_mu[i] + PK_2cmp_pop_draws$species_Vp[i])))[[1]])[,1]
  PK_goat_PO_mat[,i]  <-  ((pk_2cmp_mod_PO(time = seq(1,7200,1), 1000000, exp(PK_2cmp_pop_draws$Cl_mu[i]), exp(PK_2cmp_pop_draws$ClQ_mu[i]), exp(PK_2cmp_pop_draws$Vc_mu[i]), exp(PK_2cmp_pop_draws$Vp_mu[i]), exp(PK_2cmp_pop_draws$ka1_mu[i]), plogis(PK_2cmp_pop_draws$F_mu[i])))[[1]])[,1]
  PK_sheep_PO_mat[,i] <-  ((pk_2cmp_mod_PO(time = seq(1,7200,1), 1000000, exp(PK_2cmp_pop_draws$Cl_mu[i] + PK_2cmp_pop_draws$species_Cl[i]), exp(PK_2cmp_pop_draws$ClQ_mu[i] + PK_2cmp_pop_draws$species_ClQ[i]), exp(PK_2cmp_pop_draws$Vc_mu[i] + PK_2cmp_pop_draws$species_Vc[i]), exp(PK_2cmp_pop_draws$Vp_mu[i] + PK_2cmp_pop_draws$species_Vp[i]), exp(PK_2cmp_pop_draws$ka1_mu[i] + PK_2cmp_pop_draws$species_ka1[i]), plogis(PK_2cmp_pop_draws$F_mu[i] + PK_2cmp_pop_draws$species_F[i])))[[1]])[,1]
}
PK_goat_IV_pop   <- cbind(data.frame(time = seq(1,7200,1), dose = 1000000, species = 'goat',  route = 'IV'), rowQuantiles(PK_goat_IV_mat,  probs = c(0.05,0.25,0.50,0.75,0.95)))
PK_sheep_IV_pop  <- cbind(data.frame(time = seq(1,7200,1), dose = 1000000, species = 'sheep', route = 'IV'), rowQuantiles(PK_sheep_IV_mat, probs = c(0.05,0.25,0.50,0.75,0.95)))
PK_goat_PO_pop   <- cbind(data.frame(time = seq(1,7200,1), dose = 1000000, species = 'goat',  route = 'PO'), rowQuantiles(PK_goat_PO_mat,  probs = c(0.05,0.25,0.50,0.75,0.95)))
PK_sheep_PO_pop  <- cbind(data.frame(time = seq(1,7200,1), dose = 1000000, species = 'sheep', route = 'PO'), rowQuantiles(PK_sheep_PO_mat, probs = c(0.05,0.25,0.50,0.75,0.95)))
rm(PK_goat_IV_mat, PK_sheep_IV_mat, PK_goat_PO_mat, PK_sheep_PO_mat)
PK_pop_pred_data <- rbind(PK_goat_IV_pop, PK_sheep_IV_pop, PK_goat_PO_pop, PK_sheep_PO_pop)
PK_pop_pred_plots <- ggplot(data = PK_pop_pred_data, aes(x = time, y = `50%`)) + facet_grid(rows = vars(species), cols = vars(route)) + geom_line(aes(color = species)) + geom_ribbon(aes(ymin = `25%`, ymax = `75%`, fill = species), alpha = 0.3) + geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = species), alpha = 0.2) + theme_bw() + scale_y_log10(breaks = c(0.1,1,10,100,1000), labels = c('0.1','1','10','100','1000')) + coord_cartesian(xlim = c(0,5000), ylim = c(0.1,1000)) + geom_point(data = subject_PK_data_nozero, aes(x = time, y = Ct, shape = as.factor(CENS), color = species), fill = 'white', inherit.aes = FALSE, alpha = 0.6) + scale_shape_manual(values = c(21,19)) + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_fill_viridis(discrete = TRUE, end = 0.8) + theme(legend.position = 'none') + xlab('Time since adminstration (min)') + ylab('Drug concentration (ng/mL)')
ggsave(PK_pop_pred_plots, file = 'PK_pop_pred_plots.svg', units = 'mm', height = 240, width = 200)

# Extract the subject-level parameter estimates (all samples).
PK_2cmp_sub2_draws  <- extract(PK_2cmp_test_mod2, pars = c('Cl_VEC','ClQ_VEC','ka1_VEC','Vc_VEC','Vp_VEC','F_VEC'))
PK_2cmp_pred2_list <- vector(mode = 'list', length = dim(subject_PK_cov)[1])
for(i in 1:dim(subject_PK_cov)[1]){
  iter_data <- matrix(data = NA, nrow = 7200, ncol = 4000)
  if(subject_PK_cov$route[i] == 'IV'){
    for(h in 1:4000){
      iter_data[,h] <- (((pk_2cmp_mod_IV(time = seq(1,7200,1), 1000000, PK_2cmp_sub2_draws$Cl_VEC[h,i], PK_2cmp_sub2_draws$ClQ_VEC[h,i], PK_2cmp_sub2_draws$Vc_VEC[h,i], PK_2cmp_sub2_draws$Vp_VEC[h,i]))[[1]])[,1])
    }
  }else{
    for(h in 1:4000){
      iter_data[,h] <- (((pk_2cmp_mod_PO(seq(1,7200,1), 1000000, PK_2cmp_sub2_draws$Cl_VEC[h,i], PK_2cmp_sub2_draws$ClQ_VEC[h,i], PK_2cmp_sub2_draws$Vc_VEC[h,i], PK_2cmp_sub2_draws$Vp_VEC[h,i], PK_2cmp_sub2_draws$ka1_VEC[h,i], PK_2cmp_sub2_draws$F_VEC[h,i]))[[1]])[,1])
    }
  }
  PK_2cmp_pred2_list[[i]] <- cbind(data.frame(subject = subject_PK_cov$subject[i], time = seq(1,7200,1), route = subject_PK_cov$route[i], species = subject_PK_cov$species[i]),  as.data.frame(rowQuantiles(iter_data, probs = c(0.05,0.25,0.50,0.75,0.95))))
}
PK_2cmp_pred2_data <- list.rbind(PK_2cmp_pred2_list)

# Generate and export graphics of the subject-level concentration-time predictions and the observations.
PK_2cmp_pred2_plot_IV_log <- ggplot(data = PK_2cmp_pred2_data[PK_2cmp_pred2_data$route == 'IV',], aes(x = time, y = `50%`)) + geom_line(aes(color = species)) + geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = species), alpha = 0.2) + facet_wrap(~subject) + scale_y_log10(breaks = c(0.1,1,10,100,1000), labels = c('0.1','1','10','100','1000')) + scale_x_log10() + theme_bw() + coord_cartesian(xlim = c(1,10000), ylim = c(0.1,1000)) + geom_point(data = subject_PK_data_nozero[subject_PK_data_nozero$route == 'IV',], aes(x = time, y = Ct, shape = as.factor(CENS), color = species), fill = 'white', inherit.aes = FALSE, alpha = 0.6) + scale_shape_manual(values = c(21,19)) + guides(shape = 'none') + theme(legend.position = 'top', legend.title = element_blank()) + xlab('Time (min)') + ylab('Concentration (ng/mL)') + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_fill_viridis(discrete = TRUE, end = 0.8) + ggtitle('Intravenous administration')
ggsave(PK_2cmp_pred2_plot_IV_log, file = 'PK_2cmp_pred2_plot_IV_log.svg', units = 'mm', height = 200, width = 140)
PK_2cmp_pred2_plot_IV_lin <- ggplot(data = PK_2cmp_pred2_data[PK_2cmp_pred2_data$route == 'IV',], aes(x = time, y = `50%`)) + geom_line(aes(color = species)) + geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = species), alpha = 0.2) + facet_wrap(~subject) + scale_y_log10(breaks = c(0.1,1,10,100,1000), labels = c('0.1','1','10','100','1000')) +                   theme_bw() + coord_cartesian(xlim = c(1,5000), ylim = c(0.1,1000))  + geom_point(data = subject_PK_data_nozero[subject_PK_data_nozero$route == 'IV',], aes(x = time, y = Ct, shape = as.factor(CENS), color = species), fill = 'white', inherit.aes = FALSE, alpha = 0.6) + scale_shape_manual(values = c(21,19)) + guides(shape = 'none') + theme(legend.position = 'top', legend.title = element_blank()) + xlab('Time (min)') + ylab('Concentration (ng/mL)') + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_fill_viridis(discrete = TRUE, end = 0.8) + ggtitle('Intravenous administration')
ggsave(PK_2cmp_pred2_plot_IV_lin, file = 'PK_2cmp_pred2_plot_IV_lin.svg', units = 'mm', height = 200, width = 140)
PK_2cmp_pred2_plot_PO_log <- ggplot(data = PK_2cmp_pred2_data[PK_2cmp_pred2_data$route == 'PO',], aes(x = time, y = `50%`)) + geom_line(aes(color = species)) + geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = species), alpha = 0.2) + facet_wrap(~subject) + scale_y_log10(breaks = c(0.1,1,10,100,1000), labels = c('0.1','1','10','100','1000')) + scale_x_log10() + theme_bw() + coord_cartesian(xlim = c(1,10000), ylim = c(0.1,1000)) + geom_point(data = subject_PK_data_nozero[subject_PK_data_nozero$route == 'PO',], aes(x = time, y = Ct, shape = as.factor(CENS), color = species), fill = 'white', inherit.aes = FALSE, alpha = 0.6) + scale_shape_manual(values = c(21,19)) + guides(shape = 'none') + theme(legend.position = 'top', legend.title = element_blank()) + xlab('Time (min)') + ylab('Concentration (ng/mL)') + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_fill_viridis(discrete = TRUE, end = 0.8) + ggtitle('Oral administration')
ggsave(PK_2cmp_pred2_plot_PO_log, file = 'PK_2cmp_pred2_plot_PO_log.svg', units = 'mm', height = 200, width = 140)
PK_2cmp_pred2_plot_PO_lin <- ggplot(data = PK_2cmp_pred2_data[PK_2cmp_pred2_data$route == 'PO',], aes(x = time, y = `50%`)) + geom_line(aes(color = species)) + geom_ribbon(aes(ymin = `5%`, ymax = `95%`, fill = species), alpha = 0.2) + facet_wrap(~subject) + scale_y_log10(breaks = c(0.1,1,10,100,1000), labels = c('0.1','1','10','100','1000')) +                   theme_bw() + coord_cartesian(xlim = c(1,5000), ylim = c(0.1,1000))  + geom_point(data = subject_PK_data_nozero[subject_PK_data_nozero$route == 'PO',], aes(x = time, y = Ct, shape = as.factor(CENS), color = species), fill = 'white', inherit.aes = FALSE, alpha = 0.6) + scale_shape_manual(values = c(21,19)) + guides(shape = 'none') + theme(legend.position = 'top', legend.title = element_blank()) + xlab('Time (min)') + ylab('Concentration (ng/mL)') + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_fill_viridis(discrete = TRUE, end = 0.8) + ggtitle('Oral administration')
ggsave(PK_2cmp_pred2_plot_PO_lin, file = 'PK_2cmp_pred2_plot_PO_lin.svg', units = 'mm', height = 200, width = 140)
PK_2cmp_pred_lin_plots <- ggarrange(PK_2cmp_pred2_plot_IV_lin, PK_2cmp_pred2_plot_PO_lin, common.legend = TRUE, legend = 'bottom')
ggsave(PK_2cmp_pred_lin_plots, file = 'PK_2cmp_pred_lin_plots.svg', units = 'mm', height = 200, width = 375)

# Extract the posterior mean residuals and visualize with the data.
PK_2cmp_residuals <- data.frame(subject = subject_PK_data_nozero$subject, time = subject_PK_data_nozero$time, route = subject_PK_data_nozero$route, species = subject_PK_data_nozero$species, Ct = subject_PK_data_nozero$Ct, CENS = subject_PK_data_nozero$CENS, meanpred = colMeans((extract(PK_2cmp_test_mod2, pars = 'output_predictions'))$output_predictions))
PK_2cmp_residuals$residual <- PK_2cmp_residuals$meanpred-PK_2cmp_residuals$Ct
PK_2cmp_residuals$W_residual <- PK_2cmp_residuals$residual/PK_2cmp_residuals$meanpred
PK_2cmp_resplot  <- ggplot(data = PK_2cmp_residuals[PK_2cmp_residuals$CENS == 0,], aes(x = meanpred, y = W_residual, color = species)) + geom_point(alpha = 0.6) + facet_wrap(~route, ncol = 2) + theme_bw() + theme(legend.position = 'top', legend.title = element_blank()) + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_x_log10() + geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.6) + xlab('Model-predicted concentration (ng/mL)') + ylab('Weighted (1/Y) residual')
PK_2cmp_identity <- ggplot(data = PK_2cmp_residuals[PK_2cmp_residuals$CENS == 0,], aes(x = Ct, y = meanpred, color = species)) + geom_point(alpha = 0.6) + facet_wrap(~route, ncol = 2) + theme_bw() + theme(legend.position = 'top', legend.title = element_blank()) + scale_color_viridis(discrete = TRUE, end = 0.8) + scale_x_log10() + scale_y_log10() + geom_abline(slope = 1, intercept = 0, alpha = 0.6) + xlab('Observed concentration (ng/mL)') + ylab('Model-predicted concentration (ng/mL)')
PK_2cmp_checkplots <- ggarrange(PK_2cmp_resplot, PK_2cmp_identity, ncol = 1, common.legend = TRUE)
ggsave(PK_2cmp_checkplots, file = 'PK_2cmp_checkplots.svg', units = 'mm', height = 220, width = 180)

# Extract the individual parameter estimates.
PK_2cmp_draws_Cl_list  <- vector(mode = 'list', length = dim(subject_PK_cov[subject_PK_cov$route == 'IV',])[1])
PK_2cmp_draws_ClQ_list <- vector(mode = 'list', length = dim(subject_PK_cov[subject_PK_cov$route == 'IV',])[1])
PK_2cmp_draws_ka1_list <- vector(mode = 'list', length = dim(subject_PK_cov[subject_PK_cov$route == 'IV',])[1])
PK_2cmp_draws_Vc_list  <- vector(mode = 'list', length = dim(subject_PK_cov[subject_PK_cov$route == 'IV',])[1])
PK_2cmp_draws_Vp_list  <- vector(mode = 'list', length = dim(subject_PK_cov[subject_PK_cov$route == 'IV',])[1])
PK_2cmp_draws_F_list   <- vector(mode = 'list', length = dim(subject_PK_cov[subject_PK_cov$route == 'IV',])[1])
for(i in 1:dim(subject_PK_cov[subject_PK_cov$route == 'IV',])[1]){
  PK_2cmp_draws_Cl_list[[i]]  <- data.frame(subject = subject_PK_cov$subject[i], route = subject_PK_cov$route[i], species = subject_PK_cov$species[i], sample = PK_2cmp_sub2_draws$Cl_VEC[,i])
  PK_2cmp_draws_ClQ_list[[i]] <- data.frame(subject = subject_PK_cov$subject[i], route = subject_PK_cov$route[i], species = subject_PK_cov$species[i], sample = PK_2cmp_sub2_draws$ClQ_VEC[,i])
  PK_2cmp_draws_ka1_list[[i]] <- data.frame(subject = subject_PK_cov$subject[i], route = subject_PK_cov$route[i], species = subject_PK_cov$species[i], sample = PK_2cmp_sub2_draws$ka1_VEC[,i])
  PK_2cmp_draws_Vc_list[[i]]  <- data.frame(subject = subject_PK_cov$subject[i], route = subject_PK_cov$route[i], species = subject_PK_cov$species[i], sample = PK_2cmp_sub2_draws$Vc_VEC[,i])
  PK_2cmp_draws_Vp_list[[i]]  <- data.frame(subject = subject_PK_cov$subject[i], route = subject_PK_cov$route[i], species = subject_PK_cov$species[i], sample = PK_2cmp_sub2_draws$Vp_VEC[,i])
  PK_2cmp_draws_F_list[[i]]   <- data.frame(subject = subject_PK_cov$subject[i], route = subject_PK_cov$route[i], species = subject_PK_cov$species[i], sample = PK_2cmp_sub2_draws$F_VEC[,i])
}
PK_2cmp_draws_Cl  <- list.rbind(PK_2cmp_draws_Cl_list)
PK_2cmp_draws_ClQ <- list.rbind(PK_2cmp_draws_ClQ_list)
PK_2cmp_draws_ka1 <- list.rbind(PK_2cmp_draws_ka1_list)
PK_2cmp_draws_Vc  <- list.rbind(PK_2cmp_draws_Vc_list)
PK_2cmp_draws_Vp  <- list.rbind(PK_2cmp_draws_Vp_list)
PK_2cmp_draws_F   <- list.rbind(PK_2cmp_draws_F_list)
rm(PK_2cmp_draws_Cl_list,PK_2cmp_draws_ClQ_list,PK_2cmp_draws_ka1_list,PK_2cmp_draws_Vc_list,PK_2cmp_draws_Vp_list,PK_2cmp_draws_F_list)

# Generate interval plots of the individual parameter estimates.
sub_est_plot_Cl  <- ggplot(data = PK_2cmp_draws_Cl,  aes(x = sample, y = subject, color = species))     + stat_pointinterval(alpha = 0.8, point_interval = median_qi, .width = c(0.5,0.9)) + geom_point(data = pk_ind_theta, aes(y = subject, x = clearance, color = species), inherit.aes = FALSE, fill = 'white', shape = 21) + theme_bw() + theme(legend.title = element_blank()) + scale_y_discrete(limits = rev) + coord_cartesian(xlim = c(0,20)) + scale_color_viridis(discrete = TRUE, end = 0.8) + ylab('') + xlab('Plasma clearance (L/min/kg)')
sub_est_plot_ClQ <- ggplot(data = PK_2cmp_draws_ClQ, aes(x = sample, y = subject, color = species))     + stat_pointinterval(alpha = 0.8, point_interval = median_qi, .width = c(0.5,0.9)) + geom_point(data = pk_ind_theta, aes(y = subject, x = Q,         color = species), inherit.aes = FALSE, fill = 'white', shape = 21) + theme_bw() + theme(legend.title = element_blank()) + scale_y_discrete(limits = rev) + coord_cartesian(xlim = c(0,20)) + scale_color_viridis(discrete = TRUE, end = 0.8) + ylab('') + xlab('Intercompartmental clearance (L/min/kg)')
sub_est_plot_ka1 <- ggplot(data = PK_2cmp_draws_ka1, aes(x = sample, y = subject, color = species))     + stat_pointinterval(alpha = 0.8, point_interval = median_qi, .width = c(0.5,0.9)) + geom_point(data = pk_ind_theta, aes(y = subject, x = ka,        color = species), inherit.aes = FALSE, fill = 'white', shape = 21) + theme_bw() + theme(legend.title = element_blank()) + scale_y_discrete(limits = rev) + coord_cartesian(xlim = c(0,0.05)) + scale_color_viridis(discrete = TRUE, end = 0.8) + ylab('') + xlab('Absorption rate constant (1/min)')
sub_est_plot_Vc  <- ggplot(data = PK_2cmp_draws_Vc,  aes(x = sample, y = subject, color = species))     + stat_pointinterval(alpha = 0.8, point_interval = median_qi, .width = c(0.5,0.9)) + geom_point(data = pk_ind_theta, aes(y = subject, x = VC,        color = species), inherit.aes = FALSE, fill = 'white', shape = 21) + theme_bw() + theme(legend.title = element_blank()) + scale_y_discrete(limits = rev) + coord_cartesian(xlim = c(0,3000))    + scale_color_viridis(discrete = TRUE, end = 0.8) + ylab('') + xlab('Central volume of distribution (L/kg)')
sub_est_plot_Vp  <- ggplot(data = PK_2cmp_draws_Vp,  aes(x = sample, y = subject, color = species))     + stat_pointinterval(alpha = 0.8, point_interval = median_qi, .width = c(0.5,0.9)) + geom_point(data = pk_ind_theta, aes(y = subject, x = VP,        color = species), inherit.aes = FALSE, fill = 'white', shape = 21) + theme_bw() + theme(legend.title = element_blank()) + scale_y_discrete(limits = rev) + coord_cartesian(xlim = c(0,3000))    + scale_color_viridis(discrete = TRUE, end = 0.8) + ylab('') + xlab('Peripheral volume of distribution (L/kg)')
sub_est_plot_F   <- ggplot(data = PK_2cmp_draws_F,   aes(x = 100*sample, y = subject, color = species)) + stat_pointinterval(alpha = 0.8, point_interval = median_qi, .width = c(0.5,0.9)) + geom_point(data = pk_ind_theta, aes(y = subject, x = 100*`F`,   color = species), inherit.aes = FALSE, fill = 'white', shape = 21) + theme_bw() + theme(legend.title = element_blank()) + scale_y_discrete(limits = rev) + coord_cartesian(xlim = c(75,100)) + scale_color_viridis(discrete = TRUE, end = 0.8) + ylab('') + xlab('Bioavailability (%)')
sub_est_plots <- ggarrange(sub_est_plot_Cl, sub_est_plot_ClQ, sub_est_plot_ka1, sub_est_plot_F, sub_est_plot_Vc, sub_est_plot_Vp, nrow = 3, ncol = 2, common.legend = TRUE, legend = 'bottom')
ggsave(sub_est_plots, file = 'sub_est_plots.svg', units = 'mm', width = 180, height = 200)

# Specify a function that generates tolerance interval statements regarding AUC compared between sheep and goats.
#     Statements this complex may be overkill for many applications, but here serves as an example of the generation of uncertainty statements from arbitrary functions of parameters.
tol_int_PK <- function(sheep_dose,simdose,target_quantile){
  nsamples <- length(PK_2cmp_pop_draws$Cl_mu)
  goat_AUC_quant <- matrix(NA, nrow = length(simdose),ncol = nsamples)
  goat_AUC_dif   <- matrix(NA, nrow = length(simdose),ncol = nsamples)
  sheep_AUC      <- numeric(length(nsamples))
  dose_optim     <- numeric(length(nsamples))
  for(i in 1:nsamples){
    sim_sheep_theta <- data.frame(Cl = numeric(length = 1), F = numeric(length = 1))
    sim_sheep_theta$Cl <- exp(PK_2cmp_pop_draws$Cl_mu[i]+PK_2cmp_pop_draws$species_Cl[i])
    sim_sheep_theta$F  <- plogis(PK_2cmp_pop_draws$F_mu[i]+PK_2cmp_pop_draws$species_F[i])
    #
    sim_goat_theta  <- data.frame(Cl = numeric(length = 1000), F = numeric(length = 1000))
    sim_goat_theta$Cl <- exp(rnorm(1000,PK_2cmp_pop_draws$Cl_mu[i],PK_2cmp_pop_draws$Cl_sd[i]))
    sim_goat_theta$F  <- plogis(rnorm(1000,PK_2cmp_pop_draws$F_mu[i],PK_2cmp_pop_draws$F_sd[i]))
    sim_goat_AUC <- data.frame(dose = numeric(length(simdose)), quant = numeric(length(simdose)))
      for(h in 1:length(simdose)){
        pred_AUC <- (simdose[h]*sim_goat_theta$F)/sim_goat_theta$Cl
        sim_goat_AUC$dose[h] <- simdose[h]
        sim_goat_AUC$quant[h]  <- as.numeric(quantile(pred_AUC,target_quantile))
        sim_goat_AUC$sheep <- (sheep_dose*sim_sheep_theta$F)/sim_sheep_theta$Cl
      }
    sim_goat_AUC$diff <- sim_goat_AUC$quant-sim_goat_AUC$sheep
    goat_AUC_quant[,i] <- sim_goat_AUC$quant
    goat_AUC_dif[,i]   <- sim_goat_AUC$diff
    sheep_AUC[i] <- (sheep_dose*sim_sheep_theta$F)/sim_sheep_theta$Cl
    AUC_dose_mod <- smooth.spline(sim_goat_AUC$dose, sim_goat_AUC$diff)
    AUC_dose_smod_fun <- function(x){abs(predict(AUC_dose_mod, x)$y)}
    dose_optim[i] <- (optimize(AUC_dose_smod_fun, c(0,10000000)))$minimum
  }
  list(simdose,goat_AUC_quant, goat_AUC_dif, sheep_AUC, dose_optim)
}

# Obtain tolerance intervals for the AUC after oral dosing in goats relative to a reference oral dose in sheep, and plot the results.
tol_int_samples <- tol_int_PK(1000000,seq(0,10000000,10000),0.2)
goat_AUC_est <- cbind(data.frame(dose = tol_int_samples[[1]]), rowQuantiles(tol_int_samples[[2]], probs = c(0.05,0.25,0.5,0.75,0.95)))
diff_AUC_est <- cbind(data.frame(dose = tol_int_samples[[1]]), rowQuantiles(tol_int_samples[[3]], probs = c(0.05,0.25,0.5,0.75,0.95)))
sheep_AUC_est <- (quantile(tol_int_samples[[4]], probs = c(0.1,0.5,0.9)))
goat_AUC_perc_plot <- ggplot(goat_AUC_est, aes(x = dose, y = `50%`)) + geom_line() + geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.3) + scale_x_continuous(breaks = c(0,2500000,5000000,7500000,10000000), labels = c(0,2500000,5000000,7500000,10000000)/10^6) + theme_bw() + scale_y_continuous(breaks = c(0,500000,1000000), labels = c('0','500000','1000000')) + geom_hline(yintercept = c(sheep_AUC_est[1],sheep_AUC_est[3]), linetype = 'dashed', alpha = 0.6) + geom_hline(yintercept = c(sheep_AUC_est[2]), alpha = 0.7) + ylab('AUC (ng\U00B7min/mL)') + xlab('Dose (mg/kg)')
goat_AUC_diff_plot <- ggplot(diff_AUC_est, aes(x = dose, y = `50%`)) + geom_line() + geom_ribbon(aes(ymin = `5%`, ymax = `95%`), alpha = 0.2) + geom_ribbon(aes(ymin = `25%`, ymax = `75%`), alpha = 0.3) + scale_x_continuous(breaks = c(0,2500000,5000000,7500000,10000000), labels = c(0,2500000,5000000,7500000,10000000)/10^6) + theme_bw() + scale_y_continuous(breaks = c(-250000,0,250000,500000,750000), labels = c('-250000','0','250000','500000','750000')) + geom_hline(yintercept = 0, linetype = 'dotted', alpha = 1) + ylab('Difference in AUC') + ylab('Difference in AUC (ng\U00B7min/mL)') + xlab('Dose (mg/kg)')
goat_AUC_pred_plots <- ggarrange(goat_AUC_perc_plot, goat_AUC_diff_plot, ncol = 1)
ggsave(goat_AUC_pred_plots, file = 'goat_AUC_pred_plots.svg', units = 'mm', height = 240, width = 200)
quantile(tol_int_samples[[5]], c(0.05,0.50,0.95))
 
# Export the simulated dataset.
write.csv(subject_PK_data, file = 'simulated_pharmacokinetic_data_example.csv', row.names = FALSE)
