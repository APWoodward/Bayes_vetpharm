# AP. Woodward, University of Canberra, 2023.
# Generate a simulated dose-response experiment and analysis of the data with a nonlinear model.
# The objective of this implementation is to illustrate principles of Bayesian statistics relevant to veterinary pharmacology, especially the specification of priors and post-processing of the model to generate expressive graphics.
#     This R specification could be easily adapted to studies of common designs including other examples of dose-response system mdoeling.
# The analysis is focused on a Bayesian hierarchical generalized nonlinear model with weakly informative priors and is implemented with 'brms' (https://doi.org/10.18637/jss.v080.i01).

# Load the required packages.
library(rlist)
library(brms)
library(brmsmargins)
library(ggplot2)
library(reshape2)
library(ggdist)
library(ggpubr)
library(distributional)
library(viridis)
library(bayesplot)
library(ggh4x)

# This is a hierarchical nonlinear regression model for a binomial response, which will be set up with the beta-binomial response distribution.
# The design is loosely based on that of (https://doi.org/10.1111/jvim.16587).
#       There is between-subject variation in all the PD parameters, and a 'tumor type' covariate which influences only the EC50 ('location' shift).
#       All the subject-level distributions for these parameters are independent (there are no correlated PD parameters).
#       The prior distributions for the subject-level variability are log-normal; later, hyperpriors will be specified.

# Specify the study design.
n_subjects = 30
subject_PD_theta <- data.frame(upper_asym = numeric(length = n_subjects), lower_asym = numeric(length = n_subjects), EC50 = numeric(length = n_subjects), hill_coef = numeric(length = n_subjects))
subject_PD_theta$tumor_type <- c(rep('A',15),rep('B',15))

# Define the upper and lower asymptote; in this case specify that there are no covariate effects.
#       This relaxes the typical assumption that the upper asymptote is '100%' (or otherwise known in advance).
#       Note that in the model the upper asymptote is not directly estimated, but as upper_asym-lower_asym, to enforce ordering of the parameters.
#       The difference variable (to be estimated) is therefore normally distributed with mean 4 and SD sqrt(0.5), so on the linear predictor scale its SD is ~0.17.
subject_PD_theta$upper_asym <- plogis(rnorm(n_subjects, mean = 2, sd = 0.5))
subject_PD_theta$lower_asym <- plogis(rnorm(n_subjects, mean = -2, sd = 0.5))

# Define the EC50, which is based on a categorical covariate for tumor type.
subject_PD_theta$EC50[subject_PD_theta$tumor_type == 'A'] <- (exp(rnorm(n_subjects/2, mean = log(100), sd = 1)))
subject_PD_theta$EC50[subject_PD_theta$tumor_type == 'B'] <- (exp(rnorm(n_subjects/2, mean = log(500), sd = 1)))

# Define the Hill coefficient (steepness); in this case specify that there are no covariate effects.
subject_PD_theta$hill_coef <- exp(rnorm(n_subjects, mean = log(1), sd = 0.2))

# Specify a function containing the structural model (inhibitory 'sigmoid Emax', i.e. Hill model)
four_parameter_hill <- function(X, upper, lower, EC50, hill){
  resp <- upper - (((upper-lower)*(X^hill))/((EC50^hill) + (X^hill)))
  return(resp)
}

# Now simulate the concentration-response data; instead of zero which cannot appear on log-scale, we will set 'zero' (untreated control) to an arbitrarily-small concentration.
#       Beta-distributed noise is added, and then the binomial counts are drawn as the manifest observations.
#       A random series of observations (10% of the total) are removed (to emulate some missingness).
sim_PD_list <- vector(mode = 'list', length = n_subjects)
for (i in 1:n_subjects){
  sub_conc <- c(0,0.1,1,10,50,100,200,500,1000,2000,5000,10000,50000)
  sub_response = four_parameter_hill(X = sub_conc, subject_PD_theta$upper_asym[i], subject_PD_theta$lower_asym[i], subject_PD_theta$EC50[i], subject_PD_theta$hill_coef[i])
  sub_data <- data.frame(conc = sub_conc, response = sub_response, subject = i, tumor_type = subject_PD_theta$tumor_type[i])
  sim_PD_list[[i]] <- sub_data 
}
sim_PD_data <- list.rbind(sim_PD_list)
sim_PD_data$total <- round(runif(length(sim_PD_data$conc), 5000, 10000))
sim_PD_data$viable <- rbeta_binomial(length(sim_PD_data$total),sim_PD_data$total,sim_PD_data$response,100)
sim_PD_data$prop <- sim_PD_data$viable/sim_PD_data$total
sim_PD_data <- sim_PD_data[(-1*(sample((1:dim(sim_PD_data)[1]), dim(sim_PD_data)[1]*0.1))),]
ggplot(data = sim_PD_data, aes(x = conc, y = prop)) + geom_point() + facet_wrap(~subject) + scale_x_log10() + theme_bw()

# Specify and estimate the nonlinear hierarchical model.
count_PD_form  <- bf(viable|trials(total) ~ inv_logit(lower+exp(upperdif)) - (((inv_logit(lower+exp(upperdif))-inv_logit(lower))*(conc^(exp(hill))))/((exp(EC50)^(exp(hill))) + (conc^(exp(hill))))), upperdif + lower + EC50 + hill ~ 1 + tumor_type + (1|a|subject), nl = TRUE)
count_PD_priors <- set_prior('normal(0,1)', class = 'b', coef = 'Intercept', nlpar = 'hill') + set_prior('normal(0,0.25)', class = 'b', coef = 'tumor_typeB', nlpar = 'hill') + set_prior('normal(0,1)', class = 'sd', nlpar = 'hill') + set_prior('normal(1,1)', class = 'b', coef = 'Intercept', nlpar = 'upperdif') + set_prior('normal(0,0.5)', class = 'b', coef = 'tumor_typeB', nlpar = 'upperdif') + set_prior('normal(0,1)', class = 'sd', nlpar = 'upperdif') + set_prior('normal(0,3)', class = 'b', coef = 'Intercept', nlpar = 'lower') + set_prior('normal(0,1)', class = 'b', coef = 'tumor_typeB', nlpar = 'lower') + set_prior('normal(0,1)', class = 'sd', nlpar = 'lower') + set_prior('normal(0,5)', class = 'b', coef = 'Intercept', nlpar = 'EC50') + set_prior('normal(0,5)', class = 'b', coef = 'tumor_typeB', nlpar = 'EC50') + set_prior('normal(0,1)', class = 'sd', nlpar = 'EC50')
count_PD_inits <- function(){
  list(b_upperdif = runif(1,-2,2), 
       b_lower = runif(1,-2,2), 
       b_EC50  = runif(1,-5,5),
       b_hill = runif(1,-2,2))
}
count_PD_model <- brm(count_PD_form, data = sim_PD_data, family = beta_binomial(link = 'identity'), init = count_PD_inits, chains = 4, cores = 4, iter = 5000, control = list(adapt_delta = 0.9, max_treedepth = 12), prior = count_PD_priors)
summary(count_PD_model)
pp_check(count_PD_model, ndraws = 30)
count_PD_fixef <- fixef(count_PD_model, probs = c(0.05,0.95), robust = TRUE)
count_PD_fixef <- cbind(data.frame(parameter = rownames(count_PD_fixef)), count_PD_fixef)
rownames(count_PD_fixef) <- NULL
count_PD_varcorr <- (VarCorr(count_PD_model, probs = c(0.05,0.95), robust = TRUE))$subject$sd
count_PD_varcorr <- cbind(data.frame(parameter = rownames(count_PD_varcorr)), count_PD_varcorr)
rownames(count_PD_varcorr) <- NULL
write.csv(count_PD_fixef, file = 'count_PD_fixef.csv', row.names = FALSE)
write.csv(count_PD_varcorr, file = 'count_PD_varcorr.csv', row.names = FALSE)

# Define a function for the input concentration X; given values of the other parameters X corresponding to the specified response is obtained by rearrangement.
#   Parameters; theta[1]: Emax, theta[2]: EC50, theta[3]: E0, theta[4]: H.
EC_find_fun <- function(theta, quant){
  Y_target <- ((theta[1] - theta[3])*(1-quant)) + theta[3]
  X_target <- ((((theta[1]*(theta[2]^theta[4]))/(Y_target-theta[3])) - ((Y_target*(theta[2]^theta[4]))/(Y_target-theta[3]))))^(1/theta[4])  
  out_list <- list(X_target, Y_target)
  return(out_list)
}

# Obtain population estimates of the EC90 by tumor type.
pop_PD_A   <- data.frame(tumor_type = 'A', total = 1, conc = 1)
pop_thetaA <- matrix(ncol = 4, nrow = ndraws(count_PD_model))
pop_PD_B   <- data.frame(tumor_type = 'B', total = 1, conc = 1)
pop_thetaB <- matrix(ncol = 4, nrow = ndraws(count_PD_model))
pop_thetaA[,1] <- plogis(fitted(count_PD_model, newdata = pop_PD_A, nlpar = 'lower', summary = FALSE, re_formula = NA) + exp(fitted(count_PD_model, newdata = pop_PD_A, nlpar = 'upperdif', summary = FALSE, re_formula = NA)))
pop_thetaA[,2] <- exp(fitted(count_PD_model, newdata = pop_PD_A,    nlpar = 'EC50',  summary = FALSE, re_formula = NA))
pop_thetaA[,3] <- plogis(fitted(count_PD_model, newdata = pop_PD_A, nlpar = 'lower', summary = FALSE, re_formula = NA))
pop_thetaA[,4] <- exp(fitted(count_PD_model, newdata = pop_PD_A,    nlpar = 'hill',  summary = FALSE, re_formula = NA))
pop_thetaB[,1] <- plogis(fitted(count_PD_model, newdata = pop_PD_B, nlpar = 'lower', summary = FALSE, re_formula = NA) + exp(fitted(count_PD_model, newdata = pop_PD_B, nlpar = 'upperdif', summary = FALSE, re_formula = NA)))
pop_thetaB[,2] <- exp(fitted(count_PD_model, newdata = pop_PD_B,    nlpar = 'EC50',  summary = FALSE, re_formula = NA))
pop_thetaB[,3] <- plogis(fitted(count_PD_model, newdata = pop_PD_B, nlpar = 'lower', summary = FALSE, re_formula = NA))
pop_thetaB[,4] <- exp(fitted(count_PD_model, newdata = pop_PD_B,    nlpar = 'hill',  summary = FALSE, re_formula = NA))
pop_PD_est <- data.frame(EC90_A = numeric(length = ndraws(count_PD_model)), y90_A = numeric(length = ndraws(count_PD_model)), EC90_B = numeric(length = ndraws(count_PD_model)), y90_B = numeric(length = ndraws(count_PD_model)))
for (i in 1:ndraws(count_PD_model)){
  pop_PD_est$EC90_A[i] <- (EC_find_fun(pop_thetaA[i,], 0.9))[[1]]
  pop_PD_est$y90_A[i]  <- (EC_find_fun(pop_thetaA[i,], 0.9))[[2]]
  pop_PD_est$EC90_B[i] <- (EC_find_fun(pop_thetaB[i,], 0.9))[[1]]
  pop_PD_est$y90_B[i]  <- (EC_find_fun(pop_thetaB[i,], 0.9))[[2]]
} 
signif(quantile(pop_PD_est$EC90_A, c(0.05,0.50,0.95)),4)
signif(quantile(pop_PD_est$y90_A, c(0.05,0.50,0.95)), 4)
signif(quantile(pop_PD_est$EC90_B, c(0.05,0.50,0.95)),4)
signif(quantile(pop_PD_est$y90_B, c(0.05,0.50,0.95)), 4)

# Generate individual-level (conditional) estimates of the EC90.
sub_PD_newdata <- cbind(unique(sim_PD_data[c('subject','tumor_type')]), data.frame(total = 1, conc = 1))
sub_PD_newdata$EC90_q5  <- NA
sub_PD_newdata$EC90_q50 <- NA
sub_PD_newdata$EC90_q95 <- NA
sub_PD_newdata$y90_q5   <- NA
sub_PD_newdata$y90_q50  <- NA
sub_PD_newdata$y90_q95  <- NA
for (i in 1:length(sub_PD_newdata$subject)){
  sub_theta     <- matrix(ncol = 4, nrow = ndraws(count_PD_model))
  sub_theta[,1] <- plogis(fitted(count_PD_model, newdata = sub_PD_newdata[i,], nlpar = 'lower', summary = FALSE) + exp(fitted(count_PD_model, newdata = sub_PD_newdata[i,], nlpar = 'upperdif', summary = FALSE)))
  sub_theta[,2] <- exp(fitted(count_PD_model, newdata = sub_PD_newdata[i,], nlpar = 'EC50', summary = FALSE))
  sub_theta[,3] <- plogis(fitted(count_PD_model, newdata = sub_PD_newdata[i,], nlpar = 'lower', summary = FALSE))
  sub_theta[,4] <- exp(fitted(count_PD_model, newdata = sub_PD_newdata[i,], nlpar = 'hill', summary = FALSE))
  sub_EC90_est <- matrix(ncol = 2, nrow = ndraws(count_PD_model))
  for (j in 1:(ndraws(count_PD_model))){
    sub_EC90_est[j,1] <- (EC_find_fun(sub_theta[j,], 0.9))[[1]]
    sub_EC90_est[j,2] <- (EC_find_fun(sub_theta[j,], 0.9))[[2]]
  }
  sub_PD_newdata[i,5:7] <-  quantile(sub_EC90_est[,1], c(0.05,0.50,0.95))
  sub_PD_newdata[i,8:10] <- quantile(sub_EC90_est[,2], c(0.05,0.50,0.95))
}

# Generate individual-level predictions for visualization, including the individual-level EC90.
count_PD_newdata <- rbind(expand.grid(subject = unique(sim_PD_data$subject[sim_PD_data$tumor_type == 'A']), tumor_type = 'A', conc = 10^seq(-6,8,0.1)),expand.grid(subject = unique(sim_PD_data$subject[sim_PD_data$tumor_type == 'B']), tumor_type = 'B', conc = 10^seq(-6,8,0.1)))
count_PD_newdata$total <- 1
count_PD_pred <- cbind(count_PD_newdata, fitted(count_PD_model, newdata = count_PD_newdata, probs = c(0.05,0.95)))    
count_PD_sub_plot <- ggplot(data = sim_PD_data, aes(x = conc, y = prop, color = tumor_type)) + geom_point(alpha = 0.5) + facet_wrap(~subject) + scale_x_log10(breaks = c(1e-03,1e+01,1e+05), labels = c('0.001', '10', '100000')) + theme_bw() + theme(legend.position = 'top') + xlab('Drug concentration (arbitrary unit)') + ylab('Proportion of viable cells') + geom_line(data = count_PD_pred, aes(x = conc, y = Estimate, color = tumor_type), inherit.aes = FALSE) + geom_ribbon(data = count_PD_pred, aes(x = conc, ymin = Q5, ymax = Q95, fill = tumor_type), inherit.aes = FALSE, alpha = 0.2) + scale_fill_viridis(discrete = TRUE, end = 0.8, option = 'plasma') + scale_color_viridis(discrete = TRUE, end = 0.8, option = 'plasma') + labs(fill = 'Tumor type', color = 'Tumor type')
count_PD_sub_plot <- count_PD_sub_plot + geom_point(data = sub_PD_newdata, aes(x = EC90_q50, y = y90_q50, color = tumor_type), inherit.aes = FALSE,fill = 'white', shape = 21, size = 2)
ggsave(count_PD_sub_plot, file = 'count_PD_sub_plot.svg', units = 'mm', height = 300, width = 250)

# Generate a descriptive graphic of the observations and the apparent population relationships. 
count_PD_pop <- cbind(expand.grid(conc = 10^seq(-6,8,0.1), tumor_type = c('A','B'), total = 1), fitted(count_PD_model, newdata = expand.grid(conc = 10^seq(-6,8,0.1), tumor_type = c('A','B'), total = 1), re_formula = NA, probs = c(0.05,0.95)))
PD_fit_pop_plot   <- ggplot(data = count_PD_pop, aes(x = conc, y = Estimate)) + geom_line(aes(color = tumor_type)) + geom_ribbon(aes(ymin = Q5, ymax = Q95, fill = tumor_type), alpha = 0.2) + theme_bw() + theme(legend.position = 'top') + xlab('Drug concentration (arbitrary unit)') + ylab('Proportion of viable cells') + scale_x_log10(breaks = c(1e-03,1e+01,1e+05), labels = c('0.001', '10', '100000')) + geom_point(data = sim_PD_data, aes(x = conc, y = prop, color = tumor_type), alpha = 0.1) + scale_color_viridis(discrete = TRUE, end = 0.8, option = 'plasma') + scale_fill_viridis(discrete = TRUE, end = 0.8, option = 'plasma') + labs(fill = 'Tumor type', color = 'Tumor type')
ggsave(PD_fit_pop_plot, file = 'PD_fit_pop_plot.svg', units = 'mm', height = 120, width = 160)

# Generate graphical analyses of the residuals.
count_PD_residual <- cbind(sim_PD_data, residual = residuals(count_PD_model)[,1], fitted = fitted(count_PD_model)[,1])
PD_pop_residual   <- ggplot(data = count_PD_residual, aes(x = fitted, y = residual)) + geom_point(alpha = 0.3) + geom_hline(yintercept = 0, linetype = 'dashed') + ylim(c(-1500,1500)) + facet_wrap(~tumor_type) + theme_bw() + ylab('Residual (count of viable cells)') + xlab('Model-predicted value (count of viable cells)')
ggsave(PD_pop_residual, file = 'PD_pop_residual.svg', units = 'mm', height = 120, width = 200)
PD_sub_residual   <- ggplot(data = count_PD_residual, aes(x = fitted, y = residual, color = tumor_type)) + geom_point(alpha = 0.5) + geom_hline(yintercept = 0, linetype = 'dashed', alpha = 0.5) + ylim(c(-1500,1500)) + theme_bw() + theme(legend.position = 'top') + ylab('Residual (count of viable cells)') + xlab('Model-predicted value (count of viable cells)') + facet_wrap(~subject) + scale_color_viridis(discrete = TRUE, end = 0.8, option = 'plasma') + labs(fill = 'Tumor type', color = 'Tumor type')
ggsave(PD_sub_residual, file = 'PD_sub_residual.svg', units = 'mm', height = 300, width = 250)
PD_all_residuals <- ggarrange(PD_pop_residual, PD_sub_residual, ncol = 1, heights = c(1,2))
ggsave(PD_all_residuals, file = 'PD_all_residuals.svg', units = 'mm', height = 220, width = 180)

# Summarize the parameter estimates; the posterior distributions, and the apparent population distributions of the PD parameters.
count_PD_loc_draws <- as.data.frame(cbind(fitted(count_PD_model, nlpar = 'upperdif',  re_formula = NA, newdata = data.frame(tumor_type = c('A','B'), conc = 1, total = 1), summary = FALSE),
                                          fitted(count_PD_model, nlpar = 'lower',     re_formula = NA, newdata = data.frame(tumor_type = c('A','B'), conc = 1, total = 1), summary = FALSE),
                                          fitted(count_PD_model, nlpar = 'EC50',      re_formula = NA, newdata = data.frame(tumor_type = c('A','B'), conc = 1, total = 1), summary = FALSE),
                                          fitted(count_PD_model, nlpar = 'hill',      re_formula = NA, newdata = data.frame(tumor_type = c('A','B'), conc = 1, total = 1), summary = FALSE)))
colnames(count_PD_loc_draws) <- c('A_upperdif_location', 'B_upperdif_location', 'A_lower_location', 'B_lower_location', 'A_EC50_location', 'B_EC50_location', 'A_hill_location', 'B_hill_location')
count_PD_loc_draws <- melt(count_PD_loc_draws)
count_PD_loc_draws$tumor_type <- substr(count_PD_loc_draws$variable, 1, 1)
count_PD_loc_draws$parameter <- substr(count_PD_loc_draws$variable, 3, nchar(as.character(count_PD_loc_draws$variable)))
count_PD_loc_unique <- unique(count_PD_loc_draws$parameter)
count_PD_names <- c('Median E0 (difference from Emax, logit scale)', 'Median Emax (proportion of viable cells)', 'Median EC50 (arbitrary unit)', 'Median Hill (steepness) coefficient', 'SD log(logit(E0))', 'SD logit(Emax)', 'SD log(EC50)', 'SD log(Hill coefficient)', 'E0 (difference from Emax, logit scale)', 'Emax (proportion of viable cells)', 'EC50 (arbitrary unit)', 'Hill (steepness) coefficient')
count_PD_parms <- data.frame(parameter = c('A_upperdif_location', 'B_upperdif_location', 'A_lower_location', 'B_lower_location', 'A_EC50_location', 'B_EC50_location', 'A_hill_location', 'B_hill_location'), mean = as.data.frame(fixef(count_PD_model, robust = TRUE))[,1], sd = as.data.frame(VarCorr(count_PD_model, robust = TRUE))[,1])
count_PD_sd_draws <- as.data.frame((VarCorr(count_PD_model, summary = FALSE))$subject$sd)
count_PD_sd_draws <- melt(count_PD_sd_draws)
count_PD_sd_unique <- unique(count_PD_sd_draws$variable)
count_PD_parms <- data.frame(mean = (c(fitted(count_PD_model, nlpar = 'upperdif',  re_formula = NA, robust = TRUE, newdata = data.frame(tumor_type = c('A','B'), conc = 1, total = 1))[,1],
                                       fitted(count_PD_model, nlpar = 'lower',     re_formula = NA, robust = TRUE, newdata = data.frame(tumor_type = c('A','B'), conc = 1, total = 1))[,1],
                                       fitted(count_PD_model, nlpar = 'EC50',      re_formula = NA, robust = TRUE, newdata = data.frame(tumor_type = c('A','B'), conc = 1, total = 1))[,1],
                                       fitted(count_PD_model, nlpar = 'hill',      re_formula = NA, robust = TRUE, newdata = data.frame(tumor_type = c('A','B'), conc = 1, total = 1))[,1])))
count_PD_parms$sd <- rep(((VarCorr(count_PD_model, robust = TRUE))$subject$sd)[,1], each = 2)
rownames(count_PD_parms) <- c('A_upperdif_location', 'B_upperdif_location', 'A_lower_location', 'B_lower_location', 'A_EC50_location', 'B_EC50_location', 'A_hill_location', 'B_hill_location')
count_PD_parms$tumor_type <- substr(rownames(count_PD_parms), 1, 1)
count_PD_parms$parameter <- substr(rownames(count_PD_parms), 3, nchar(as.character(rownames(count_PD_parms))))
count_PD_subject <- cbind(sim_PD_data[sim_PD_data$conc == 1,], data.frame(hill = (fitted(count_PD_model, robust = TRUE, nlpar = 'hill', newdata = sim_PD_data[sim_PD_data$conc == 1,])[,1]), upperdif = (fitted(count_PD_model, robust = TRUE, nlpar = 'upperdif', newdata = sim_PD_data[sim_PD_data$conc == 1,])[,1])), lower = (fitted(count_PD_model, robust = TRUE, nlpar = 'lower', newdata = sim_PD_data[sim_PD_data$conc == 1,])[,1]), EC50 = (fitted(count_PD_model, robust = TRUE, nlpar = 'EC50', newdata = sim_PD_data[sim_PD_data$conc == 1,])[,1]))

# Generate distribution plots of the hyperposteriors, and the apparent population distributions at the posterior medians of the hyperparameters.
#     Care to correctly express the scale for 'upperdif', because the modelled parameter is on logarithmic scale, and it is then exponentialized in the model to enforce positivity (log-log transformation).
count_PD_plots <- vector(mode = 'list', length = 12)
count_PD_plots[[1]]  <-  ggplot(data = count_PD_loc_draws[count_PD_loc_draws$parameter == count_PD_loc_unique[1],], aes(x = exp(value), y = tumor_type))              + stat_halfeye(alpha = 0.5) + geom_vline(linetype = 'dashed', xintercept = 4, alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Location')              + xlab(count_PD_names[1])  + ylab('')
count_PD_plots[[2]]  <-  ggplot(data = count_PD_loc_draws[count_PD_loc_draws$parameter == count_PD_loc_unique[2],], aes(x = plogis(value), y = tumor_type))           + stat_halfeye(alpha = 0.5) + geom_vline(linetype = 'dashed', xintercept = exp(-2), alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Location')                   + xlab(count_PD_names[2])  + ylab('')
count_PD_plots[[3]]  <-  ggplot(data = count_PD_loc_draws[count_PD_loc_draws$parameter == count_PD_loc_unique[3],], aes(x = exp(value), y = tumor_type))              + stat_halfeye(alpha = 0.5) + geom_vline(linetype = 'dashed', xintercept = c(100,500), alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Location*')                   + xlab(count_PD_names[3])  + ylab('') + coord_cartesian(xlim = c(0,1200))
count_PD_plots[[4]]  <-  ggplot(data = count_PD_loc_draws[count_PD_loc_draws$parameter == count_PD_loc_unique[4],], aes(x = exp(value), y = tumor_type))              + stat_halfeye(alpha = 0.5) + geom_vline(linetype = 'dashed', xintercept = 1, alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Location')                   + xlab(count_PD_names[4])  + ylab('')
count_PD_plots[[5]]  <-  ggplot(data = count_PD_sd_draws[count_PD_sd_draws$variable ==    count_PD_sd_unique[1],],  aes(x = value))                                   + stat_halfeye(alpha = 0.5) + geom_vline(linetype = 'dashed', xintercept = sqrt(log(1+(0.5/(4^2)))), alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Scale')                      + xlab(count_PD_names[5])  + ylab('') + scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = seq(0,2,0.2)) + coord_cartesian(xlim = c(0,2))
count_PD_plots[[6]]  <-  ggplot(data = count_PD_sd_draws[count_PD_sd_draws$variable ==    count_PD_sd_unique[2],],  aes(x = value))                                   + stat_halfeye(alpha = 0.5) + geom_vline(linetype = 'dashed', xintercept = 0.5, alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Scale')                      + xlab(count_PD_names[6])  + ylab('') + scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = seq(0,2,0.2)) + coord_cartesian(xlim = c(0,2))
count_PD_plots[[7]]  <-  ggplot(data = count_PD_sd_draws[count_PD_sd_draws$variable ==    count_PD_sd_unique[3],],  aes(x = value))                                   + stat_halfeye(alpha = 0.5) + geom_vline(linetype = 'dashed', xintercept = 1, alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Scale')                      + xlab(count_PD_names[7])  + ylab('') + scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = seq(0,2,0.2)) + coord_cartesian(xlim = c(0,2))
count_PD_plots[[8]]  <-  ggplot(data = count_PD_sd_draws[count_PD_sd_draws$variable ==    count_PD_sd_unique[4],],  aes(x = value))                                   + stat_halfeye(alpha = 0.5) + geom_vline(linetype = 'dashed', xintercept = 0.2, alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Scale')                      + xlab(count_PD_names[8])  + ylab('') + scale_y_continuous(breaks = NULL) + scale_x_continuous(breaks = seq(0,2,0.2)) + coord_cartesian(xlim = c(0,2))
count_PD_plots[[9]]  <-  ggplot(data = count_PD_parms[count_PD_parms$parameter == count_PD_loc_unique[1],], aes(xdist = exp(dist_normal(mean, sd)), y = tumor_type))  + stat_slab(alpha = 0.6)    + geom_vline(linetype = 'dashed', xintercept = 4, alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Population distribution') + xlab(count_PD_names[9])  + ylab('')                                     + scale_x_continuous()                                                                                            + geom_point(data = count_PD_subject, aes(x = exp(upperdif), y = tumor_type), inherit.aes = FALSE, alpha = 1, shape = '|', size = 2)
count_PD_plots[[10]] <-  ggplot(data = count_PD_parms[count_PD_parms$parameter == count_PD_loc_unique[2],], aes(xdist = (dist_normal(mean, sd)), y = tumor_type))     + stat_slab(alpha = 0.6)    + geom_vline(linetype = 'dashed', xintercept = qlogis(exp(-2)), alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Population distribution')    + xlab(count_PD_names[10]) + ylab('')                                     + scale_x_continuous(breaks = qlogis(seq(0.05,0.30,0.05)), labels = function(x)(round(plogis(x), digits = 2)))    + geom_point(data = count_PD_subject, aes(x = (lower),       y = tumor_type), inherit.aes = FALSE, alpha = 1, shape = '|', size = 2)
count_PD_plots[[11]] <-  ggplot(data = count_PD_parms[count_PD_parms$parameter == count_PD_loc_unique[3],], aes(xdist = exp(dist_normal(mean, sd)), y = tumor_type))  + stat_slab(alpha = 0.6)    + geom_vline(linetype = 'dashed', xintercept = c(100,500), alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Population distribution*') + xlab(count_PD_names[11]) + ylab('')                                     + scale_x_log10(breaks = c(10,100,1000,10000))                                                                    + geom_point(data = count_PD_subject, aes(x = exp(EC50),     y = tumor_type), inherit.aes = FALSE, alpha = 1, shape = '|', size = 2)
count_PD_plots[[12]] <-  ggplot(data = count_PD_parms[count_PD_parms$parameter == count_PD_loc_unique[4],], aes(xdist = exp(dist_normal(mean, sd)), y = tumor_type))  + stat_slab(alpha = 0.6)    + geom_vline(linetype = 'dashed', xintercept = 1, alpha = 0.3) + theme_bw() + theme(panel.grid.minor = element_blank()) + ggtitle('Population distribution') + xlab(count_PD_names[12]) + ylab('')                                     + scale_x_log10(breaks = c(0.6,0.8,1,1.2,1.4,1.6, 1.8, 2.0))                                                      + geom_point(data = count_PD_subject, aes(x = exp(hill),     y = tumor_type), inherit.aes = FALSE, alpha = 1, shape = '|', size = 2)
count_PD_lattice <- ggarrange(count_PD_plots[[1]], count_PD_plots[[5]], count_PD_plots[[9]], count_PD_plots[[2]], count_PD_plots[[6]], count_PD_plots[[10]], count_PD_plots[[3]], count_PD_plots[[7]], count_PD_plots[[11]], count_PD_plots[[4]], count_PD_plots[[8]], count_PD_plots[[12]], ncol = 3, nrow = 4)
ggsave(count_PD_lattice, file = 'count_PD_lattice.svg', units = 'mm', width = 270, height = 360)

# Export the simulated dataset.
write.csv(sim_PD_data, file = 'simulated_pharmacodynamic_data_example.Rdata', row.names = FALSE)