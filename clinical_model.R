# AP. Woodward, University of Canberra, 2023.
# Generate a simulated multi-centre clinical trial of a bovine mastitis therapy.
#     The study is parallel and placebo-controlled with two active interventions, and two disease categories.
#     For simplicity, simple randomization is assumed.
#     The outcome is binomial (clinical cure or not), with the simulation and modelling conducted on logit scale.
# The analysis is focused on a Bayesian hierarchical generalized linear model with weakly informative priors and is implemented with 'brms' (https://doi.org/10.18637/jss.v080.i01).

# Load the required packages.
library(brms)
library(brmsmargins)
library(reshape2)
library(ggplot2)
library(ggdist)
library(ggpubr)
library(viridis)
library(scales)

# Set the total number of subjects.
n_subjects <- 200

# Implement the design matrix, with n_subjects rows and 9 columns.
# An intercept:
design_table <- data.frame(intercept = numeric(length = n_subjects))
design_table$intercept <- 1
theta_table  <- data.frame(parameter = character(length = 12), value = numeric(length = 12))

# A binary effect for infection type:
design_table$infection <- NA
design_table$infection <- rbinom(length(design_table$infection),1,2/3)

# Two binary effects for treatment type, and their interaction with infection type:
design_table$treatment1 <- NA
design_table$treatment2 <- NA
design_table$infection1_treatment1 <- NA
design_table$infection1_treatment2 <- NA

# Three treatments, each of equal allocation probability.
design_table[,3:4] <- (t(rmultinom(length(design_table$treatment1),1,c(1/3,1/3,1/3))))[,2:3]
design_table[,5]   <- as.numeric((design_table[2] == 1) & (design_table[3] == 1))
design_table[,6]   <- as.numeric((design_table[2] == 1) & (design_table[4] == 1))

# Six binary effects, representing six farms.
design_table[,7:12] <- NA
colnames(design_table)[7:12] <- c('farm1','farm2','farm3','farm4','farm5','farm6')
design_table[,7:12] <- (t(rmultinom(length(design_table$farm1),1,rep(1/6,6))))
design_matrix <- as.matrix(design_table)

# Set the coefficients.
# The intercept represents directly the control and infection status A.
logit_coefficients <- numeric(length = 12) 
logit_coefficients[1] <-  qlogis(0.3)

# Infection B has lower chance of success than A.
logit_coefficients[2] <- -0.5

# For Infection A, both treatments are effective, but Treatment1 is less effective than Treatment2.
logit_coefficients[3] <- 1
logit_coefficients[4] <- 1.5

# For Infection B, Treatment1 is similarly effective, but Treatment2 is less effective (actually a bit worse than control).
logit_coefficients[5] <- 0
logit_coefficients[6] <- -2

# There is a small difference in the outcome across farms (SD 0.3 on logit scale).
logit_coefficients[7:12] <- rnorm(6, 0, 0.3)

# Build the working simulated dataset.
sim_mastitis_data <-  data.frame(latent_prob = rowSums(sweep(design_matrix, 2, logit_coefficients, '*')))
sim_mastitis_data$infection <- 'A'
sim_mastitis_data$infection[design_table$infection == 1] <- 'B'
sim_mastitis_data$treatment <- 'control'
sim_mastitis_data$treatment[design_table$treatment1 == 1] <- 'treatment1'
sim_mastitis_data$treatment[design_table$treatment2 == 1] <- 'treatment2'
sim_mastitis_data$farm <- apply(design_table[,7:12], 1, function (x) which(x == 1))
sim_mastitis_data$cure <- rbinom(length(sim_mastitis_data$latent_prob),1,plogis(sim_mastitis_data$latent_prob))

# Generate a hierarchical GLM for cure, using weakly informative priors.
#     Obtain the posterior samples for the expected probability of cure by condition.
#     Note that for the Bayesian model, it's easy to obtain conditional risk and risk ratios just by working from the samples.
# Priors for the coefficients are set to N(0,2), which suggests that treatment and infection effects are most plausibly within 4 units of 0 on the logit scale.
#     The prior for the intercept is N(0,2) which suggests that the probability of cure on average over conditions is most plausibly 0.5, and probably not lower than about 0.12 or higher than 0.88.
cure_brm_fixef <- brm(cure ~ 1 + infection*treatment, family = bernoulli(), data = sim_mastitis_data, prior = set_prior('normal(0,2)', class = 'b') + set_prior('normal(0,2)', class = 'Intercept'))
summary(cure_brm_fixef)
cure_fix_design <- expand.grid(infection = c('A','B'), treatment = c('control', 'treatment1', 'treatment2'))
cure_fixef_marg <- brmsmargins(cure_brm_fixef, at = cure_fix_design, effects = 'fixedonly', CI = 0.9, seed = 5678)
cure_fixef_draws <- as.data.frame(cure_fixef_marg$Posterior)
colnames(cure_fixef_draws) <- paste(cure_fix_design$infection, cure_fix_design$treatment, sep = '_')
cure_fixef_draws$`A_treatment1:A_control`    <- cure_fixef_draws$A_treatment1/cure_fixef_draws$A_control
cure_fixef_draws$`A_treatment2:A_control`    <- cure_fixef_draws$A_treatment2/cure_fixef_draws$A_control
cure_fixef_draws$`A_treatment1:A_treatment2` <- cure_fixef_draws$A_treatment1/cure_fixef_draws$A_treatment2
cure_fixef_draws$`B_treatment1:B_control`    <- cure_fixef_draws$B_treatment1/cure_fixef_draws$B_control
cure_fixef_draws$`B_treatment2:B_control`    <- cure_fixef_draws$B_treatment2/cure_fixef_draws$B_control
cure_fixef_draws$`B_treatment1:B_treatment2` <- cure_fixef_draws$B_treatment1/cure_fixef_draws$B_treatment2
cure_fixef_long <- melt(cure_fixef_draws[,1:6])
cure_fixef_long$infection <- substr(cure_fixef_long$variable, 1, 1)
cure_fixef_long$treatment  <- 'NA'
cure_fixef_long$treatment[((cure_fixef_long$variable == 'A_control') | (cure_fixef_long$variable == 'B_control'))] <- 'control'
cure_fixef_long$treatment[((cure_fixef_long$variable == 'A_treatment1') | (cure_fixef_long$variable == 'B_treatment1'))] <- 'treatment1'
cure_fixef_long$treatment[((cure_fixef_long$variable == 'A_treatment2') | (cure_fixef_long$variable == 'B_treatment2'))] <- 'treatment2'
cure_fixef_plot <- ggplot(data = cure_fixef_long, aes(x = value, y = treatment)) + stat_halfeye(alpha = 0.6) + theme_bw() + facet_wrap(~infection) + scale_y_discrete(breaks = c('treatment1','treatment2','control'), labels = c('Treatment 1', 'Treatment 2', 'Control')) + coord_cartesian(xlim = c(0,1)) + xlab('Cure Probability') + ylab('') + ggtitle('No group (`random`) effects')
ggsave(cure_fixef_plot, file = 'cure_fixef_plot.svg', units = 'mm', height = 200, width = 250)
cure_fixef_contrasts <- melt(cure_fixef_draws[,7:12])
cure_fixef_contrasts$infection <- substr(cure_fixef_contrasts$variable, 1, 1)
cure_fixef_contrasts$contrast  <- 'NA'
cure_fixef_contrasts$contrast[((cure_fixef_contrasts$variable == 'A_treatment1:A_control') | (cure_fixef_contrasts$variable == 'B_treatment1:B_control'))] <- 'treatment1:control'
cure_fixef_contrasts$contrast[((cure_fixef_contrasts$variable == 'A_treatment2:A_control') | (cure_fixef_contrasts$variable == 'B_treatment2:B_control'))] <- 'treatment2:control'
cure_fixef_contrasts$contrast[((cure_fixef_contrasts$variable == 'A_treatment1:A_treatment2') | (cure_fixef_contrasts$variable == 'B_treatment1:B_treatment2'))] <- 'treatment1:treatment2'
cure_fixef_contrasts$contrast <- factor(cure_fixef_contrasts$contrast, levels = c('treatment1:control','treatment2:control','treatment1:treatment2'))
cure_fixef_contrasts_plot <- ggplot(data = cure_fixef_contrasts, aes(x = value, y = contrast)) + stat_halfeye(aes(fill = after_stat(x < 1)), alpha = 0.6) + theme_bw() + facet_wrap(~infection) + scale_y_discrete(breaks = c('treatment1:control','treatment2:control','treatment1:treatment2'), labels = c('Treatment 1: Control', 'Treatment 2: Control', 'Treatment 1: Treatment 2')) + scale_x_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) + coord_cartesian(xlim = c(0.01,100)) + geom_vline(xintercept = 1, linetype = 'dashed', alpha = 0.5) + scale_fill_viridis(discrete = TRUE, end = 0.9, option = 'cividis') + xlab('Cure Probability Ratio') + ylab('') + ggtitle('No group (`random`) effects') + theme(legend.position = 'none')
ggsave(cure_fixef_contrasts_plot, file = 'cure_fixef_contrasts_plot.svg', units = 'mm', height = 200, width = 250)
write.csv(round(fixef(cure_brm_fixef, robust = TRUE, probs = c(0.05,0.95)),2), file = 'cure_brm_fixef_linmod.csv')
fixef_theta_table <- unique(cure_fixef_long[,3:4])
fixef_theta_table$q5  <- NA
fixef_theta_table$q50 <- NA
fixef_theta_table$q95 <- NA
for(i in 1:(dim(fixef_theta_table)[1])){
  fixef_theta_table[i,3:5] <- round(quantile(cure_fixef_long$value[(cure_fixef_long$infection == fixef_theta_table$infection[i]) & (cure_fixef_long$treatment == fixef_theta_table$treatment[i])], c(0.05,0.50,0.95)),3)
}
write.csv(fixef_theta_table, file = 'fixef_cure_table.csv', row.names = FALSE)

# Generate a hierarchical GLM for cure, using weakly informative priors.
#     Obtain the posterior samples for the expected probability of cure by condition.
#     Note that for the Bayesian model, it's easy to obtain conditional risk and risk ratios just by working from the samples.
# Priors for the coefficients are set to N(0,2), which suggests that treatment and infection effects are most plausibly within 4 units of 0 on the logit scale.
#     The prior for the intercept is N(0,2) which suggests that the probability of cure on average over conditions is most plausibly 0.5, and probably not lower than about 0.12 or higher than 0.88.
#     This model additionally has a prior for the between-farm standard deviation on the logit scale; HN(0,1) specifies that it is most plausibly less than 1.
cure_brm_ranef <- brm(cure ~ 1 + infection*treatment + (1|farm), family = bernoulli(), data = sim_mastitis_data, prior = set_prior('normal(0,2)', class = 'b') + set_prior('normal(0,2)', class = 'Intercept') + set_prior('normal(0,0.5)', class = 'sd'))
summary(cure_brm_ranef)
cure_ranef_marg <- brmsmargins(cure_brm_ranef, at = cure_fix_design, effects = 'integrateoutRE', CI = 0.9, seed = 3456)
cure_ranef_draws <- as.data.frame(cure_ranef_marg$Posterior)
colnames(cure_ranef_draws) <- paste(cure_fix_design$infection, cure_fix_design$treatment, sep = '_')
cure_ranef_draws$`A_treatment1:A_control`    <- cure_ranef_draws$A_treatment1/cure_ranef_draws$A_control
cure_ranef_draws$`A_treatment2:A_control`    <- cure_ranef_draws$A_treatment2/cure_ranef_draws$A_control
cure_ranef_draws$`A_treatment1:A_treatment2` <- cure_ranef_draws$A_treatment1/cure_ranef_draws$A_treatment2
cure_ranef_draws$`B_treatment1:B_control`    <- cure_ranef_draws$B_treatment1/cure_ranef_draws$B_control
cure_ranef_draws$`B_treatment2:B_control`    <- cure_ranef_draws$B_treatment2/cure_ranef_draws$B_control
cure_ranef_draws$`B_treatment1:B_treatment2` <- cure_ranef_draws$B_treatment1/cure_ranef_draws$B_treatment2
cure_ranef_long <- melt(cure_ranef_draws[,1:6])
cure_ranef_long$infection <- substr(cure_ranef_long$variable, 1, 1)
cure_ranef_long$treatment  <- 'NA'
cure_ranef_long$treatment[((cure_ranef_long$variable == 'A_control') | (cure_ranef_long$variable == 'B_control'))] <- 'control'
cure_ranef_long$treatment[((cure_ranef_long$variable == 'A_treatment1') | (cure_ranef_long$variable == 'B_treatment1'))] <- 'treatment1'
cure_ranef_long$treatment[((cure_ranef_long$variable == 'A_treatment2') | (cure_ranef_long$variable == 'B_treatment2'))] <- 'treatment2'
cure_ranef_plot <- ggplot(data = cure_ranef_long, aes(x = value, y = treatment)) + stat_halfeye(alpha = 0.6) + theme_bw() + facet_wrap(~infection) + scale_y_discrete(breaks = c('treatment1','treatment2','control'), labels = c('Treatment 1', 'Treatment 2', 'Control')) + coord_cartesian(xlim = c(0,1)) + xlab('Cure Probability') + ylab('') + ggtitle('Group (`random`) effects integrated out')
ggsave(cure_ranef_plot, file = 'cure_ranef_plot.svg', units = 'mm', height = 200, width = 250)
cure_ranef_constrasts <- melt(cure_ranef_draws[,7:12])
cure_ranef_constrasts$infection <- substr(cure_ranef_constrasts$variable, 1, 1)
cure_ranef_constrasts$contrast  <- 'NA'
cure_ranef_constrasts$contrast[((cure_ranef_constrasts$variable == 'A_treatment1:A_control') | (cure_ranef_constrasts$variable == 'B_treatment1:B_control'))] <- 'treatment1:control'
cure_ranef_constrasts$contrast[((cure_ranef_constrasts$variable == 'A_treatment2:A_control') | (cure_ranef_constrasts$variable == 'B_treatment2:B_control'))] <- 'treatment2:control'
cure_ranef_constrasts$contrast[((cure_ranef_constrasts$variable == 'A_treatment1:A_treatment2') | (cure_ranef_constrasts$variable == 'B_treatment1:B_treatment2'))] <- 'treatment1:treatment2'
cure_ranef_constrasts$contrast <- factor(cure_ranef_constrasts$contrast, levels = c('treatment1:control','treatment2:control','treatment1:treatment2'))
cure_ranef_contrasts_plot <- ggplot(data = cure_ranef_constrasts, aes(x = value, y = contrast)) + stat_halfeye(aes(fill = after_stat(x < 1)), alpha = 0.6) + theme_bw() + facet_wrap(~infection) + scale_y_discrete(breaks = c('treatment1:control','treatment2:control','treatment1:treatment2'), labels = c('Treatment 1: Control', 'Treatment 2: Control', 'Treatment 1: Treatment 2')) + scale_x_log10(breaks = c(0.01,0.1,1,10,100), labels = c(0.01,0.1,1,10,100)) + coord_cartesian(xlim = c(0.01,100)) + geom_vline(xintercept = 1, linetype = 'dashed', alpha = 0.5) + scale_fill_viridis(discrete = TRUE, end = 0.9, option = 'cividis') + xlab('Cure Probability Ratio') + ylab('') + ggtitle('Group (`random`) effects integrated out') + theme(legend.position = 'none')
ggsave(cure_ranef_contrasts_plot, file = 'cure_ranef_contrasts_plot.svg', units = 'mm', height = 200, width = 250)
write.csv(round(rbind(fixef(cure_brm_ranef, robust = TRUE, probs = c(0.05,0.95)),VarCorr(cure_brm_ranef, robust = TRUE, probs = c(0.05,0.95))$farm$sd),2), file = 'cure_brm_ranef_linmod.csv')
VarCorr(cure_brm_ranef, probs = c(0.05,0.95))
ranef_theta_table <- unique(cure_ranef_long[,3:4])
ranef_theta_table$q5  <- NA
ranef_theta_table$q50 <- NA
ranef_theta_table$q95 <- NA
for(i in 1:(dim(ranef_theta_table)[1])){
  ranef_theta_table[i,3:5] <- round(quantile(cure_ranef_long$value[(cure_ranef_long$infection == ranef_theta_table$infection[i]) & (cure_ranef_long$treatment == ranef_theta_table$treatment[i])], c(0.05,0.50,0.95)),3)
}
write.csv(ranef_theta_table, file = 'ranef_cure_table.csv', row.names = FALSE)

# Aggregate and save the posterior distribution graphics.
cure_model_plots <- ggarrange(cure_fixef_plot, cure_ranef_plot, nrow = 2)
ggsave(cure_model_plots, file = 'cure_model_plots.svg', units = 'mm', height = 270, width = 220)
contrasts_model_plots <- ggarrange(cure_fixef_contrasts_plot, cure_ranef_contrasts_plot, nrow = 2)
ggsave(contrasts_model_plots, file = 'contrasts_model_plots.svg', units = 'mm', height = 270, width = 220)

# Summarize the probability of direction (https://doi.org/10.3389/fpsyg.2019.02767), as the simple probability that one treatment is superior (contrast >1).
length(which(cure_ranef_constrasts$value[((cure_ranef_constrasts$contrast == 'treatment1:control')    & (cure_ranef_constrasts$infection == 'A'))]>1))/ndraws(cure_brm_ranef)
length(which(cure_ranef_constrasts$value[((cure_ranef_constrasts$contrast == 'treatment2:control')    & (cure_ranef_constrasts$infection == 'A'))]>1))/ndraws(cure_brm_ranef)
length(which(cure_ranef_constrasts$value[((cure_ranef_constrasts$contrast == 'treatment1:treatment2') & (cure_ranef_constrasts$infection == 'A'))]>1))/ndraws(cure_brm_ranef)
length(which(cure_ranef_constrasts$value[((cure_ranef_constrasts$contrast == 'treatment1:control')    & (cure_ranef_constrasts$infection == 'B'))]>1))/ndraws(cure_brm_ranef)
length(which(cure_ranef_constrasts$value[((cure_ranef_constrasts$contrast == 'treatment2:control')    & (cure_ranef_constrasts$infection == 'B'))]>1))/ndraws(cure_brm_ranef)
length(which(cure_ranef_constrasts$value[((cure_ranef_constrasts$contrast == 'treatment1:treatment2') & (cure_ranef_constrasts$infection == 'B'))]>1))/ndraws(cure_brm_ranef)

# Obtain posterior samples of the absolute risk reduction (ARR), transform to the number-needed-to-treat (NNT), and plot the posterior distributions.
cure_fixef_ARR <- rbind(data.frame(ARR = cure_fixef_long$value[((cure_fixef_long$treatment == 'treatment1') & (cure_fixef_long$infection == 'A'))] - cure_fixef_long$value[((cure_fixef_long$treatment == 'control') & (cure_fixef_long$infection == 'A'))], treatment = 'treatment1', infection = 'A', model = 'fixef'),
                        data.frame(ARR = cure_fixef_long$value[((cure_fixef_long$treatment == 'treatment2') & (cure_fixef_long$infection == 'A'))] - cure_fixef_long$value[((cure_fixef_long$treatment == 'control') & (cure_fixef_long$infection == 'A'))], treatment = 'treatment2', infection = 'A', model = 'fixef'),
                        data.frame(ARR = cure_fixef_long$value[((cure_fixef_long$treatment == 'treatment1') & (cure_fixef_long$infection == 'B'))] - cure_fixef_long$value[((cure_fixef_long$treatment == 'control') & (cure_fixef_long$infection == 'B'))], treatment = 'treatment1', infection = 'B', model = 'fixef'),
                        data.frame(ARR = cure_fixef_long$value[((cure_fixef_long$treatment == 'treatment2') & (cure_fixef_long$infection == 'B'))] - cure_fixef_long$value[((cure_fixef_long$treatment == 'control') & (cure_fixef_long$infection == 'B'))], treatment = 'treatment2', infection = 'B', model = 'fixef'))
cure_ranef_ARR <- rbind(data.frame(ARR = cure_ranef_long$value[((cure_ranef_long$treatment == 'treatment1') & (cure_ranef_long$infection == 'A'))] - cure_ranef_long$value[((cure_ranef_long$treatment == 'control') & (cure_ranef_long$infection == 'A'))], treatment = 'treatment1', infection = 'A', model = 'ranef'),
                        data.frame(ARR = cure_ranef_long$value[((cure_ranef_long$treatment == 'treatment2') & (cure_ranef_long$infection == 'A'))] - cure_ranef_long$value[((cure_ranef_long$treatment == 'control') & (cure_ranef_long$infection == 'A'))], treatment = 'treatment2', infection = 'A', model = 'ranef'),
                        data.frame(ARR = cure_ranef_long$value[((cure_ranef_long$treatment == 'treatment1') & (cure_ranef_long$infection == 'B'))] - cure_ranef_long$value[((cure_ranef_long$treatment == 'control') & (cure_ranef_long$infection == 'B'))], treatment = 'treatment1', infection = 'B', model = 'ranef'),
                        data.frame(ARR = cure_ranef_long$value[((cure_ranef_long$treatment == 'treatment2') & (cure_ranef_long$infection == 'B'))] - cure_ranef_long$value[((cure_ranef_long$treatment == 'control') & (cure_ranef_long$infection == 'B'))], treatment = 'treatment2', infection = 'B', model = 'ranef'))
cure_pooled_ARR <- rbind(cure_fixef_ARR, cure_ranef_ARR)
cure_pooled_ARR$NNT <- 1/cure_pooled_ARR$ARR
trans_sign_log <- trans_new(name = 'signed log10', transform = function(x) sign(x)*log10(sign(x)*x), inverse = function(x) sign(x)*(10^(sign(x)*x)))
NTT_label_fun <- labeller(treatment = c(treatment1 = 'Treatment 1', treatment2 = 'Treatment 2'), infection = c(A = 'Infection A', B = 'Infection B'))
NTT_ranef_plot <- ggplot(data = cure_pooled_ARR[cure_pooled_ARR$model == 'ranef',], aes(x = NNT)) + stat_slab(aes(fill = after_stat(x < 0)), alpha = 0.6) + facet_grid(cols = vars(treatment), rows = vars(infection), labeller = NTT_label_fun) + theme_bw() + scale_fill_viridis(discrete = TRUE, end = 0.9, option = 'cividis') + scale_x_continuous(trans = trans_sign_log, breaks = c(-1000,-100,-10,-2,2,10,100,1000)) + theme(panel.grid.minor = element_blank(), legend.position = 'none') + scale_y_continuous(breaks = NULL) + ylab('') + xlab('Number needed to treat') + geom_vline(xintercept = 1, alpha = 0.4, linetype = 'dashed')
ARR_ranef_plot <- ggplot(data = cure_pooled_ARR[cure_pooled_ARR$model == 'ranef',], aes(x = ARR)) + stat_slab(aes(fill = after_stat(x < 0)), alpha = 0.6) + facet_grid(cols = vars(treatment), rows = vars(infection), labeller = NTT_label_fun) + theme_bw() + scale_fill_viridis(discrete = TRUE, end = 0.9, option = 'cividis') + xlim(-1,1) + theme(panel.grid.minor = element_blank(), legend.position = 'none') + scale_y_continuous(breaks = NULL) + ylab('') + xlab('Absolute risk reduction') + geom_vline(xintercept = 0, alpha = 0.4, linetype = 'dashed')
ARR_NTT_plot <- ggarrange(ARR_ranef_plot, NTT_ranef_plot, ncol = 1)
ggsave(ARR_NTT_plot, file = 'ARR_NTT_plot.svg', units = 'mm', height = 220, width = 200)

# Export the simulated dataset.
write.csv(sim_mastitis_data, file = 'simulated_clinical_data_example.Rdata', row.names = FALSE)