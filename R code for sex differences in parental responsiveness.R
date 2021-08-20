########### Sex Differences in Parental Responsiveness ###########
########### Shana Caro (shana_caro@post.harvard.edu)

rm(list=ls()) # clear R environment

#### Load packages ####

library(MCMCglmm)
library(mcmcplots)
library(ape)
library(coda)
library(mcgibbsit)
library(corrplot)
library(metafor)
library(lattice)
library(phylogram)
library(ggtree)
library(tidyverse)
library(lme4)

############## Import and Clean Data #############

# Download and save the data files "Full Dataset Excel" and the two nex files with the phylogenies. 
# Save the first sheet of the dataset as an excel or csv file, as you prefer
# Make sure you've put the data in your correct working directory

full_data <- read.csv("full_dataset.csv", header=T, na.strings="", strip.white=T) #import dataset
eric_30_trees <- read.nexus("eric_30_output.nex") #import random phylogenetic trees with Erikson backbone
hackett_30_trees <- read.nexus("eric_30_output.nex") #import random phylogenetic trees with Hackett backbone

full_data$study <- as.factor(full_data$study)
full_data$species <- as.factor(full_data$species)
full_data$reduction_y_n <- as.factor(full_data$brood_reduction_strategy)
full_data$animal <- as.factor(full_data$animal)
full_data$cooperative_breeding <- as.factor(full_data$cooperative_breeding)
full_data$social_bond_strength <- as.factor(full_data$social_bond_strength)
full_data$which_parent <- as.factor(full_data$which_parent)
full_data$environment <- factor(full_data$environment)
full_data$tarsus_dimorphism <- as.numeric(full_data$tarsus_dimorphism)
full_data$plumage_dimorphism <- as.numeric(full_data$plumage_dimorphism)
full_data$z_tarsus_dimorphism <- scale(full_data$tarsus_dimorphism, center=T, scale=T)
full_data$z_plumage_dimorphism <- scale(full_data$plumage_dimorphism, center=T, scale=T)
full_data$z_mean_perEPP <- scale(full_data$mean_perEPP, center=T, scale=T) #extrapair paternity
full_data$z_divorce_rate <- scale(full_data$divorce_rate, center=T, scale=T) #divorce rate

#convert effect size (correlation coefficient) to Z 
full_data$Z_whole_beg <- 0.5 * log((1+full_data$R_whole_beg)/(1-full_data$R_whole_beg), base=exp(1))
#calculate variance
full_data$variance <- 1 /(full_data$sample_size - 3)

length(full_data$Z_whole_beg) #156 effect sizes
length(levels(full_data$study)) #48 studies
length(levels(full_data$species)) #30 species

## exclude the species where we only have data from one sex

excl_data <- subset(full_data, species != "Pandion haliaetus" & species != "Strix aluco") 
  excl_data <- droplevels(excl_data)
length(excl_data$Z_whole_beg) #153 effect sizes
length(levels(excl_data$study)) #46 studies
length(levels(excl_data$species)) #28 species

#create a new dataframe by with species-sex averages
by_sex_species <- full_data %>%
  group_by(species, animal, which_parent, reduction_y_n, cooperative_breeding, social_bond_strength) %>%
  summarise(Z_whole_beg = mean(Z_whole_beg), R_whole_beg = mean(R_whole_beg), 
            mean.variance = mean(variance),  max.variance = max(variance), species.sample.size=sum(sample_size),
            z_divorce_rate = mean(z_divorce_rate, na.rm = TRUE), 
            z_mean_perEPP = mean(z_mean_perEPP, na.rm = TRUE),
            z_tarsus_dimorphism = mean(z_tarsus_dimorphism, na.rm = TRUE),
            z_plumage_dimorphism = mean(z_plumage_dimorphism, na.rm = TRUE),
            mean_male_baseline_corticosterone_during_young_care = mean(mean_male_baseline_corticosterone_during_young_care, na.rm = TRUE),
            mean_female_baseline_corticosterone_during_young_care = mean(mean_female_baseline_corticosterone_during_young_care, na.rm = TRUE))
data_sex_species <- as.data.frame(by_sex_species)
data_sex_species$variance <- 1/(data_sex_species$species.sample.size-3)

#create a new dataframe by with species' averages, only species with both sexes
by_species <- excl_data %>%
  group_by(species, animal, reduction_y_n, cooperative_breeding, social_bond_strength) %>%
  summarise(Z_whole_beg = mean(Z_whole_beg), R_whole_beg = mean(R_whole_beg), 
            mean.variance = mean(variance), species.sample.size=sum(sample_size),
            z_divorce_rate = mean(z_divorce_rate, na.rm = TRUE), 
            z_mean_perEPP = mean(z_mean_perEPP, na.rm = TRUE),
            z_tarsus_dimorphism = mean(z_tarsus_dimorphism, na.rm = TRUE),
            z_plumage_dimorphism = mean(z_plumage_dimorphism, na.rm = TRUE),
            mean_male_baseline_corticosterone_during_young_care = mean(mean_male_baseline_corticosterone_during_young_care, na.rm = TRUE),
            mean_female_baseline_corticosterone_during_young_care = mean(mean_female_baseline_corticosterone_during_young_care, na.rm = TRUE))
data_species_level <- as.data.frame(by_species)
data_species_level$variance <- 1/(data_species_level$species.sample.size-3)

#create columns for sex... I know there is probably a more elegant way to code this, but it gets the job done :)
one_row_species <- excl_data %>%
  group_by(species, which_parent) %>%
  summarise(Z_whole_beg = mean(Z_whole_beg))
data_one_row_species <- as.data.frame(one_row_species)
one_row_species_spread <- one_row_species  %>% spread(which_parent, Z_whole_beg)
data_one_row_species_spread <- as.data.frame(mutate(one_row_species_spread, diff.male.minus.female = male - female))

data_species_level$diff.male.minus.female <- data_one_row_species_spread$diff.male.minus.female
data_species_level$female.Z_whole_beg <- one_row_species_spread$female
data_species_level$male.Z_whole_beg <- one_row_species_spread$male

############## Metafor analyses #############

speciesSex_rmaZ <- rma(yi=Z_whole_beg, vi=variance ,  data=data_sex_species, slab = species)
  summary(speciesSex_rmaZ)
  regtest(speciesSex_rmaZ)
  funnel(speciesSex_rmaZ)

full_rmaZ <- rma(yi=Z_whole_beg, vi=variance, data=full_data, slab = species)
  summary(full_rmaZ)
  regtest(full_rmaZ)
  funnel(full_rmaZ, xlim=c(-1.6,1.6), cex=0.75, xlab="Z-transformed correlation coefficient")
  
  
diff_rmaZ <- rma.uni(yi=diff.male.minus.female, vi=variance , data=data_species_level, slab = species)
  summary(diff_rmaZ)
  regtest(diff_rmaZ)
  funnel(diff_rmaZ, xlim=c(-.8,.8), xlab="Within-species sex difference")


############## MCMCglmm analyses #############
# priors
  ## models reported in paper use informative priors, but uninformative priors can be seen below
  prior.informative1 <- list(R=list(V = 1, nu = 0.002), G=list(G1 = list(V=1,nu=2))) #high belief in phylogeny
  prior.informative3 <- list(R=list(V = 1, nu = 0.002), G=list(G1 = list(V=1,nu=2), G2 = list(V=1,nu=2), G3=list(V=1,nu=2))) #high belief in phylogeny, species and study
  
  prior.uninformative1 <- list(R=list(V = 1, nu = 0.002), G=list(G1 = list(V=1,nu=0.002))) #low belief in  phylogeny
  prior.uninformative3 <- list(R=list(V = 1, nu = 0.002), G=list(G1 = list(V=1,nu=0.002), G2 = list(V=1,nu=0.002), G3=list(V=1,nu=0.002))) #low belief in phylogeny, species and study

# The models below show the basic code for the models, using one random phylogenetic tree. In the manuscript, the reported results are the average over 20 random phylogenetic trees.   
# you can run a faster version of this model by reducing to (nitt=500000, burnin=200000, thin=100)
# you can change which random phylogenetic tree you use by eric_30_trees[[n]], where n < 100 , or by hackett_30_trees[[n]], where n < 100 

### overall sex difference ####

sex_diff_mod_erikson_1 <- MCMCglmm(Z_whole_beg ~ which_parent,
                                     random = ~ animal + study + species, 
                                     prior = prior.informative3, 
                                     pedigree = eric_30_trees[[1]], 
                                     mev = full_data$variance ,
                                     data = full_data,
                                     family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                     nitt=3000000, burnin=1000000, thin=1000)
  summary(sex_diff_mod_erikson_1)
sex_diff_mod_hackett_1 <- MCMCglmm(Z_whole_beg ~ which_parent -1, #this will show the mean for males and for females
                                     random = ~ animal, 
                                     prior = prior.informative1, 
                                     pedigree = hackett_30_trees[[47]], 
                                     mev = full_data$variance ,
                                     data = full_data,
                                     family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                     nitt=3000000, burnin=1000000, thin=1000)
  summary(sex_diff_mod_hackett_1)

### social bond models ####  
  sex_diff_mod_full <- MCMCglmm(Z_whole_beg ~ which_parent*social_bond_strength,
                                random = ~ animal + study + species, 
                                prior = prior.informative3, 
                                pedigree = hackett_30_trees[[73]], 
                                mev = full_data$variance,
                                data = full_data,
                                family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                nitt=3000000, burnin=1000000, thin=1000)
  summary(sex_diff_mod_full)
within_species_model <- MCMCglmm(diff.male.minus.female ~ social_bond_strength,
                                   random = ~ animal, 
                                   prior = prior.informative1, 
                                   pedigree = hackett_30_trees[[73]], 
                                   mev = data_species_level$variance,
                                   data = data_species_level,
                                   family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                   nitt=3000000, burnin=1000000, thin=1000)
  summary(within_species_model)
### uninformative prior models ####
    sex_diff_mod_full_uninform <- MCMCglmm(Z_whole_beg ~ which_parent*social_bond_strength,
                                         random = ~ animal + study + species, 
                                         prior = prior.uninformative3, 
                                         pedigree = hackett_30_trees[[73]], 
                                         mev = full_data$variance,
                                         data = full_data,
                                         family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                         nitt=3000000, burnin=1000000, thin=1000)
  summary(sex_diff_mod_full_uninform)  
  within_species_model_uninform <- MCMCglmm(diff.male.minus.female ~ social_bond_strength,
                                   random = ~ animal, 
                                   prior = prior.uninformative1, 
                                   pedigree = hackett_30_trees[[73]], 
                                   mev = data_species_level$variance,
                                   data = data_species_level,
                                   family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                   nitt=3000000, burnin=1000000, thin=1000)
  summary(within_species_model_uninform)
### environment models ####
    full_data$environment.ord <- factor(full_data$environment, ordered=TRUE)
  enviro_sex_diff_mod_full <- MCMCglmm(Z_whole_beg ~ which_parent*environment.ord + which_parent*reduction_y_n,
                                       random = ~ animal + study + species, 
                                       prior = prior.informative3, 
                                       pedigree = hackett_30_trees[[6]], 
                                       mev = full_data$variance,
                                       data = full_data,
                                       family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                       nitt=3000000, burnin=1000000, thin=1000)
  summary(enviro_sex_diff_mod_full)
enviro_sex_diff_mod_full2 <- MCMCglmm(Z_whole_beg ~ which_parent*social_bond_strength + which_parent*environment.ord + which_parent*reduction_y_n,
                                        random = ~ animal + study + species, 
                                        prior = prior.informative3, 
                                        pedigree = hackett_30_trees[[6]], 
                                        mev = full_data$variance,
                                        data = full_data,
                                        family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                        nitt=3000000, burnin=1000000, thin=1000)
  summary(enviro_sex_diff_mod_full2)
### tarsus models ####
tarsus_sex_diff_mod_full <- MCMCglmm(Z_whole_beg ~ which_parent*tarsus_dimorphism,
                                       random = ~ animal + study + species, 
                                       prior = prior.informative3, 
                                       pedigree = hackett_30_trees[[73]], 
                                       mev = subset(full_data, tarsus_dimorphism!="")$variance,
                                       data = subset(full_data, tarsus_dimorphism!=""),
                                       family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                       nitt=3000000, burnin=1000000, thin=1000)
  summary(tarsus_sex_diff_mod_full)
  tarsus_sex_diff_mod_within <- MCMCglmm(diff.male.minus.female ~ z_tarsus_dimorphism,
                                       random = ~ animal , 
                                       prior = prior.informative1, 
                                       pedigree = hackett_30_trees[[73]], 
                                       mev = subset(data_species_level, z_tarsus_dimorphism!="NaN")$variance,
                                       data = subset(data_species_level, z_tarsus_dimorphism!="NaN"),
                                       family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                       nitt=3000000, burnin=1000000, thin=1000)
  summary(tarsus_sex_diff_mod_within) 
  
### plumage models ####
### divorce model ####
  divorce_sex_diff_mod_full <- MCMCglmm(Z_whole_beg ~ which_parent*z_divorce_rate,
                                        random = ~ animal + study + species, 
                                        prior = prior.informative3, 
                                        pedigree = hackett_30_trees[[73]], 
                                        mev = subset(full_data, z_divorce_rate!="")$variance,
                                        data = subset(full_data, z_divorce_rate!=""),
                                        family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                        nnitt=3000000, burnin=1000000, thin=1000)
  summary(divorce_sex_diff_mod_full)
  within_divorce <- MCMCglmm(diff.male.minus.female ~ z_divorce_rate,
                             random = ~ animal, 
                             prior = prior.informative1, 
                             pedigree = hackett_30_trees[[6]], 
                             mev = subset(data_species_level, z_divorce_rate!="NaN")$variance,
                             data = subset(data_species_level, z_divorce_rate!="NaN"),
                             family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                             nitt=3000000, burnin=1000000, thin=1000)
  summary(within_divorce)
  
  
  
### epp models ####
  epp_sex_diff_mod_full <- MCMCglmm(Z_whole_beg ~ which_parent*z_mean_perEPP,
                                    random = ~ animal + study + species, 
                                    prior = prior.informative3, 
                                    pedigree = hackett_30_trees[[73]], 
                                    mev = subset(full_data, z_mean_perEPP!="")$variance,
                                    data = subset(full_data, z_mean_perEPP!=""),
                                    family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                    nitt=3000000, burnin=1000000, thin=1000)
  summary(epp_sex_diff_mod_full)
  eppBR_sex_diff_mod_full <- MCMCglmm(diff.male.minus.female ~ z_mean_perEPP,
                                      random = ~ animal, 
                                      prior = prior.informative1, 
                                      pedigree = hackett_30_trees[[73]], 
                                      mev = subset(data_species_level, z_mean_perEPP!="NaN")$variance,
                                      data = subset(data_species_level, z_mean_perEPP!="NaN"),
                                      family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                      nnitt=3000000, burnin=1000000, thin=1000)
  summary(eppBR_sex_diff_mod_full)

  eppdiovroce_sex_diff_mod_full <- MCMCglmm(Z_whole_beg ~ which_parent*z_mean_perEPP + which_parent*z_divorce_rate,
                                            random = ~ animal + study + species, 
                                            prior = prior.informative3, 
                                            pedigree = hackett_30_trees[[73]], 
                                            mev = subset(full_data, z_mean_perEPP!="" & z_divorce_rate!="")$variance,
                                            data = subset(full_data, z_mean_perEPP!="" & z_divorce_rate!=""),
                                            family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                            nitt=3000000, burnin=1000000, thin=1000)
  summary(eppdiovroce_sex_diff_mod_full)
  
  
  within_eppdiovroce_sex_diff_mod_full <- MCMCglmm(diff.male.minus.female ~ z_divorce_rate + z_mean_perEPP,
                                                   random = ~ animal, 
                                                   prior = prior.informative1, 
                                                   pedigree = hackett_30_trees[[73]], 
                                                   mev = subset(data_species_level, z_mean_perEPP!="NaN" & z_divorce_rate!="NaN")$variance,
                                                   data = subset(data_species_level, z_mean_perEPP!="NaN" & z_divorce_rate!="NaN"),
                                                   family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                                   nitt=3000000, burnin=1000000, thin=1000)
  summary(within_eppdiovroce_sex_diff_mod_full)
  
  INTER_within_eppdiovroce_sex_diff_mod_full <- MCMCglmm(diff.male.minus.female ~ z_divorce_rate * z_mean_perEPP,
                                                   random = ~ animal, 
                                                   prior = prior.informative1, 
                                                   pedigree = hackett_30_trees[[73]], 
                                                   mev = subset(data_species_level, z_mean_perEPP!="NaN" & z_divorce_rate!="NaN")$variance,
                                                   data = subset(data_species_level, z_mean_perEPP!="NaN" & z_divorce_rate!="NaN"),
                                                   family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                                   nitt=3000000, burnin=1000000, thin=1000)
  summary(INTER_within_eppdiovroce_sex_diff_mod_full)
  
### hormone models ####
  male_full_data <- subset(full_data, which_parent=="male")
  female_full_data <- subset(full_data, which_parent=="female")
  
  male_BC_model <- MCMCglmm(Z_whole_beg ~ mean_male_baseline_corticosterone_during_young_care,
                            random = ~ animal, 
                            prior = prior.informative1, 
                            pedigree = hackett_30_trees[[73]], 
                            mev = subset(male_full_data, mean_male_baseline_corticosterone_during_young_care!="")$variance,
                            data = subset(male_full_data, mean_male_baseline_corticosterone_during_young_care!=""),
                            family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                            nitt=3000000, burnin=1000000, thin=1000)
  summary(male_BC_model)
  male_BC_model_lm <- lmer(Z_whole_beg ~ mean_male_baseline_corticosterone_during_young_care + (1|common_name), 
                           data = subset(male_full_data, mean_male_baseline_corticosterone_during_young_care!=""))
  summary(male_BC_model_lm)
  
  female_BC_model <- MCMCglmm(female.Z_whole_beg ~ mean_female_baseline_corticosterone_during_young_care +
                                I(mean_female_baseline_corticosterone_during_young_care^2),
                              random = ~ animal, 
                              prior = prior.informative1, 
                              pedigree = hackett_30_trees[[73]], 
                              mev = subset(data_species_level, mean_female_baseline_corticosterone_during_young_care!="NaN")$variance,
                              data = subset(data_species_level, mean_female_baseline_corticosterone_during_young_care!="NaN"),
                              family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                              nitt=3000000, burnin=1000000, thin=1000)
  summary(female_BC_model)
  
  female_BC_model_lm <- lm(Z_whole_beg ~ mean_female_baseline_corticosterone_during_young_care 
                           +I(mean_female_baseline_corticosterone_during_young_care^2), 
                           data = subset(female_full_data,mean_female_baseline_corticosterone_during_young_care!="NaN"))
  summary(female_BC_model_lm)
  
#### phylogenetic signal  ####
  null_mod_phylo <- MCMCglmm(Z_whole_beg ~ 1,
                             random = ~ animal + study + species, 
                             prior = prior.informative3,
                             pedigree = hackett_30_trees[[73]], 
                             mev = full_data$variance,
                             data = full_data,
                             family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                             nitt=3000000, burnin=1000000, thin=1000)
  summary(null_mod_phylo)
###### Calculate heritability ####
  # I^2 values
  posterior.mode(null_mod_phylo$VCV)	
  animalVAR <- posterior.mode(null_mod_phylo$VCV[,1])			# phylogenetic variance 
  studyVAR <- posterior.mode(null_mod_phylo$VCV[,2])			# study variance
  speciesVAR <- posterior.mode(null_mod_phylo$VCV[,3])				# species variance
  mevVAR <- posterior.mode(null_mod_phylo$VCV[,4])				# mev variance
  unitsVAR <- posterior.mode(null_mod_phylo$VCV[,5])	 # units variance
  I2.phylo <- animalVAR / (studyVAR + speciesVAR  + animalVAR + mevVAR + unitsVAR) * 100	
  
#### potential methodological confounding factors models ####
  table(full_data$exp_obs) #sample sizes 
  table(full_data$estimated_test_statistic )
  table(full_data$beg_variable)
  table(full_data$feed_variable)
  table(full_data$beg_mode)
  table(full_data$brood_manipulation)
  table(full_data$deprivation)
  table(full_data$supplemented)
  table(full_data$playback)
  
exp_obs_1 <- MCMCglmm(Z_whole_beg ~ exp_obs,
                        random = ~ animal + study + species, 
                        prior = prior.informative3, 
                        pedigree = eric_30_trees[[1]], 
                        mev = full_data$variance ,
                        data = full_data,
                        family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                        nitt=3000000, burnin=1000000, thin=1000)
  summary(exp_obs_1)
  
beg_variable_1 <- MCMCglmm(Z_whole_beg ~ beg_variable,
                             random = ~ animal + study + species, 
                             prior = prior.informative3, 
                             pedigree = eric_30_trees[[1]], 
                             mev = full_data$variance ,
                             data = full_data,
                             family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                             nitt=3000000, burnin=1000000, thin=1000)
  round(summary(beg_variable_1)$solutions,2)
  
estimated_1 <- MCMCglmm(Z_whole_beg ~ estimated_test_statistic,
                          random = ~ animal + study + species, 
                          prior = prior.informative3, 
                          pedigree = eric_30_trees[[1]], 
                          mev = full_data$variance ,
                          data = full_data,
                          family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                          nitt=3000000, burnin=1000000, thin=1000)
  round(summary(estimated_1)$solutions,2)
  
beg_mode_1 <- MCMCglmm(Z_whole_beg ~ beg_mode,
                         random = ~ animal + study + species, 
                         prior = prior.informative3, 
                         pedigree = eric_30_trees[[1]], 
                         mev = full_data$variance ,
                         data = full_data,
                         family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                         nitt=3000000, burnin=1000000, thin=1000)
  round(summary(beg_mode_1)$solutions,2)
  
feed_variable_1 <- MCMCglmm(Z_whole_beg ~ feed_variable,
                              random = ~ animal + study + species, 
                              prior = prior.informative3, 
                              pedigree = eric_30_trees[[1]], 
                              mev = full_data$variance ,
                              data = full_data,
                              family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                              nitt=3000000, burnin=1000000, thin=1000)
  round(summary(feed_variable_1)$solutions,2)

deprivation_1 <- MCMCglmm(Z_whole_beg ~ deprivation,
                            random = ~ animal + study + species, 
                            prior = prior.informative3, 
                            pedigree = eric_30_trees[[1]], 
                            mev = full_data$variance ,
                            data = full_data,
                            family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                            nnitt=3000000, burnin=1000000, thin=1000)
  round(summary(deprivation_1)$solutions,2)
  
supplemented_1 <- MCMCglmm(Z_whole_beg ~ supplemented,
                             random = ~ animal + study + species, 
                             prior = prior.informative3, 
                             pedigree = eric_30_trees[[1]], 
                             mev = full_data$variance ,
                             data = full_data,
                             family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                             nitt=3000000, burnin=1000000, thin=1000)
  round(summary(supplemented_1)$solutions,2)

brood_manipulation_1 <- MCMCglmm(Z_whole_beg ~ brood_manipulation,
                                   random = ~ animal + study + species, 
                                   prior = prior.informative3, 
                                   pedigree = eric_30_trees[[1]], 
                                   mev = full_data$variance ,
                                   data = full_data,
                                   family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                                   nnitt=3000000, burnin=1000000, thin=1000)
  round(summary(brood_manipulation_1)$solutions,2)
  
playback_1 <- MCMCglmm(Z_whole_beg ~ playback,
                         random = ~ animal + study + species, 
                         prior = prior.informative3, 
                         pedigree = eric_30_trees[[1]], 
                         mev = full_data$variance ,
                         data = full_data,
                         family = "gaussian", verbose=FALSE, pr=TRUE, slice=TRUE,
                         nitt=3000000, burnin=1000000, thin=1000)
  round(summary(playback_1)$solutions,2)
  
