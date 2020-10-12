# fit null models: intercepts only
# we can't use the intecpt only model from the model comparison set 
# as all our trait models use years known as a fixed covariate

load("Phytophthora_covariance.rData")
load("traits.rData")
traits$obs <- 1:nrow(traits)

# set priors
prior <- c(prior(normal(0, 50), "Intercept"),
           prior(student_t(4, 0, 1), "sd")) 


# fit the null models with only random effects
# our null models include fixed effect of years_known_to_science


# binomial model with additive dispersion
null_model_geographic_extent <-
  brm(impact_geographic_extent | trials(179) ~ 1 + (1 | species_name) + (1 | obs), 
      data = traits, 
      family = binomial(link = "logit"),
      cov_ranef = list(species_name = Phytophthora_covariance),
      prior = prior,
      control = list(adapt_delta = 0.999,
                     max_treedepth = 15),
      iter = 5000,
      warmup = 4000)
save(null_model_geographic_extent,
     file = "null_model_geographic_extent.rData")

# poisson model with additive dispersion
null_model_host_range <-
  brm(impact_host_range ~ 1 + (1 | species_name) + (1 | obs), 
      data = traits, 
      family = poisson(link = "log"),
      cov_ranef = list(species_name = Phytophthora_covariance),
      prior = prior,
      control = list(adapt_delta = 0.999,
                     max_treedepth = 15),
      iter = 5000,
      warmup = 4000)
save(null_model_host_range,
     file = "null_model_host_range.rData")

# gaussian model 
null_model_lat_max <-
  brm(impact_lat_max ~ 1 + (1 | species_name) + (1 | obs), 
      data = traits, 
      family = gaussian(),
      cov_ranef = list(species_name = Phytophthora_covariance),
      prior = prior,
      control = list(adapt_delta = 0.999,
                     max_treedepth = 15),
      iter = 5000,
      warmup = 4000)
save(null_model_lat_max,
     file = "null_model_lat_max.rData")

