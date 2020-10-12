
get_R2_ICC <- function(model_f,        # fitted model with fixed effects (traits and years known)
                       model_null){    # null model with random effects only
  # calculation of the variance in fitted values
  # Nakagawa example seems to be variance on the linear predictor scale
  # this value includes the variance explained by years known to science
  require(brms)
  
  sigma2_fixed <- data.frame(apply(fitted(model_f, 
                                          scale = "linear", # linear predictor scale
                                          robust = TRUE, 
                                          re_formula = NA, 
                                          summary = FALSE), 1, var))
  # it would be useful to pull out the variance explained by traits only
  # the code below does this by fixing years known as a constant (the mean): marginl effects of traits
  # before calculating the variance
  sigma2_traits <- data.frame(apply(fitted(model_f, 
                                           scale = "linear", 
                                           robust = TRUE, 
                                           re_formula = NA, 
                                           summary = FALSE,
                                           newdata = data.frame(years_known = mean(model_f$data$years_known),
                                                                model_f$data[,grep("years_known",
                                                                                   colnames(model_f$data),
                                                                                   invert = TRUE)])), 1, var))
  
  # Vt is the total variance of the null model: includes phylogenetic random effect 
  # and observation-level random effect (overdispersion relative to binomial / Poisson)
  Vt <- (posterior_samples(model_null, pars = "sd_species_name")^2 ) + (posterior_samples(model_null, pars = "sd_obs")^2 )
  sigma2_phylo <- posterior_samples(model_f, pars = "sd_species_name")^2
  sigma2_phylo_null <- posterior_samples(model_null, pars = "sd_species_name")^2
  
  # Intercept of the null model: expected grand mean on the linear predictor scale
  # accounting for phylogenetically structured random intercept
  # and observation-level random intercepts (additive dispersion parameter)
  b0 <- as.numeric(fixef(model_null, summary = FALSE))
  
  # if the model is binomial (geographic extent)
  if(model_f$family$family == "binomial" ){
    n <- 159
    q <- b0 - 0.5 * Vt * tanh(b0 * (1 + 2 * exp(-0.5 * Vt))/6)
    # convert to the response scale using Eqn. 6.8 (Nakagawa et al. 2017)
    pmean <- plogis(q[,1])
    
    # residual variance 
    # of null model
    sigma2_epsilon <- as.numeric(1/(n*pmean * (1 - pmean)))
  }
  if(model_f$family$family == "poisson"){ # see model 5 in the paper
    # estimate lambda: grand mean (= variance) via Eqn. 5.8 
    lambda <- exp(b0 + 0.5*Vt)
    # sigma2_epsilon is ln(1 + 1/lambda)
    sigma2_epsilon <- log(1 + 1/lambda)
  }
  if(model_f$family$family == "gaussian"){ # see model 5 in the paper
    # then sigma_epsilon is just the total residual variance of the null model
    sigma2_epsilon <- Vt
  }

  
  #empty array for results (chose this format for easy binding to model results array)
  results <-
    array(dim = c(1,7,3),
          dimnames = list(NULL, 
                          c("R2_fixed_m",
                            "R2_traits_m",
                            "R2_fixed_c",
                            "R2_traits_c",
                            "R2_phylo",
                            "ICC",
                            "ICCadj"),
                          c("est", "lower", "upper")))
    
  # marginal (m) and conditional (c) R2
  # raw (ICC_phylo) and adjusted (ICCadj_phylo) intra-class correlation
  # return the summaries with the median estimate and lower and upper 95 credible intervals
  results[,"R2_fixed_m",]  <-  brms:::get_summary(robust = TRUE, sigma2_fixed/(sigma2_fixed + sigma2_phylo + sigma2_epsilon))[,c(1,3,4)]
  results[,"R2_traits_m",] <-  brms:::get_summary(robust = TRUE, sigma2_traits / (sigma2_fixed + sigma2_phylo + sigma2_epsilon))[,c(1,3,4)]
  results[,"R2_fixed_c",]  <-  brms:::get_summary(robust = TRUE, (sigma2_fixed + sigma2_phylo)/(sigma2_fixed + sigma2_phylo + sigma2_epsilon))[,c(1,3,4)]
  results[,"R2_traits_c",] <-  brms:::get_summary(robust = TRUE, (sigma2_traits + sigma2_phylo)/(sigma2_fixed + sigma2_phylo + sigma2_epsilon))[,c(1,3,4)]
  results[,"R2_phylo",]    <-  brms:::get_summary(robust = TRUE, sigma2_phylo/(sigma2_fixed + sigma2_phylo + sigma2_epsilon))[,c(1,3,4)]
  results[,"ICC",]         <-  brms:::get_summary(robust = TRUE, sigma2_phylo_null / (sigma2_phylo_null + sigma2_epsilon))[,c(1,3,4)] # raw proportion of all unexplained varinace that is phylogentically correlated 
  results[,"ICCadj",]      <-  brms:::get_summary(robust = TRUE, sigma2_phylo / (sigma2_phylo + sigma2_epsilon))[,c(1,3,4)] # proportion of unexplained variance that phylogenetically correlated, after adjusting for fixed effects

  return(results)
  
}

load_and_R2 <- function(model_name, model_null){
  load(model_name)
  get_R2_ICC(model_f = model_fit,
             model_null)
}


