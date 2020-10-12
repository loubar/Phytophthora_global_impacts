###############################################################################
# This script fits alternative Bayesian hierarchical models of global impact  #
# of 102 Phytophthora species, ready for model comparison.                    #
# Each model is fitted, and the                                               #
# results are saved as an R object to the cirrus directory                    #
# /home/pywell/loubar/Phytothreats/no_disease_trait_3LVs                      #
# The information criteria fro model comparison are calculated in a  second   #
# script (host_range_model_comparison.R, geographic_extent_model_comparsion.R)#   
# in the scame directory                                                      #
###############################################################################


rm(list = ls())
set.seed(42)




#file.remove(grep("[[:digit:]].rData", list.files(), value = TRUE))
library(brms)
library(snowfall)
library(Matrix)

# 1.  Load the required data
#nLVs <- 2
load("traits.rData")          # the processed trait database
load("Phytophthora_covariance.rData")  # Phylogenetic covariance matrix
load("all_predictors_disease.rData")            # list of predictors for each candidate mode
n_trials <- traits$n_trials[1]

# 2. get a list of all models to fit
##### missing models ################################
# the programme still trips up occsaionally, so for easy restarts
# get the list of models to fit based on what models are missing
# from the saved folder.  This should fit all models if the folder has no models saved in it.

missing_lat <- setdiff(1:length(all_predictors_disease), as.integer(gsub("\\D", "" ,grep("impact_lat_max", grep("\\.rData", list.files(), value=TRUE), value = TRUE))))
missing_geo <- setdiff(1:length(all_predictors_disease), as.integer(gsub("\\D", "" ,grep("impact_geographic_extent", grep("\\.rData", list.files(), value=TRUE), value = TRUE))))
missing_host <- setdiff(1:length(all_predictors_disease), as.integer(gsub("\\D", "" ,grep("impact_host_range", grep("\\.rData", list.files(), value=TRUE), value = TRUE))))

to_fit_lat <- expand.grid("impact_lat_max",
                      all_predictors_disease[missing_lat],
                      stringsAsFactors = FALSE)

to_fit_geo <- expand.grid("impact_geographic_extent",
                      all_predictors_disease[missing_geo],
                      stringsAsFactors = FALSE)


to_fit_host <- expand.grid("impact_host_range",
                      all_predictors_disease[missing_host],
                      stringsAsFactors = FALSE)


to_fit <- rbind(to_fit_lat,
		to_fit_geo,
                to_fit_host)

names(to_fit) <-
  c("response_var",
    "predictor_var")
 

# 3. Define a function to model impact as a function of the specified set of predictors in the dataframe to_fit
# model_i is an integer from 1:nrow(to_fit) 
fit_brms_global_model <- 
  function(model_i){
    response_var <- to_fit[model_i,"response_var"]
    predictor_var <- to_fit[model_i,"predictor_var"]
    predictor_index <- match(predictor_var, all_predictors_disease)
        # define the name of the model
    
	model_name <- 
    	paste0(response_var,
	       "_",
	       "disease",
               "_",
               predictor_index) 
     
        traits_disease <- traits[complete.cases(traits[,grep("\\:", c(unique(unlist(strsplit(all_predictors_disease[-1], split = " \\+ "))),
			    response_var), invert = TRUE, value = TRUE)]),]
    # define the LHS of the formula:
    # only geographic extent is a binomial variable
    # brms requires that the number of trials is specified
    # for binomial reponse variables
    if (grepl("geographic_extent", 
              response_var)){
      
      impact_variable <- 
        paste(response_var,
              " | trials(",
	      n_trials,
	      ")",
              collapse = ""
        )
      error_structure <- binomial(link = "logit")
    }
    # the other impact metrics are all modelled as 
    # negative binomial responses
    if (grepl("host_range", response_var)){
      impact_variable <- response_var      
	error_structure <- poisson(link = "log")

    }

    if (grepl("lat_max", response_var)){
      impact_variable <- response_var      
	error_structure <- gaussian()

    }


    model_formula <-
      as.formula(
        paste(
          impact_variable,
          "~",
	  "years_known",
	  "+",
	  predictor_var,
          "+",
          "(1|species_name)",
          "+",
          "(1|obs)",
          collapse = " "
        )
      )
      if(grepl("^1$", predictor_var)){
      	prior <- c(prior(normal(0, 50), "Intercept"),prior(student_t(4, 0, 1), "sd")) 
      } 
      if(!grepl("^1$", predictor_var)){
      	prior <- c(prior(normal(0, 10), "b"),prior(normal(0, 50), "Intercept"),prior(student_t(4, 0, 1), "sd")) 
      } 

   traits_disease$obs <- 1:nrow(traits_disease) 
 
    model_fit <-
    	brm(model_formula, 
            data = traits_disease, 
            family = error_structure[[1]],
            cov_ranef = list(species_name = Phytophthora_covariance),
	    prior = prior,
            control = list(adapt_delta = 0.999,
                           max_treedepth = 15),
            #cores = 4,
            iter = 5000,
	    warmup = 4000,
	    save_dso = TRUE)

          # get lambda estimate and significance (is it significantly greater than 0)
     hyp <- hypothesis(model_fit,
                       "sd_species_name__Intercept^2 / (sd_species_name__Intercept^2 + sd_obs__Intercept^2) = 0",
                       class = NULL)

     parests <-
          c(as.numeric(fixef(model_fit)[,"Estimate"]),
          brms:::get_summary(posterior_samples(model_fit, pars = "^sd_"))[,"Estimate"],
          hyp$hypothesis[,"Estimate"]
        )

     names(parests) <-
        c(grep("^b_|^sd_", 
               parnames(model_fit), 
               value = TRUE),
          "lambda")
     parests <- as.data.frame(t(parests))

     upper <-
         c(as.numeric(fixef(model_fit)[,"97.5%ile"]),
           brms:::get_summary(posterior_samples(model_fit, pars = "^sd_"))[,"97.5%ile"],
           hyp$hypothesis[,"u-95% CI"]
           )
     names(upper) <-
       names(parests)
     upper <- as.data.frame(t(upper))

     lower <-
        c(as.numeric(fixef(model_fit)[,"2.5%ile"]),
        brms:::get_summary(posterior_samples(model_fit, pars = "^sd_"))[,"2.5%ile"],
        hyp$hypothesis[,"l-95% CI"]
        )
     names(lower) <-
        names(parests) 
        lower <- as.data.frame(t(lower))

       model_info <-
       data.frame(
          object_name = model_name,
          model_index = predictor_index,
	  LOOIC = as.numeric(LOO(model_fit)$looic),
	  kfold10IC = as.numeric(kfold(model_fit)$kfoldic),
	  obs = nrow(model_fit$data),
	  stringsAsFactors = FALSE
      )
  
  # define a template data.frame which will be used to store results from each model before 
  # rbinding together for model ranking 
  # name the parameters as they will be called in brms
  parameters <- c("b_Intercept",
		"b_years_known", 
		"sd_obs__Intercept", 
		"sd_species_name__Intercept", 
		paste0("b_", unique(unlist(strsplit(all_predictors_disease[-1], split = " \\+ "))))
		)

par_names <- rep(0, length(parameters))
names(par_names) <- parameters


template_df <- as.data.frame(t(par_names))[0,]
	
  write.csv(
    merge(cbind(model_info, parests), 
          template_df,
          all.x = TRUE),
    file = paste0("est_",
                    model_name,
                    ".csv"),
    row.names = FALSE
  )
  
  write.csv(
    merge(cbind(model_info, upper), 
          template_df,
          all.x = TRUE),
    file = paste0("upper_",
                    model_name,
                    ".csv"),
    row.names = FALSE
  )

  write.csv(
    merge(cbind(model_info, lower), 
          template_df,
          all.x = TRUE),
    file = paste0("lower_",
                    model_name,
                    ".csv"),
    row.names = FALSE
  )

    
  save(model_fit,
       file = paste0(model_name,
                     ".rData"))
   rm(list = c("model_info", "model_fit"))
    
  }

# check it works
#fit_brms_global_model(model_i = 1)


# 9. set up a connection to the cluster
hosts<-as.character(read.table(Sys.getenv('PBS_NODEFILE'),header=FALSE)[,1]) # read the nodes to use
sfSetMaxCPUs(length(hosts)) # ensure that snowfall can cope with this many hosts
sfInit(parallel=TRUE, type="MPI", cpus=length(hosts), useRscript=TRUE) # initialise the connection


# 10. export to all processors the libraries, function and all data needed to run the model on each processor

sfLibrary(brms) 
sfExportAll()
sfLapply(1:nrow(to_fit), 
         fun=fit_brms_global_model)
sfRemoveAll()
sfStop()
q(save = "no")