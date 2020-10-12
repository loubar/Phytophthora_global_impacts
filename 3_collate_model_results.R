

rm(list = ls())

########################## geographic extent #######################


load("all_predictors_disease.rData")

length(missing_lat <- setdiff(1:length(all_predictors_disease), as.integer(gsub("\\D", "" ,grep("impact_lat_max", grep("\\.rData", list.files(), value=TRUE), value = TRUE)))))
length(missing_geo <- setdiff(1:length(all_predictors_disease), as.integer(gsub("\\D", "" ,grep("impact_geographic_extent", grep("\\.rData", list.files(), value=TRUE), value = TRUE)))))
length(missing_host <- setdiff(1:length(all_predictors_disease), as.integer(gsub("\\D", "" ,grep("impact_host_range", grep("\\.rData", list.files(), value=TRUE), value = TRUE)))))



# first check all the models have been fitted - if all models are present, both host_missing and geo_missing should be length 0
# If not, go back and rerun 3_fit_global_impact_models_disease.R
# It will automatically check which models have results already and only rereun the missing models
# Sometimes a node throws up an error I can't find the root of and the job is exited before all models have been fitted
# I think snowfall is a bit temperamental with any functions that have a built-in "cores = " argument, even if it isn't used
#length(gsub("\\D", "" ,grep("\\.rData", list.files("latitude_models"), value=TRUE)))
#setdiff(1:length(all_predictors_disease), as.integer(gsub("\\D", "" ,grep("\\.rData", list.files("latitude_models"), value=TRUE))))

#create an empty array for host range and geographic extent and latitudinal range to populate with the point estimates and lower and upper 95% credible intervals of the 
# parameter estimatesingle_gene
geographic_extent_model_comparison_disease <-
	array(dim = c(length(all_predictors_disease),
		      length(grep("object_name", colnames(read.csv("est_impact_geographic_extent_disease_1.csv")), value = TRUE, invert = TRUE)
),
		      3),
               dimnames = list(rep(NULL,length(all_predictors_disease)) , grep("object_name", colnames(read.csv("est_impact_geographic_extent_disease_1.csv")), value = TRUE, invert = TRUE)
, c("est", "lower", "upper")))
lat_max_model_comparison_disease <-host_range_model_comparison_disease <- geographic_extent_model_comparison_disease

########################### geographic extent #######################


for (i in c("est", "lower", "upper")){
	x <- do.call(rbind, lapply(paste0(grep(paste0(i, "_impact_geographic_extent_"), 
                      grep("\\.csv", 
                           list.files(), 
                           value = TRUE), 
                      value = TRUE)), 
                 read.csv,
		 stringsAsFactors = FALSE))
	z <- x[order(x$model_index),]
	z$object_name <- NULL
	geographic_extent_model_comparison_disease[,,i] <- as.matrix(z)[,dimnames(geographic_extent_model_comparison_disease)[[2]]]
}

save(geographic_extent_model_comparison_disease,
     file = "geographic_extent_model_comparison_disease.rData")


############################ host range ###############################

for (i in c("est", "lower", "upper")){
	x <- do.call(rbind, lapply(paste0(grep(paste0(i, "_impact_host_range_"), 
                      grep("\\.csv", 
                           list.files(), 
                           value = TRUE), 
                      value = TRUE)), 
                 read.csv,
		 stringsAsFactors = FALSE))
	z <- x[order(x$model_index),]
	z$object_name <- NULL
	host_range_model_comparison_disease[,,i] <- as.matrix(z)[,dimnames(host_range_model_comparison_disease)[[2]]]
}

save(host_range_model_comparison_disease,
     file = "host_range_model_comparison_disease.rData")

############################ latitude max ###############################

for (i in c("est", "lower", "upper")){
	x <- do.call(rbind, lapply(paste0(grep(paste0(i, "_impact_lat_max_"), 
                      grep("\\.csv", 
                           list.files(), 
                           value = TRUE), 
                      value = TRUE)), 
                 read.csv,
		 stringsAsFactors = FALSE))
	z <- x[order(x$model_index),]
	z$object_name <- NULL
	lat_max_model_comparison_disease[,,i] <- as.matrix(z)[,dimnames(lat_max_model_comparison_disease)[[2]]]
}

save(lat_max_model_comparison_disease,
     file = "lat_max_model_comparison_disease.rData")




