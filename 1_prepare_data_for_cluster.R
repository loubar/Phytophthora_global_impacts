rm(list=ls())


library(dplyr)
library(MuMIn)
library(ape)
library(RNeXML)
library(phytools)


############################################################################################
# prepare the global spread and host range data
############################################################################################
# read in the geographic extent data  
geographic_extent <- read.csv("2018-04-17_impact_geographic_extent.csv",
                              stringsAsFactors = FALSE)

lat_range <- read.csv("impact_lat_range.csv",
                      stringsAsFactors = FALSE) %>% rename(species_name = species, impact_lat_max = max_lat)
# host range data
host_range <- read.csv("2018-04-18_impact_host_range.csv",
                       stringsAsFactors = FALSE)

# pathogenicity index from Liew
path_index <- read.csv("pathogenicity_index.csv",
                       stringsAsFactors = FALSE)



# disease symptom data
disease_symptoms <- read.csv("2019-01-28_disease_symptoms.csv",
                             stringsAsFactors = FALSE)[,1:3]



# trait database
traits <- read.csv("Phytophthora_trait_database.csv",
                   stringsAsFactors = FALSE,
                   na.strings = c("NA", "#N/A"))

# merge each of these with the trait database
# check species names match first
all(geographic_extent$species_name %in% traits$species_name)
all(lat_range$species_name %in% traits$species_name)
all(host_range$species_name %in% traits$species_name)
all(disease_symptoms$species_name %in% traits$species_name)



traits <- 
  left_join(traits, merge(lat_range, merge(geographic_extent, merge(disease_symptoms, host_range))))
# 1 + years known to science 
traits$years_known <- 2019 - traits$date
# number of countries sampled for binomial model of geographic extent
traits$n_trials <- length(unique(read.csv("2019-08-06_country_level_occurrence_database.csv",
                                          stringsAsFactors = FALSE)$iso3)) 


# choose the traits to model impact with
included_traits <- c("proliferation",
                     "caduceus",
                     "oospore_wall_index",
                     "chlamydospores",
                     "hyphal_swelling",
                     "root_disease",
                     "foliar_disease",
                     "oospores",
                     "gr_at_opt",
                     "temp_opt",
                     "temp_min")

impact_metrics <- c("impact_geographic_extent",
                    "impact_host_range",
                    "impact_lat_max")


cont_variables <- 
  names(which(sapply(traits[,included_traits], is.double)))

binary_variables <- 
  c(names(which(sapply(traits[,included_traits], is.character))), "root_disease", "foliar_disease")




# convert traits to binary variables for binomial modelling
for(i in grep("disease",
              binary_variables,
              value = TRUE,
              invert=TRUE)){
  traits[,i] <- ifelse(grepl("^N", traits[,i]), 
                       yes = 0,
                       no = 1)
}

raw_trait_data_subset <- traits[,c("species_name", "phylo_clade","years_known", included_traits, impact_metrics, "n_trials")]

write.csv(raw_trait_data_subset, file = "raw_trait_data.csv", row.names = FALSE)

# set the oospore wall index for species without oospores (sterile) to the mean value
traits$oospore_wall_index[which(traits$ho_he_s == "S")] <- mean(traits$oospore_wall_index, 
                                                                na.rm = TRUE) 


############################################################################################
# prepare the phylogenetic covariance matrix
############################################################################################
Phytoph_tree <-
  read.tree("updated - Posterior output.newick")
# for the multi-gene phylogeny, this would be 

# clean the tip labels to match the trait database
Phytoph_tree$tip.label <- gsub("Phytophthora_", "Phytophthora ", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("x_", "x ", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("austrocedrae",
                               "austrocedri",
                               Phytoph_tree$tip.label)
Phytoph_tree$tip.label[grep(" sp", Phytoph_tree$tip.label)] <- 
  sapply(grep(" sp", Phytoph_tree$tip.label),
         function(i)
           paste(strsplit(Phytoph_tree$tip.label, split = " sp")[[i]][1], 
                 paste0("'", strsplit(Phytoph_tree$tip.label, split = " sp")[[i]][2], "'"), 
                 collapse = " "))
Phytoph_tree$tip.label <- gsub("gondwanense", "gondwanensis", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("_", "", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("emzansi", "emanzi", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("bisheria", "bishii", Phytoph_tree$tip.label)

# check all species in the trait database are represented in the phylogeny
setdiff(traits$species_name, Phytoph_tree$tip.label)

Phytoph_tree <- multi2di(Phytoph_tree)
Phytoph_tree$edge.length[Phytoph_tree$edge.length==0] <- max(nodeHeights(Phytoph_tree))*0.001

Phytoph_tree <-
  chronopl(phy = Phytoph_tree,
           lambda = 1)
is.ultrametric(Phytoph_tree)

inv.phylo <- 
  MCMCglmm::inverseA(Phytoph_tree, nodes = "TIPS", scale = TRUE)
Phytophthora_covariance <- solve(inv.phylo$Ainv)
rownames(Phytophthora_covariance) <- rownames(inv.phylo$Ainv)

save(Phytophthora_covariance,
     file = "Phytophthora_covariance.rData")


#################################################################################
# also prepare the Martin et al. multi-gene phylogeny also for comparison of results #
#################################################################################

Phytoph_phylogeny <-
  nexml_read(x = "http://purl.org/phylo/treebase/phylows/tree/TB2:Tr65776?format=nexml")

# coerce the NeXML object to a phylo object recognised by
# ape, ade4, phylobase, adephylo
Phytoph_tree <- get_trees(Phytoph_phylogeny)

# prune away one of the citrophthora isolates as duplicated tip labels won't work when matching to the trait data
Phytoph_tree <- drop.tip(phy = Phytoph_tree, tip = c("Phytophthora citrophthora P10341", "Phytophthora palmivora arecae P10213"))

# remove isolate numbers
Phytoph_tree$tip.label <- gsub(" P$", "", gsub("[[:digit:]]", "", Phytoph_tree$tip.label))



# clean the tip labels to match the trait database
#Phytoph_tree$tip.label <- gsub("Phytophthora ", "P. ", Phytoph_tree$tip.label)

Phytoph_tree$tip.label <- gsub("alni", "x alni", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("andina", "x andina", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("katsurae", "castaneae", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("kelmania", "'kelmania'", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("lagoariana", "'lagoariana'", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("ohioensis", "'ohioensis'", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("personii", "'personii'", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("palmivora arecae", "palmivora", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("cinnamomi parvispora", "parvispora", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("salixsoil", "lacustris", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("sansomea", "sansomeana", Phytoph_tree$tip.label)
Phytoph_tree$tip.label <- gsub("bisheria", "bishii", Phytoph_tree$tip.label)

is.rooted(Phytoph_tree) # FALSE
#There are two outgroups, *Pythium undulata* and *Pythium vexans*,
# forming a polytomy at the node with the Phytophthora genus.
# If we drop the *Pythium vexans* tip, the tree should be properly rooted.
Phytoph_tree <- drop.tip(phy = Phytoph_tree, tip = "Pythium vexans")


is.ultrametric(Phytoph_tree) # FALSE
# No, we'll have to estimate the ages of the nodes by assuming branch lengths
# are proportional to the number of substitutions.  Is this a reasonable approach?
# Most of the trait evolution tools assume the phylogeny is ultrametric.
Phytoph_tree <- chronopl(phy = Phytoph_tree,
                         lambda = 1)


is.rooted(Phytoph_tree) # TRUE
is.ultrametric(Phytoph_tree) # TRUE

setdiff(traits$species_name, Phytoph_tree$tip.label)

# prepare the covariance matrix
inv.phylo <-
  MCMCglmm::inverseA(Phytoph_tree, nodes = "TIPS", scale = TRUE)
Phytophthora_covariance <- solve(inv.phylo$Ainv)
rownames(Phytophthora_covariance) <- rownames(inv.phylo$Ainv)

# how many species can we include in the models using multi-gene phylogeny?
length(intersect(Phytoph_tree$tip.label, traits$species_name))
# 78

save(Phytophthora_covariance,
     file = "Phytophthora_covariance_multigene.rData")

############################################################################################
# set up a list of candidate models for model comparison                                   #
############################################################################################

# For the individual trait models
# rescale continuous variables by 2 SDs as per Gelman 
# so they have variance 
# approximately the same as the binary variables 
# this means the parameter estimates will be comparable between binary and continuous variables
# the parameter estimates can then be interpreted as the number of units increase
# in impact (hosts, countries) with a 2SD increase in trait values

scale2sd <- function(var){
  return(var/(2*sd(var, na.rm = TRUE)))
}

# save the scaling factor so we can unscale the predictors later
get_2sd <- function(var){2 * sd(var, na.rm = TRUE)}
scaled_by <- sapply(traits[,c(cont_variables, "years_known")], get_2sd)
save(scaled_by, file = "scaled_by.rData")

traits[,c(cont_variables, "years_known")] <- 
  sapply(traits[,c(cont_variables, "years_known")], scale2sd)

# set up the global model predictors
# we'll need two sets, one with temp opt and one with temp min

predictors <- c("years_known", included_traits, "foliar_disease:root_disease")
# use dredge to get a list of all possible candidate models (subsets of predictors) from the global models with each set of predictors (individual traits, clusters and axes of trait-space)



#model_sets <- list()
dummy_traits <- as.data.frame(sapply(1:length(c("impact", included_traits, "years_known")),
                                     function(i) 
                                       rnorm(n = 20)))
colnames(dummy_traits) <- c("impact", included_traits, "years_known")
model_traits <- 
    lm(formula = as.formula(
      paste("impact",
            "~", 
            paste(predictors, 
                  collapse = " + "), 
            collapse = " ")),
      data = dummy_traits,
      na.action = "na.fail")

# model_traits <- 
#   lm(formula = as.formula(
#     paste("impact",
#           "~", 
#           paste(included_traits, 
#                 collapse = " + "), 
#           collapse = " ")),
#     data = dummy_traits,
#     na.action = "na.fail")


dredge_global_model <- 
    dredge(global.model = model_traits,
           evaluate = FALSE)
  
  
model_sets <-
    gsub(" \\+ 1", 
         "", 
         sapply(1:length(dredge_global_model),
                function(x) 
                  strsplit(split = "\\ ~ ",
                           as.character(dredge_global_model[[x]]))[[2]][2]))

  

all_predictors <-
  unique(unlist(model_sets))
#all_predictors <- all_predictors[grep("^1$", 
#                                      all_predictors,
#                                      invert = TRUE)]

# how many models to fit for model selection (with and without disease symptoms)
length(all_predictors_disease <- all_predictors)

# we'll always inlcude years known to science as a covariate, so exclude models that
# drop this variable
all_predictors_disease <- grep("years_known", 
                               all_predictors_disease, 
                               value = TRUE, 
                               invert = TRUE)
length(all_predictors_disease) # 2560 trait-based models
# 2^11 = 2048 subsets of individial traits
# plus 512 models (where both root and foliar disease predictors are present) with an interaction term between these variables
# 2048+512 = 2560 models in total

traits <- traits[,names(raw_trait_data_subset)]

save(traits,
     file = "traits.rData")
save(all_predictors_disease,
     file = "all_predictors_disease.rData")

# check variance inflation factors
library(car)
library(gpairs)
testmodel <- lm(impact_lat_max ~ years_known + caduceus + chlamydospores + foliar_disease + oospore_wall_index + oospores + proliferation + root_disease + temp_min + temp_opt + foliar_disease:root_disease,
                 data = traits)

vif_all_predictors <- vif(testmodel)
data.frame(predictor = names(vif_all_predictors), VIF = vif_all_predictors, row.names = NULL)

gpairs(traits[,included_traits])



# traits$obs <- 1:nrow(traits)
# prior <- c(prior(normal(0, 10), "b"),prior(normal(0, 50), "Intercept"),prior(student_t(4, 0, 1), "sd")) 
# testmodel <- brm(impact_lat_max ~ years_known + caduceus + chlamydospores + foliar_disease + oospore_wall_index + oospores + proliferation + root_disease + temp_opt + foliar_disease:root_disease + (1|species_name) + (1|obs),
#                   data = traits,
#                   prior = prior,
#                   cov_ranef = list(species_name = Phytophthora_covariance),
#                   iter = 5000,
#                   warmup = 4000, cores = 4,
#                  family = gaussian())
# 


testmodel <- lm(impact_lat_max ~ years_known + gr_at_opt + oospores + root_disease + temp_min,
                data = traits)




