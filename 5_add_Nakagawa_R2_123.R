
# caluclate the Nakagawa R2 values for each models and bind it to the model results
rm(list=ls())
library(brms)
library(abind)
library(snowfall)
# load the function to calculate r2glmm
source("get_Nakagawa_R2glmm.R")



#for each of the impact responses 
addR2 <- function(i){
  # load the model comparison object to bind the R2 results to
  load(grep("model_comparison", list.files(), value = TRUE)[i])
  # load the null model required by Nakagawa methods for variance partitioning
  load(paste0("null_model_",
              paste0(strsplit(grep("model_comparison_disease", 
                                   list.files(), 
                                   value = TRUE)[i], 
                              split = "_")[[1]][1:2], 
                     collapse = "_"),
              ".rData"))
  
  # get the model objects

  model_list <- grep("\\.csv", grep(paste0(strsplit(grep("model_comparison_disease",
                                          list.files(),
                                          value = TRUE)[i], 
                                     split = "_")[[1]][1:2],
                            collapse = "_"),
                     list.files(list.dirs()[2]),
                     value = TRUE), invert = TRUE, value = TRUE)
  # reorder to match the order of the model index in the model comparison table
  model_list <- 
    paste0(list.dirs()[2],
           "/",
           gsub("[[:digit:]]|\\.rData", 
                "",  
                model_list), 
           1:length(model_list),
           ".rData")
                           
  # collect the data in an array for binding to the model results
  R2_results <- 
    array(dim = c(length(model_list),7,3),
          dimnames = list(NULL, 
                          c("R2_fixed_m",
                            "R2_traits_m",
                            "R2_fixed_c",
                            "R2_traits_c",
                            "R2_phylo",
                            "ICC",
                            "ICCadj"),
                          c("est", "lower", "upper")))
  
  for(j in 1:length(model_list)){
    R2_results[j,,] <- 
      load_and_R2(model_name = model_list[j],
                  model_null = get(grep("null_model", ls(), value = TRUE)))
    print(j)
  }
  # bind it to the model results table
  assign(grep("model_comparison_disease", 
              ls(), 
              value = TRUE),
         abind(get(grep("model_comparison_disease", 
                        ls(), 
                        value = TRUE)), 
               R2_results, 
               along =2)
         )
  
  save(list =grep("model_comparison_disease", 
                ls(), 
                value = TRUE),
       file = paste0(grep("model_comparison_disease", 
                          ls(), 
                          value = TRUE), ".rData")
  )
  rm(list = c(grep("model_comparison_disease", 
                 ls(), 
                 value = TRUE),
              grep("null_model", ls(), value = TRUE),
              "model_list"))
}

# set up a connection to the cluster
hosts<-as.character(read.table(Sys.getenv('PBS_NODEFILE'),header=FALSE)[,1]) # read the nodes to use
sfSetMaxCPUs(length(hosts)) # ensure that snowfall can cope with this many hosts
sfInit(parallel=TRUE, type="MPI", cpus=length(hosts), useRscript=TRUE) # initialise the connection


# export to all processors the libraries, function and all data needed to run the model on each processor

sfLibrary(brms) 
sfLibrary(abind) 
sfExportAll()
sfLapply(1:3, 
         fun=addR2)
sfRemoveAll()
sfStop()
q(save = "no")
