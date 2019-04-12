# R script to fit models of artifact class diversity
# across multiple geographic regions.

# Ryan Breslawski, rbreslawski@smu.edu
# Last edited on Mar 26, 2019

#    ++... ..... RUNNING THE SCRIPT ..... ...++
# 
# 1. Ensure that csv files geog_data.csv and 
#    point_data.csv are in the working directory.
#    
# 2. Run the script either line-by-line in
#    the source pane or by executing
#    'source("filename.R")' in the console pane.
#    


# Load libraries
library(parallel)
library(reshape2)
library(rstan)
# Allow Stan to save compiled model for use on each core
rstan_options(auto_write=TRUE)

###################################################
############ Set main arguments ##################
#################################################

# Set number of dataset permutations.
n_permutations <- 6

# Set number of iterations for each Stan model
# HMC chain. Position 1 is sampling iterations
# and position 2 is warmup iterations.
itrtns <- c(5000, 5000)

# Set adapt_delta parameter for Stan model
# HMC sampling.
adelta <- 0.95

# Define diversity index to be used:
# "H" Shannon Index
# "D" Simpson Index (1-D)
# "EH" Shannon evenness
# "Inv_D" Inverse Simpson Index (1/D)
div_arg <- "EH"


###################################################
############ Read and prepare data ###############
#################################################

geog_distr <- read.csv("geog_data.csv", header=TRUE, 
                        stringsAsFactors=FALSE)
point_distr <- read.csv("point_data.csv", header=TRUE, 
                        stringsAsFactors=FALSE)

# Convert sq km to 1000 sq km
geog_distr[,2:ncol(geog_distr)] <- geog_distr[,2:ncol(geog_distr)]/1000

# Find any cases where FIPS leading zeros were removed
geog_distr$county <- sapply(geog_distr$county, function(x){
  if(nchar(x)==1 | nchar(x)==4){
    return(paste0(0, x))
  }else{
    return(x)
  }
})
point_distr$county <- sapply(point_distr$county, function(x){
  if(nchar(x)==1 | nchar(x)==4){
    return(paste0(0, x))
  }else{
    return(x)
  }
})

 # Get zone and poitn class names
geog_zones <- colnames(geog_distr)[2:ncol(geog_distr)]
point_types <- colnames(point_distr)[2:ncol(point_distr)]

# Reformulate geog_distr such that it only has FIPS codes
# in the point data, and they are ordered the same way
geog_distr <- lapply(point_distr$county, function(x){
  geog_distr[which(geog_distr$county==x),]
})
geog_distr <- do.call("rbind", geog_distr)

# Simplify datasets. Now that they are aligned by county,
# remove the first column.
geog_distr <- as.matrix(geog_distr[,2:ncol(geog_distr)])
point_distr <- point_distr[,2:ncol(point_distr)]

## Create a function to calculate a user-specified 
## diversity index.
if(div_arg=="EH"){
  
  calc_div <- function(abundances){
  
    p <- abundances/sum(abundances)
    i <- -1*sum(sapply(p, function(j) {
      ifelse(j>0, j*(log(j)), 0)
    }))/log(length(p))
    
    return(i)
  }
    
}else if(div_arg=="H"){
  
  calc_div <- function(abundances){
    
    p <- abundances/sum(abundances)
    i <- -1*sum(sapply(p, function(j) {
      ifelse(j>0, j*(log(j)), 0)
    }))
    
    return(i)
  }
  
}else if(div_arg=="D"){
  
  calc_div <- function(abundances){
    
    p <- abundances/sum(abundances)  
    i <- 1 - sum(p^2)
    
    return(i)
  }
  
}else if(div_arg=="Inv_D"){
  
  calc_div <- function(abundances){
    
    p <- abundances/sum(abundances)  
    i <- 1/sum(p^2)
    
    return(i)
  }
}


########################################################
## Fit multinomial Stan model to only those counties ##
###### contained entirely within a geographic zone ###
#####################################################

# Get data frame of point counts for those points contained
# entirely within a geographic zone
point_distr_contained <- t(sapply(1:ncol(geog_distr), function(x){
  
  tempg <- geog_distr[, -x]
  
  if(ncol(geog_distr)>2){
    inds <- which(rowSums(tempg)==0)
  }else{
    inds <- which(tempg==0)
  }
  
  tempp <- point_distr[inds, ]
  
  return(colSums(tempp))
  
}))

# Create a summary list for the points contained entirely
# within geographic zones
summary_contained <- lapply(1:nrow(point_distr_contained), 
                                 function(u){
  
  # Subset points matrix
  sub_df <- point_distr_contained[u, ]
  points_in_dataset <- sum(point_distr_contained)
  classes_in_dataset <- length(which(colSums(point_distr_contained)>0))
  
  # Store in data.frame
  sdf <- data.frame(dataset="Counties within zones",
                    geog_zone=geog_zones[u],
                    diversity=calc_div(sub_df),
                    points_in_zone=sum(sub_df),
                    classes_in_zone=length(which(sub_df>0)),
                    points_in_dataset=points_in_dataset,
                    classes_in_dataset=classes_in_dataset,
                    stringsAsFactors=FALSE)
  colnames(sdf)[which(colnames(sdf)=="diversity")] <- div_arg
  
  return(sdf)
  
})

# Collapse summary list into data frame
summary_contained <- do.call("rbind", summary_contained)



###########################################################
######## Compile and fit Stan model ######################
#########################################################

point_distr_contained <- unname(point_distr_contained)

meds <- sapply(1:nrow(point_distr_contained), function(x){
  ind2 <- 2:(ncol(point_distr_contained)-1)
  m <- median(point_distr_contained[x, ])
  position <- which.min(abs(point_distr_contained[x, 
                                          ind2] - m))

  return(position + 1)
})

# Create data list for Stan model code compilation
d <- list(N_zones=length(geog_zones), 
          N_classes=length(point_types),
          p_counts=point_distr_contained,
          alpha=rep(0.05, length(point_types)),
          m=meds)

# Compile and fit Stan model
crs <- ifelse(detectCores() > 4, 4, detectCores())
Stan.model <- stan(file="point_diversity.stan", data=d, 
                   chains=4, iter=sum(itrtns), 
                   warmup=itrtns[2], cores=crs, 
                   control=list(max_treedepth=12,
                                adapt_delta=adelta))

# Get diagnostics
Stan.model.diagnostics <- summary(Stan.model)$summary
# Samples from model
Stan.model.samples <- extract(Stan.model)

rm(d, crs) # Remove uneeded variables from environment
cat("\nStan model compiled and fitted to point data from counties\n")
cat("contained entirely within individual geographic zones.\n\n")



######################################################################
###### Use geographic sampling distribution to simulate possible ####
###### point frequency distributions across the research question ##
###### defined geographical zones. ################################
##################################################################

# Create cluster
cl <- makeCluster(detectCores()-1)

# Export variables and packages to cluster
clusterExport(cl, varlist=c("point_types", "geog_zones",
                            "point_distr", "geog_distr",
                            "Stan.model.samples", 
                            "calc_div", "div_arg",
                            "itrtns", "adelta"))
clusterEvalQ(cl, {
  library(reshape2)
  library(rstan)
  })

# Sample points from counties to geographic zones until
# points are exhausted. Iterate n_permutations times.
sim_point_distrs <- parLapply(cl, 1:n_permutations, function(f){
  
  # Retrieve posterior sample
  p <- Stan.model.samples$theta[f, , ]
  
  # Simulate point counts across all counties
  all_points <- sapply(1:ncol(p), function(x){
    
    samples <- sapply(1:nrow(geog_distr), function(y){
      
      prbs <- geog_distr[y,]*p[,x]
      prbs <- prbs/sum(prbs)
      
      s <- sample(1:nrow(p), size=point_distr[y,x],
                  prob=prbs, replace=TRUE)
      s_table <- sapply(1:nrow(p), function(z){
        length(which(s==z))
      })
        
      return(s_table)

    })
    
    return(rowSums(samples))
    
  })
  
  # Create a summary list for the points contained entirely
  # within geographic zones
  summary_all_points <- lapply(1:nrow(all_points), 
                              function(u){
    
    # Subset points matrix
    sub_df <- all_points[u, ]
    points_in_dataset <- sum(all_points)
    classes_in_dataset <- length(which(colSums(all_points)>0))
    
    # Store in data.frame
    sdf <- data.frame(dataset=paste("Permutation", f),
                      geog_zone=geog_zones[u],
                      diversity=calc_div(sub_df),
                      points_in_zone=sum(sub_df),
                      classes_in_zone=length(which(sub_df>0)),
                      points_in_dataset=points_in_dataset,
                      classes_in_dataset=classes_in_dataset,
                      stringsAsFactors=FALSE)
    colnames(sdf)[which(colnames(sdf)=="diversity")] <- div_arg
    
    return(sdf)
    
  })
  
  # Collapse summary list into data frame
  summary_all_points <- do.call("rbind", summary_all_points)
  
  return(list(all_points, summary_all_points))
  
})

# Extract summaries of data permutations from list
datasets_summary <- lapply(1:length(sim_point_distrs),
                           function(x) sim_point_distrs[[x]][[2]])
datasets_summary <- do.call("rbind", datasets_summary)
datasets_summary <- rbind(summary_contained, datasets_summary)

# Reduce list to simulated point distributions
sim_point_distrs <- lapply(1:length(sim_point_distrs),
                           function(x) sim_point_distrs[[x]][[1]])

cat(paste(n_permutations, "dataset permutations generated.\n\n"))

###########################################################
######## Fit Stan model to data permutations #############
#########################################################

cat(paste0("Fitting model to ", n_permutations, " datasets (", 
           Sys.time(), ")...\n"))

# Fit a model to each simulated dataset
clusterExport(cl, varlist=c("sim_point_distrs", "Stan.model"))
Sim_list <- parLapply(cl, 1:n_permutations, function(t){

  
  meds <- sapply(1:nrow(sim_point_distrs[[t]]), function(x){
    ind2 <- 2:(ncol(sim_point_distrs[[t]])-1)
    m <- median(sim_point_distrs[[t]][x, ])
    position <- which.min(abs(sim_point_distrs[[t]][x, ind2] - m))
    
    return(position + 1)
  })
  
  # Create data list for Stan model code compilation
  d <- list(N_zones=length(geog_zones), 
            N_classes=length(point_types),
            p_counts=sim_point_distrs[[t]],
            alpha=rep(0.05, length(point_types)),
            m=meds)
  
  # Fit stan model to data d using compiled stan code
  rstan_options(auto_write=FALSE)
  returned_model <- stan(fit=Stan.model, data=d, 
                         chains=4, verbose=FALSE, 
                         iter=sum(itrtns), warmup=itrtns[2],
                         cores=1, save_warmup=FALSE,
                         control=list(adapt_delta=adelta, 
                                      max_treedepth=12))
  
  # Get model diagnostics
  ms <- summary(returned_model)$summary
  simpars <- suppressWarnings(get_sampler_params(returned_model))
  simpars <- do.call("rbind", simpars)
  Model_Diagnostics <- data.frame(permutationID=t,
                       Rhat_max=max(ms[which(!is.nan(ms[,10])),10]),
                       neff_min_perc=round(25*min(ms[which(!is.nan(ms[,9])),
                                                     9])/itrtns[1], 2),
                       max_treedepth=max(simpars[,3]),
                       Rhat_over_1.01=length(which(ms[,10]>1.01)),
                       neff_under_10perc=length(which(ms[,9]<itrtns[1]*0.4)),
                       n_exceed_max_treedepth=length(which(simpars[,3]>12)),
                       n_divergent=sum(simpars[,5]))
  
  # Extract 500 posterior point class probs from model
  dx <- extract(returned_model)$theta[1:500, , ]
  
  # Get point diversity for each zone
  ze <- lapply(1:length(geog_zones), function(z){
    
    # Extract posterior probs for zone z
    dxz <- dx[, z, ]
    
    ### Calculate diversity across posterior samples
    He <- sapply(1:nrow(dxz), function(i) calc_div(dxz[i,]))
    
    # Store as a data frame with geographic zone and dataset
    # permutation information.
    He_DF <- data.frame(diversity=He, 
                        Geog_Zone=rep(geog_zones[z], length(He)),
                        Dataset_Permutation=rep(t, length(He)),
                        stringsAsFactors=FALSE)
    colnames(He_DF)[which(colnames(He_DF)=="diversity")] <- div_arg
    return(He_DF)
  })
  
  div <- do.call("rbind", ze)
  
  return(list(returned_model, div, Model_Diagnostics))
  
})

# End cluster
stopCluster(cl)

cat(paste0("All model fits successful (", Sys.time(), ").\n\n"))



# Save Stan models in separate chunks for memory
# memory limits on personal computers. Each file 
# name ends with an integer indexing the first
# model in the saved object.
sll <- length(Sim_list)
for(i in seq(1, n_permutations, 100)){
  ifelse(i==1, l0s <- "_000", 
         ifelse(i<1000, l0s <- "_0", l0s <- "_"))
  if(i + 49 < sll){
    M <- lapply(seq(i, i+99, 1), 
                function(x) Sim_list[[x]][[1]])
    save(M, file=paste0("Mdls", l0s, i, ".RData"))    
  }
  else{
    M <- lapply(seq(i, sll, 1), 
                function(x) Sim_list[[x]][[1]])
    save(M, file=paste0("Mdls", l0s, i, ".RData"))
  }
}
rm(i, M, l0s, sll) # Remove extraneous variables

cat("Model fits saved to disk.\n")



# Extract diversity values calculated from posterior samples
Post_Div_Samples <- lapply(1:length(Sim_list), 
                           function(x) Sim_list[[x]][[2]])
Post_Div_Samples <- do.call("rbind", Post_Div_Samples)

# Extract model diagnostics
Model_diagnostics <- lapply(1:length(Sim_list), 
                           function(x) Sim_list[[x]][[3]])
Model_diagnostics <- do.call("rbind", Model_diagnostics)
Model_diagnostics$Potential_Problems <- 
  sapply(1:nrow(Model_diagnostics), function(x){
    ifelse(sum(Model_diagnostics[x,5:ncol(Model_diagnostics)])>0,
           1, 0)    
  }) 

# Sample posterior model parameters for one data permutation
sample_permutation <- sample(size=1, 1:n_permutations)
sample_model <- list(paste("Data permutation", sample_permutation),
                     extract(Sim_list[[sample_permutation]][[1]]))

# Save output
save(Model_diagnostics, Post_Div_Samples, sim_point_distrs,
     sample_model, Stan.model, Stan.model.diagnostics,
     Stan.model.samples, datasets_summary, div_arg,
     calc_div, point_types, file="Output_and_Diagnostics.RData")

cat("Posterior diversity values calculated from samples and saved to disk.\n")
