# He aprendido mucho pero la movida es que: igual no merece la pena
# hacerlo así, sino que es mejor dejar el while y DENTRO DEL WHILE llamar a
# una funcíón que chequee el if y haga el draw   --->   crear una función c++
#                                                       con while = inf loop
# .
# Y precisamente el sample lo he dejado con código de R...

setwd("/home/silvia/AAA/2021-06-28_my_null_model/2022-10-04_my_null_model_V2/")
library("untb") # simulate_timeseries uses "as.count"
library("tidyverse")
source("/home/silvia/Apps/my_functions.R")
library("gsubfn")
source("/home/silvia/Apps/functions_for_neutral_modelling.R")
library("parallel")
library("optparse")
library("dplyr")
library("Rcpp")


Rcpp::sourceCpp("./growth.cpp")

load("borrar.RData") # PCG_counts

# function without c++
simulate_timeseries_old <- function (counts_data,
                                 dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                 no_of_dil=12,
                                 abun_total=NULL,
                                 grow_step=1,
                                 keep_all_timesteps=FALSE,
                                 keep_going=FALSE){
  # This function simulates the changes in abundances for a given initial 
  # abundance data. The model has a dilution and growth system. After each 
  # dilution, a random organism is duplicated. This happens, by default, until 
  # the total abundance reaches the value previous to the dilution. But it can 
  # be changed with the option abun_total. "grow_step" is the number of 
  # individuals that grow each timestep.
  
  # keep_all_timesteps==FALSE saves a bit of RAM if there are many dilutions.
  # Only advised to make it TRUE for plotting.
  
  # If keep_going==TRUE, just pick any random bug in case the dilution factor is
  # so strong everything goes extinct. Prints a warning. If FALSE, simply exit.
  
  start <- as.count(my_transpose(counts_data))
  start <- start[order(names(start))] # this order thing is to keep indices consistent
  
  # Preparo la matriz que iré rellenando
  len <- length (start)    # número de especies
  if (is.null(abun_total)){ # establish default value.
    abun_total <- sum(start) # número de individuos al que dejaremos crecer antes de volver a diluir cada vez
  }
  
  # initialize data.frame with starting abundances
  if (keep_all_timesteps){
    trajectory <- matrix(NA, ncol=len, nrow = no_of_dil+1)
    rownames(trajectory) <- 0:no_of_dil
    colnames(trajectory) <- names(start)
    trajectory["0",]=start
  } 
  
  # y el "zero counter"
  empty <- start
  empty[1:length(empty)] <- 0
  
  this_timestep <- start
  
  message(paste0("Original abundance: ", sum(this_timestep)))
  message(paste0("After dilution: ", dilution*sum(this_timestep)))
  
  for (i in 1:no_of_dil){
    # Hago la dilución (diluimos el inóculo inicial también, como en el experimento original)
    
    # new_bugs <-c(1:length(this_timestep))[.Internal(sample(
    #   x = length(this_timestep),
    #   size = sum(this_timestep)*dilution, 
    #   replace = TRUE,
    #   prob = this_timestep))]
    # 
    # for (i in new_bugs) {
    #   this_timestep[new_bugs] <- this_timestep[new_bugs] + 1 # añado al que ha nacido al censo
    # }
    
    ## CASE 1 -- All taxa extinct or almost
    if (trunc(sum(this_timestep)*dilution)==0) {
      # feb 13 2022: step to 1, and pick a random single bug anyway
      message("WARNING: less than 0.5 bugs were left after diluting! Consider changing your dilution factor.")
      if (keep_going) {
        message("Picking one random bug and running simulation anyway...")   
        i <- c(1:length(this_timestep))[.Internal(sample(
          x = length(this_timestep),
          size = 1,
          replace = FALSE,
          prob = this_timestep))]
        this_timestep[i] <- this_timestep[i] + 1
      } else {
        stop("EXIT: all bugs extinct, nothing to simulate.")
      }
      
      ## If not CASE 1, dilute normally
    } else {
      this_timestep <- table(names(this_timestep)[.Internal(sample(
        x = length(this_timestep),
        size = sum(this_timestep)*dilution,
        replace = TRUE,
        prob = this_timestep))])
      # indicamos los ceros que hayan aparecido
      temp <- empty                                                                #0
      temp[names(this_timestep)] <- this_timestep                                  #0 # aquí se desordena, pero solo hasta la siguiente línea.
      this_timestep <- temp[order(names(temp))]                                    #0
    }
    
    # (1) hago que los bichos se dupliquen al azar mediante sampling hasta 
    # llegar a la cantidad inicial (while loop); (2) run the substitution
    while (sum(this_timestep) < abun_total) {# si aún debe seguir creciendo
      if (sum(this_timestep) < grow_step) { # the step is too big, reduce it
        step <- max(trunc(sum(this_timestep)/2))
        # warning(paste("Step reduced to",step))
      } else {
        step <- grow_step
      }
      
      ## Avoid growing too much (when step>1)
      if ((sum(this_timestep)+step) > abun_total) {
        step=(abun_total-sum(this_timestep))}
      
      ## Once the step is decided, sample
      new_bugs <- c(1:length(this_timestep))[.Internal(sample(
        x = length(this_timestep),
        size = step, 
        replace = TRUE,
        prob = this_timestep))]
      for (i in new_bugs) {
        this_timestep[i] <- this_timestep[i] + 1
      } # añado al que ha nacido al censo
      
    }
    
    # una vez crecidos, se puede diluir de nuevo: se repite el bucle
    # pero primero "secuenciamos" (guardamos las abundancias actualizadas)
    if (keep_all_timesteps){
      trajectory[as.character(i),]=this_timestep[colnames(trajectory)] # order by "trajectory" column names just in case
    }
  }
  
  message(paste0("Reached final abundance of: ", sum(this_timestep)))
  
  if (keep_all_timesteps){
    return(trajectory)
  } else {
    return(this_timestep%>%setNames(names(start)))
  }
}

# functions with c++

simulate_timeseries_new <- function (counts_data,
                                     dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                     no_of_dil=12,
                                     fixation_at=1,
                                     abun_total=NULL,
                                     grow_step=1,
                                     keep_all_timesteps=FALSE,
                                     force_continue=FALSE) { 
  
  start <- as.count(my_transpose(counts_data))
  start <- start[order(names(start))] # this order thing is to keep indices consistent
  
  # Preparo la matriz que iré rellenando
  len <- length (start)    # número de especies
  if (is.null(abun_total)){ # establish default value.
    abun_total <- sum(start) # número de individuos al que dejaremos crecer antes de volver a diluir cada vez
  }
  
  # y el "zero counter"
  empty <- start
  empty[1:length(empty)] <- 0
  
  this_timestep <- start
  dil <- 0
  while (dil < no_of_dil &
         (max(this_timestep)/sum(this_timestep) < fixation_at)){
    ## CASE 1 -- All taxa extinct or almost
    if (trunc(sum(this_timestep)*dilution)==0) {
      if (force_continue) {
        message("WARNING: less than 1 bug left after diluting! Consider changing your dilution factor.")
        message("Picking one random bug and running simulation anyway...")   
        i <- c(1:length(this_timestep))[.Internal(sample(
          x = length(this_timestep),
          size = 1,
          replace = FALSE,
          prob = this_timestep))]
        this_timestep    <- empty
        this_timestep[i] <- this_timestep[i] + 1
      } else {
        write("Cancelled the simulations. Less than 1 bug left after diluting! Consider changing your dilution factor.", stdout())
        stop("Cancelled the simulations. Less than 1 bug left after diluting! Consider changing your dilution factor.")
      }
      ## If not CASE 1, dilute normally
    } else {
      ## CASE 2 -- throw a warning if not many bugs left
      if (sum(this_timestep) <= 3) {
        write("WARNING: 3 or less bugs were left after diluting! Consider changing your dilution factor.", stdout())  # TODO 3 es ridículo, pensar uno más grande
        write("Running simulation anyway...", stdout())
      }
      this_timestep <- table(names(this_timestep)[.Internal(sample(
        x = length(this_timestep),
        size = sum(this_timestep)*dilution,
        replace = TRUE,
        prob = this_timestep))]) 
      
      temp <- empty                                                                #0
      temp[names(this_timestep)] <- this_timestep                                  #0 # aquí se desordena, pero solo hasta la siguiente línea.
      this_timestep <- temp[order(names(temp))]                                    #0
    }
    
    # (1) hago que los bichos se dupliquen al azar mediante sampling hasta 
    # llegar a la cantidad inicial (while loop); (2) run the substitution
    ns <- names(this_timestep)
    this_timestep <- as.vector(this_timestep)
    while (sum(this_timestep) < abun_total) {
      this_timestep <- growth(this_timestep, abun_total, grow_step)               # Rcpp function
    }
    names(this_timestep) <- ns
    dil <- dil + 1
  }
  if (max(this_timestep)/sum(this_timestep) >= fixation_at) {
    write(paste0(fixation_at*100, "% fixation of ",
                 names(this_timestep)[this_timestep!=0], " after ", dil, " dilutions."), stdout())
  }
  return(this_timestep)
}

simulate_timeseries_newer <- function (counts_data,
                                 dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                 no_of_dil=12,
                                 fixation_at=1,
                                 abun_total=NULL,
                                 grow_step=1,
                                 keep_all_timesteps=FALSE,
                                 force_continue=FALSE) { 
  
  start <- as.count(my_transpose(counts_data))
  start <- start[order(names(start))] # this order thing is to keep indices consistent
  
  # Preparo la matriz que iré rellenando
  len <- length (start)    # número de especies
  if (is.null(abun_total)){ # establish default value.
    abun_total <- sum(start) # número de individuos al que dejaremos crecer antes de volver a diluir cada vez
  }
  
  # y el "zero counter"
  empty <- start
  empty[1:length(empty)] <- 0
  
  this_timestep <- start
  dil <- 0
  while (dil < no_of_dil &
         (max(this_timestep)/sum(this_timestep) < fixation_at)){
    ## CASE 1 -- All taxa extinct or almost
    if (trunc(sum(this_timestep)*dilution)==0) {
      if (force_continue) {
        message("WARNING: less than 1 bug left after diluting! Consider changing your dilution factor.")
        message("Picking one random bug and running simulation anyway...")   
        i <- c(1:length(this_timestep))[.Internal(sample(
          x = length(this_timestep),
          size = 1,
          replace = FALSE,
          prob = this_timestep))]
        this_timestep    <- empty
        this_timestep[i] <- this_timestep[i] + 1
      } else {
        write("Cancelled the simulations. Less than 1 bug left after diluting! Consider changing your dilution factor.", stdout())
        stop("Cancelled the simulations. Less than 1 bug left after diluting! Consider changing your dilution factor.")
      }
      ## If not CASE 1, dilute normally
    } else {
      ## CASE 2 -- throw a warning if not many bugs left
      if (sum(this_timestep) <= 3) {
        write("WARNING: 3 or less bugs were left after diluting! Consider changing your dilution factor.", stdout())  # TODO 3 es ridículo, pensar uno más grande
        write("Running simulation anyway...", stdout())
      }
      this_timestep <- table(names(this_timestep)[.Internal(sample(
        x = length(this_timestep),
        size = sum(this_timestep)*dilution,
        replace = TRUE,
        prob = this_timestep))]) 
      
      temp <- empty                                                                #0
      temp[names(this_timestep)] <- this_timestep                                  #0 # aquí se desordena, pero solo hasta la siguiente línea.
      this_timestep <- temp[order(names(temp))]                                    #0
    }
    
    # (1) hago que los bichos se dupliquen al azar mediante sampling hasta 
    # llegar a la cantidad inicial (while loop); (2) run the substitution
    ns <- names(this_timestep)
    this_timestep <- as.vector(this_timestep)
    this_timestep <- full_growth(this_timestep, abun_total, grow_step)
    names(this_timestep) <- ns
    dil <- dil + 1
  }
  if (max(this_timestep)/sum(this_timestep) >= fixation_at) {
    write(paste0(fixation_at*100, "% fixation of ",
                 names(this_timestep)[this_timestep!=0], " after ", dil, " dilutions."), stdout())
  }
  return(this_timestep)
}





# print("Old:")
# start_time <- Sys.time()
# a <- simulate_timeseries_old(PCG_counts,
#                     dilution=8*10**(-3), # format: 0.1 instead of 10(%)
#                     no_of_dil=1200,
#                     abun_total=NULL,
#                     grow_step=1)
# end_time <- Sys.time()
# print(end_time - start_time)


print("New:")
start_time <- Sys.time()
b <- simulate_timeseries_new(PCG_counts,
                     dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                     no_of_dil=1200,
                     abun_total=NULL,
                     grow_step=1)
end_time <- Sys.time()
print(end_time - start_time)


print("Newer:")
start_time <- Sys.time()
c <- simulate_timeseries_newer(PCG_counts,
                             dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                             no_of_dil=1200,
                             abun_total=NULL,
                             grow_step=1)
end_time <- Sys.time()
print(end_time - start_time)


# print("Old:")
# print(a)
# print("\nNew:")
# print(b)
# print("\nNewer:")
# print(c)





# PCG_abund_temp <- mclapply(X = 1:4,
#                            FUN = function(iter) {
#                              # 1) simulation
# 
#                              trajectory <- simulate_timeseries_new(PCG_counts,
#                                                                    dilution=8*10**(-3), # format: 0.1 instead of 10(%)
#                                                                    no_of_dil=1200,
#                                                                    abun_total=NULL,
#                                                                    grow_step=1)
#                              
#                              print(paste("Simulation", iter, "finished for", s, " PCG",  PCG_name[PCG]))
#                              
#                              return(trajectory)
#                              
#                            }, mc.cores = 4)







# PCG_abund_temp <- mclapply(X = 1:no_of_simulations,
#                            FUN = function(iter) {
#                              # 1) simulation
#                              abun_total <- round(total_counts * percs[PCG_name[PCG]])
#                              if (abun_total == 0) { # we can't simulate anything if it's 0
#                                start <- as.count(my_transpose(PCG_counts))
#                                start <- start[order(names(start))] # this order thing is to keep indices consistent
#                                trajectory <- matrix(0, ncol=length (start), nrow = no_of_dil+1)
#                                rownames(trajectory) <- 0:no_of_dil
#                                colnames(trajectory) <- names(start)
#                                trajectory["0",]=start
#                              } else {
#                                trajectory <- simulate_timeseries(PCG_counts,
#                                                                  dilution = dilution,
#                                                                  no_of_dil = no_of_dil,
#                                                                  grow_step = grow_step,
#                                                                  abun_total = round(
#                                                                    total_counts *
#                                                                      percs[PCG_name[PCG]]),
#                                                                  keep_all_timesteps = save_all)
#                              }
#                              
#                              print(paste("Simulation", iter, "finished for", s, " PCG",  PCG_name[PCG]))
#                              
#                              return(trajectory)
#                              
#                            }, mc.cores = cores)