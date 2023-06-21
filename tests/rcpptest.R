library("parallel")
library("dilgrowth")

counts_data <- read.csv("test-counts.csv", row.names = 1)

# function without c++
test_simulate_timeseries_old <- function (counts_data,
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

  start <- untb::as.count(.my_transpose(counts_data))
  start <- start[order(names(start))] # this order thing is to keep indices consistent

  # Prepare empty matrix
  len <- length (start)       # number of species
  if (is.null(abun_total)) {  # establish default value.
    abun_total <- sum (start) # total community abundance to be reached before each dilution
  }

  # initialize data.frame with starting abundances
  if (keep_all_timesteps){
    trajectory <- matrix(NA, ncol=len, nrow = no_of_dil+1)
    rownames(trajectory) <- 0:no_of_dil
    colnames(trajectory) <- names(start)
    trajectory["0",]=start
  }

  # and the "zero counter"
  empty <- start
  empty[1:length(empty)] <- 0

  this_timestep <- start

  message(paste0("Original abundance: ", sum(this_timestep)))
  message(paste0("After dilution: ", dilution*sum(this_timestep)))

  for (i in 1:no_of_dil){
    # Dilute (the initial inoculum too)

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
      # some zeroes might have appeared
      temp <- empty                                                                #0
      temp[names(this_timestep)] <- this_timestep                                  #0 # order is altered until following line
      this_timestep <- temp[order(names(temp))]                                    #0
    }

    # (1) bugs replicate randomly until reaching the initial abundance; (2) run
    # the substitution
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
      } # add newborn

    }

    # once the growth's finished, the next dilution can happen. But if
    # keep_all_timesteps, we have to save the abundances first
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

test_simulate_timeseries_new <- function (counts_data,
                                     dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                     no_of_dil=12,
                                     fixation_at=1,
                                     abun_total=NULL,
                                     grow_step=1,
                                     keep_all_timesteps=FALSE,
                                     force_continue=FALSE) {

  start <- untb::as.count(.my_transpose(counts_data))
  start <- start[order(names(start))] # this order thing is to keep indices consistent

  # Prepare empty matrix
  len <- length (start)       # number of species
  if (is.null(abun_total)) {  # establish default value.
    abun_total <- sum (start) # total community abundance to be reached before each dilution
  }

  # and the "zero counter"
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
      temp[names(this_timestep)] <- this_timestep                                  #0 # order is altered until following line
      this_timestep <- temp[order(names(temp))]                                    #0
    }

    # (1) bugs replicate randomly until reaching the initial abundance; (2) run
    # the simulation
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

test_simulate_timeseries_newer <- function (counts_data,
                                 dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                 no_of_dil=12,
                                 fixation_at=1,
                                 abun_total=NULL,
                                 grow_step=1,
                                 keep_all_timesteps=FALSE,
                                 force_continue=FALSE) {

  start <- untb::as.count(.my_transpose(counts_data))
  start <- start[order(names(start))] # this order thing is to keep indices consistent

  # Prepare empty matrix
  len <- length (start)       # number of species
  if (is.null(abun_total)) {  # establish default value.
    abun_total <- sum (start) # total community abundance to be reached before each dilution
  }

  # and the "zero counter"
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
      temp[names(this_timestep)] <- this_timestep                                  #0 # order is altered until following line
      this_timestep <- temp[order(names(temp))]                                    #0
    }

    # (1) bugs replicate randomly until reaching the initial abundance; (2) run
    # the substitution
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


print("New:")
start_time <- Sys.time()
b <- test_simulate_timeseries_new(counts_data,
                     dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                     no_of_dil=1200,
                     abun_total=NULL,
                     grow_step=1)
end_time <- Sys.time()
print(end_time - start_time)


print("Newer:")
start_time <- Sys.time()
c <- test_simulate_timeseries_newer(counts_data,
                             dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                             no_of_dil=1200,
                             abun_total=NULL,
                             grow_step=1)
end_time <- Sys.time()
print(end_time - start_time)
