#' Main function. Simulates the changes in abundances for a given initial
#' abundance data. The model has a dilution and growth system. After each
#' dilution, a random organism is duplicated
#'
#' Uses c++ functions for speed
#'
#' @param counts_data
#' @param dilution
#' @param no_of_dil
#' @param fixation_at
#' @param abun_total Total abundance for the community to grow in each dilution-growth
#' cycle. When this value is reached, the cycle ends and another dilution (if needed)
#' happens. If NULL, \code{abun_total} is set to the initial total abundance
#' @param grow_step Number of individuals that grow each timestep (Default: 1)
#' @param keep_all_timesteps
#' @param force_continue
#'
#' @importFrom untb as.count
#'
#' @return
#' @export
simulate_timeseries <- function (counts_data,
                                 carrying_capacities=NULL,
                                 dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                 no_of_dil=12,
                                 fixation_at=1,
                                 abun_total=NULL,
                                 grow_step=1,
                                 keep_all_timesteps=FALSE,
                                 force_continue=FALSE) {
  # counts_data puede ser un string o un dataframe de una columna como en el wrapper de /simuls
  if (is.atomic(counts_data)) {
    start <- counts_data
    if (!is.null(names(start))) {
      names(start) <- paste0("sp", 1:length(start))
    }
  } else {
    start <- as.count(.my_transpose(counts_data))
  }
  start <- start[order(names(start))] # this order thing is to keep indices consistent

  # Preparo la matriz que iré rellenando
  len <- length (start)    # número de especies

  if (is.null(abun_total)) { # establish default value.
    abun_total <- sum (start) # número de individuos al que dejaremos crecer antes de volver a diluir cada vez
  }

  # if keep_all_timesteps, initialize data.frame with starting abundances
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

      # We do this trick to keep all species "entries", even if they're 0 after
      # diluting. It's important for being able to put all timesteps together
      # intro a trajectory variable/file (keep_all_timesteps==TRUE) and for
      # saving multiple simulations together into a table
      temp <- empty
      temp[names(this_timestep)] <- this_timestep
      this_timestep <- temp[order(names(temp))]
    }

    # (1) hago que los bichos se dupliquen al azar mediante sampling hasta
    # llegar a la cantidad inicial (while loop)
    ns <- names(this_timestep)
    if (is.null(carrying_capacities)) {
      while (sum(this_timestep) < abun_total) {
        this_timestep <- growth(this_timestep,
                                abun_total,
                                grow_step)
      }
    } else {
        while (round(sum(this_timestep))+1 < abun_total) { # "round" to avoid infinitesimally small differences
         this_timestep <- growth_log(x = this_timestep,
                                     carrying_capacities = carrying_capacities)
        }
    }
    names(this_timestep) <- ns
    dil <- dil + 1

    if (keep_all_timesteps){
      # una vez crecidos, se puede diluir de nuevo: siguiente iteración del bucle
      # pero si keep_all(...), antes "secuenciamos" (guardamos las abundancias)
      trajectory[as.character(dil),]=this_timestep[colnames(trajectory)] # order by "trajectory" column names just in case
    }
  }
  if (max(this_timestep)/sum(this_timestep) >= fixation_at) {
    write(paste0(fixation_at*100, "% fixation of ",
                 names(this_timestep)[this_timestep!=0], " after ", dil, " dilutions."), stdout())
  }
  if (keep_all_timesteps){
    return(trajectory)
  } else {
    return(this_timestep)
  }}
