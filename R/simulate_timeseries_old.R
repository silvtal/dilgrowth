#' OLD version for the main function. Does not use c++ functions. No
#' \code{fixation_at} parameter. Accepts \code{growth_rate}.
#'
#' @param counts_data
#' @param dilution
#' @param no_of_dil
#' @param abun_total Total abundance for the community to grow in each dilution-growth
#' cycle. When this value is reached, the cycle ends and another dilution (if needed)
#' happens. If NULL, \code{abun_total} is set to the initial total abundance
#' @param grow_step Number of individuals that grow each timestep (Default: 1)
#' @param growth_rate Indicates the name of an optional column that determines how much a bug
#' grows. Each time it grows, it won't grow by 1 but by this factor. It
#' indicates the grow rate in percent mass.
#' WARNING: It can't be a decimal number because of how the sampling works.
#' That means the count_data data.table, which includes the growth rates,
#' must be multiplied by a number that turns all rates into an integer. The
#' abun_total AND grow_step arguments must be multiplied by this number, too
#' @param keep_all_timesteps
#'
#' @importFrom untb as.count
#' @importFrom stats setNames
#'
#' @return
#' @export
simulate_timeseries_old <- function (counts_data,
                                     dilution = 8 * 10 ** (-3),
                                     # format: 0.1 instead of 10(%)
                                     no_of_dil = 12,
                                     abun_total = NULL,
                                     grow_step = 1,
                                     keep_all_timesteps = FALSE,
                                     growth_rate = NULL) {
  start <- as.count(.my_transpose(counts_data))
  start <- start[order(names(start))] # this order thing is to keep indices consistent

  # Preparo la matriz que iré rellenando
  len <- length (start)    # número de especies
  if (is.null(abun_total)) {
    # establish default value.
    abun_total <-
      sum(start) # número de individuos al que dejaremos crecer antes de volver a diluir cada vez
  }

  # initialize data.frame with starting abundances
  if (keep_all_timesteps) {
    trajectory <- matrix(NA, ncol = len, nrow = no_of_dil + 1)
    rownames(trajectory) <- 0:no_of_dil
    colnames(trajectory) <- names(start)
    trajectory["0", ] = start
  }

  # y el "zero counter"
  empty <- start
  empty[1:length(empty)] <- 0

  this_timestep <- start

  message(paste0("Original abundance: ", sum(this_timestep)))
  message(paste0("After dilution: ", dilution * sum(this_timestep)))

  for (i in 1:no_of_dil) {
    ## CASE 1 -- All taxa extinct or almost
    if (trunc(sum(this_timestep) * dilution) == 0) {
      # feb 13 2022: step to 1, and pick a random single bug anyway
      message("WARNING: less than 1 bug left after diluting! Consider changing your dilution factor.")
      message("Picking one random bug and running simulation anyway...")
      i <- c(1:length(this_timestep))[.Internal(sample(
        x = length(this_timestep),
        size = 1,
        replace = FALSE,
        prob = this_timestep
      ))]
      this_timestep    <- empty
      this_timestep[i] <- this_timestep[i] + 1
      print(this_timestep)
      ## If not CASE 1, dilute normally
    } else {
      ## CASE 2 -- throw a warning if not many bugs left
      if (sum(this_timestep) <= 3) {
        message("WARNING: 3 or less bugs were left after diluting! Consider changing your dilution factor.")
        message("Running simulation anyway...")
        ### Changed feb 1 2022: i had a "skip simuls for this case" but I changed
        ###                     it so it's "change step to one and throw WARNING"
      }
      this_timestep <- table(names(this_timestep)[.Internal(sample(
        x = length(this_timestep),
        size = sum(this_timestep) * dilution,
        replace = TRUE,
        prob = this_timestep
      ))])
      # indicamos los ceros que hayan aparecido
      temp <-empty                                                                #0
      temp[names(this_timestep)] <- this_timestep                                 #0 # aquí se desordena, pero solo hasta la siguiente línea.
      this_timestep <- temp[order(names(temp))]                                   #0
    }

    # (1) hago que los bichos se dupliquen al azar mediante sampling hasta
    # llegar a la cantidad inicial (while loop); (2) run the substitution
    while (sum(this_timestep) < abun_total) {
      # si aún debe seguir creciendo
      if (sum(this_timestep) < grow_step) {
        # the step is too big, reduce it
        step <- max(trunc(sum(this_timestep) / 2))
        # warning(paste("Step reduced to",step))
      } else {
        step <- grow_step
      }

      ## Avoid growing too much (when step>1)
      if ((sum(this_timestep) + step) > abun_total) {
        step = (abun_total - sum(this_timestep))
      }

      ## Once the step is decided, sample
      new_bugs <- c(1:length(this_timestep))[.Internal(sample(
        x = length(this_timestep),
        size = step,
        replace = TRUE,
        prob = this_timestep
      ))]
      for (i in new_bugs) {
        this_timestep[i] <- this_timestep[i] + 1
      } # añado al que ha nacido al censo

    }

    # una vez crecidos, se puede diluir de nuevo: se repite el bucle
    # pero primero "secuenciamos" (guardamos las abundancias actualizadas)
    if (keep_all_timesteps) {
      trajectory[as.character(i), ] = this_timestep[colnames(trajectory)] # order by "trajectory" column names just in case
    }
  }
  message(paste0("Reached abundance of: ", sum(this_timestep)))

  if (keep_all_timesteps) {
    return(trajectory)
  } else {
    return(this_timestep %>% setNames(names(start)))
  }
}
