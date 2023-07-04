#' simulate_timeseries
#'
#' Main function. Simulates the changes in abundances for a given initial
#' abundance data. The model has a dilution and growth system. After each
#' dilution, a random organism is duplicated. Uses c++ functions for speed.
#'
#' @param counts_data
#' @param carrying_capacities
#' @param interactions
#' @param logistic
#' @param dilution
#' @param no_of_dil
#' @param fixation_at
#' @param grow_step Number of individuals that grow each timestep (Default: 1)
#' @param keep_all_timesteps
#' @param allow_group_extinctions If TRUE, simulations will continue even if one
#' or more groups go extinct, and the function will try to reach fixation in all
#' groups. Only applicable when carrying_capacities is not NULL (when there are
#' multiple functional groups) Also, this being FALSE does NOT affect groups
#' that were not in the community from the start (if there are missing groups
#' from the start, there will be a warning).
#' @param force_continue If TRUE, continue when ALL bugs go extinct because of a
#' strong dilution. This takes one random bug in order to continue, even when it
#' was actually diluted out. Affects the entire community, not each individual
#' functional group.
#'
#' @return
#' @export
simulate_timeseries <- function (counts_data,
                                 carrying_capacities=NULL,
                                 interactions=NULL,
                                 logistic=FALSE,
                                 dilution=8*10**(-3), # format: 0.1 instead of 10(%)
                                 no_of_dil=12,
                                 fixation_at=1,
                                 abun_total=NULL,
                                 grow_step=1,
                                 is_grow_step_a_perc=FALSE,
                                 keep_all_timesteps=FALSE,
                                 allow_group_extinctions=FALSE,
                                 force_continue=FALSE) {
  already_excluded <- c() # we'll exclude groups with total abundance of 0
  if (!is.null(already_excluded)) {warning(paste("Group(s) absent from starting data:", already_excluded))}

  # counts_data puede ser un string o un dataframe de una columna como en el wrapper de /simuls
  if (is.atomic(counts_data)) {
    start <- counts_data
    if (!is.null(names(start))) {
      names(start) <- paste0("sp", 1:length(start))
    }
  } else {
    start <- .my_transpose(counts_data)
  }

  # Prepare empty matrix
  len <- length (start)       # number of species
  if (is.null(abun_total)) {  # establish default value.
    abun_total <- sum (start) # total community abundance to be reached before each dilution
  }
  ns <-names(start)
  # if keep_all_timesteps, initialize data.frame with starting abundances
  if (keep_all_timesteps){
    trajectory <- data.frame(matrix(NA, ncol=len, nrow = no_of_dil+1))
    rownames(trajectory) <- 0:no_of_dil
    colnames(trajectory) <- ns
    trajectory["0",]=start
  }

  # this empty vector allows us to keep all species names even when not present
  # after a dilution
  empty <- start
  empty[1:length(empty)] <- 0

  this_timestep <- start
  dil <- 0

  # Check if fixation has occured already (if it has, simulations stop)
  # When checking for fixation, if there are multiple groups we must check for
  # fixation in each group separately
  not_fixated <- !check_for_fixation(this_timestep, carrying_capacities, fixation_at)

  while ((dil < no_of_dil) & not_fixated) {
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
      # into a trajectory variable/file (keep_all_timesteps==TRUE) and for
      # saving multiple simulations together into a table
      temp <- empty
      temp[names(this_timestep)] <- this_timestep
      this_timestep <- temp[names(empty)]
    }

    # Prepare the growth loop
    if (!(is.null(interactions))) {
      interactions <- as.matrix(interactions)
    }

    this_timestep <- as.vector(this_timestep)
    this_timestep <- as.numeric(this_timestep)
    if (is.null(carrying_capacities)) {
      # ==========================================
      # Growth without groups (growth_one_group())
      # ==========================================
      while (sum(this_timestep) < abun_total) {
        step          <- check_step(this_timestep, abun_total, grow_step, is_grow_step_a_perc)
        this_timestep <- growth_one_group(this_timestep,
                                          step,
                                          interactions)
      }
    } else {
      # =========================================
      # Growth by group (growth() / growth_log())
      # =========================================
      # Check if any group has no abundance
      sum_by_group <- c()
      groups <- unique(names(carrying_capacities))
      for (group in groups) {
        this_timestep[names(carrying_capacities)==group]
        sum_by_group <- c(sum_by_group, sum(this_timestep[names(carrying_capacities)==group]))
      }
      zero_groups <-  groups[sum_by_group == 0]

      # If there are zero groups
      if (length(zero_groups) > 0) {
        if (!allow_group_extinctions) {# If not allowed, stop()
          stop(paste0("Some group(s) have gone extinct (", zero_groups, "). Stopping the simulations..."))
        } else { # If they are allowed, redefine abun_total
          for (zg in zero_groups) {
            if (!(zg %in% already_excluded)) {
              abun_total <- abun_total - carrying_capacities[zg][[1]]
              already_excluded <- c(already_excluded, zg)
              message(paste0("WARNING: The following group has no abundance and was excluded: ", zg, ". Adjusting abun_total to ", abun_total, " accordingly."))
            }
          }
        }
      }
      # Choose growth function; we could use "empty" groups for speed reasons;
      # that is, redefining carrying_capacity and interactions so they don't
      # contain species from zero groups. But it's not necessary.
      # The important bit is to redefine abun_total to avoid an infinite loop.
      # (As a side note: imagine groups G1(70%), G2(20%) and G3(10%). If G2
      # becomes a zero group after dilution, group_growth in growth_per_group()
      # would go from 70/100=0.7 to 70/80=0.875 for G1 and from 10/100=0.1 to
      # 10/80=0.125 for G3. However, the G1/G3 ratio stays the same: 0.7/0.1=
      # =0.875/0.125=7. So it doesn't matter if we don't filter out zero groups)
      if (logistic) {
        while (round(sum(this_timestep)) < abun_total) { # "round" to avoid infinitesimally small differences
          this_timestep <- growth_log(
            x = this_timestep,
            carrying_capacities = carrying_capacities,
            interactions = interactions
          )
        }
      } else {
        while (round(sum(this_timestep)) < abun_total) {
          step          <- check_step(this_timestep, abun_total, grow_step, is_grow_step_a_perc)
          this_timestep <- growth(
            x = this_timestep,
            carrying_capacities = carrying_capacities,
            grow_step = step,
            interactions = interactions
          )
        }
      }
    }

    names(this_timestep) <- ns

    dil <- dil + 1

    if (keep_all_timesteps){
      # once the growth's finished, the next dilution can happen. But if
      # keep_all_timesteps, we have to save the abundances first
      trajectory[as.character(dil),] <- this_timestep
    }

    # Check again for fixation before next iteration
    not_fixated <- !check_for_fixation(this_timestep, carrying_capacities, fixation_at)
  }

  # If the loop has stopped because of fixation
  if (not_fixated == FALSE) {
    write(paste0(fixation_at*100, "% fixation of ",
                 ns[this_timestep!=0], " after ", dil, " dilutions."), stdout())
  }

  if (keep_all_timesteps){
    return(trajectory)
  } else {
    return(this_timestep)
  }}
