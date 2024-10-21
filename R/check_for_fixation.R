#' check_for_fixation
#'
#' Checks for fixation. If there's one group (no PCG table), fixation occurs
#' when the relative abundance of one of the elements/species on the abundance
#' vector reaches a given fixation threshold (value from 0 to 1). If there are
#' multiple groups (PCG table given), fixation is checked for each group
#' individually.
#' Total fixation occurs when ALL groups reach fixation, meaning there's one
#' species in every group that has reached the fixation threshold within that
#' group. Total fixation is not evaluated in this function
#'
#' @param this_timestep abundance vector
#' @param carrying_capacities information about the groups each species belongs to
#' @param fixation_at fixation threshold
#' @param extinction_equals_fixation if TRUE, groups with an abundance of 0 will be considered fixated
#'
#' @return
#' @export
#'
#' @examples
#' this_timestep=c(20, 20, 20, 10, 10, 80)
#' carrying_capacities=c(100, 100, 100, 100, 100, 100)
#' names(carrying_capacities)=c("g1", "g1", "g1", "g2", "g2", "g2")
#' fixation_at <- 0.5
#' check_for_fixation(this_timestep, carrying_capacities, fixation_at)
check_for_fixation <- function(this_timestep, carrying_capacities, fixation_at, extinction_equals_fixation=FALSE) {
  if (is.null(carrying_capacities)) {
    fixated <- max(this_timestep)/sum(this_timestep) >= fixation_at
  } else {
    df <- data.frame(groups=names(carrying_capacities), abundances=as.numeric(this_timestep))
    total_abundances <- aggregate(abundances ~ groups, df, sum)
    max_abundances <- aggregate(abundances ~ groups, df, max)
    # Check if any element in the groups vector surpasses fixation_at% of the within-group total abundance
    # Check for each group
    fixated <- (max_abundances$abundances / total_abundances$abundances) >= fixation_at
    names(fixated) <- max_abundances$groups
  }
  if (extinction_equals_fixation) {
    fixated[is.na(fixated)] <- T
  } else {
    fixated[is.na(fixated)] <- F
  }
  return(fixated)
}
