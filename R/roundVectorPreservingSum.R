#' roundVectorPreservingSum
#'
#' This function rounds the input vector and then calculates the difference in
#' the sum of the original vector and the rounded vector. To fix this difference,
#' it randomly selects elements from the rounded vector and adjusts their values
#' to preserve the sum of the original vector.
#'
#' If carrying_capacities is not null, this is done in a group by group fashion.
#'
#' This function is tailored for this specific application, and our abundance
#' vectors will always sum an integer. But in the case it was given a vector
#' which sum totals a decimal number (i. e. 20.2), the resulting vector will not
#' have the same total value but a rounded one (i. e. 20)]
#'
#' This is needed to avoid decimals in the final simulated communities when
#' simulating multiple functional groups.
#'
#' @param vector
#'
#' @return
#' @export
#'
#' @examples
roundVectorPreservingSum <- function(vector, carrying_capacities = NULL) {

  rounded_vector <- round(vector) # start by rounding

  if (is.null(carrying_capacities)) { # only one group
    groups <- "gr"
    all_group_names <- rep("gr", length(vector))
  } else {
    all_group_names <- names(carrying_capacities)
    groups <- unique(all_group_names)
  }


  for (group in groups) {
    # then we check the difference (for this group only)
    diff_sum <- sum(vector[all_group_names==group]) - sum(rounded_vector[all_group_names==group])

    # then adjust
    if (diff_sum > 0) { # we need to sum
      indices <- which((all_group_names==group))
      indices <- sample(indices, diff_sum, replace = T) # (won't do anything if < 1)
      rounded_vector[indices] <- rounded_vector[indices] + 1
    } else if (diff_sum < 0) { # we need to substract
      indices <- which((all_group_names==group) & vector>0) # keep in mind we cannot substract from 0
      indices <- sample(indices, -diff_sum, replace = T) # (won't do anything if < 1)
      rounded_vector[indices] <- rounded_vector[indices] - 1
    }
  }
  return(rounded_vector)
}
