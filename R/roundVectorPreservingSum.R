#' roundVectorPreservingSum
#'
#' This function rounds the input vector and then calculates the difference in
#' the sum of the original vector and the rounded vector. To fix this difference,
#' it randomly selects elements from the rounded vector and adjusts their values
#' to preserve the sum of the original vector.
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
roundVectorPreservingSum <- function(vector) {
  rounded_vector <- round(vector)  # Round the vector
  diff_sum <- sum(vector) - sum(rounded_vector)  # Calculate the difference in sum

  # Adjust the rounded vector to preserve the sum
  if (diff_sum > 0) {
    indices <- sample(length(vector), diff_sum)
    rounded_vector[indices] <- rounded_vector[indices] + 1
  } else if (diff_sum < 0) {
    indices <- sample(length(vector), -diff_sum)
    rounded_vector[indices] <- rounded_vector[indices] - 1
  }

  return(rounded_vector)
}
