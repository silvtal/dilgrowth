#' Transpose a data.frame
#'
#' @param df A data.frame
#'
#' @importFrom data.table transpose
#'
#' @return The same dataframe but transposed keeping column/row names
#' @keywords internal
#' @export
.my_transpose <- function(df){
  t_df <- transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  return(t_df)
}


#' Create empty data.frame
#' Source: https://stackoverflow.com/questions/10689055/create-an-empty-data-frame
#' @param num_rows
#' @param num_cols
#'
#' @return Empty data.frame object
#' @keywords internal
#' @export
.create_empty_table <- function(num_rows, num_cols) {
  frame <- data.frame(matrix(NA, nrow = num_rows, ncol = num_cols))
  return(frame)
}


#' Create factor from taxonomy data. Deprecated in v4. Entero == 1; Pseudo == 2; others == 3
#'
#' @param tax
#' @param f_
#'
#' @return Factor with levels Entero, Pseudo and others
#' @keywords internal
#' @export
.pcgcodemaker <- function(tax, f_=FALSE) {
  if (f_ == TRUE) {
    e = " f__Enterobacteriaceae"
    p = " f__Pseudomonadaceae"
  } else {
    e = "Enterobacteriaceae"
    p = "Pseudomonadaceae"
  }
  entero000 = as.numeric(tax["Family"] != e) * 2 + 1
  pseudo000 = as.numeric(tax["Family"] == p)
  pcgcode = factor(entero000 - pseudo000, ordered = TRUE)
  rm(entero000)
  rm(pseudo000)
  return(pcgcode)
}


#' Create factor from taxonomy data. Deprecated in v4. Bulkholderiales, Pseudomonadales, others, Xanthomonadales
#'
#' @param tax
#' @param f_
#'
#' @return Factor with levels 1-4 (Bulkholderiales, Pseudomonadales, others, Xanthomonadales)
#' @keywords internal
#' @export
.pcgcodemaker2 <- function(tax, together=TRUE) {

  '%!in%' <- function(x, y) !('%in%'(x, y))

  if (together == TRUE) {
    b = c(" o__Burkholderiales", " o__Xanthomonadales")
    p = " o__Pseudomonadales"
  } else {
    b = " o__Burkholderiales"
    x = " o__Xanthomonadales"
    p = " o__Pseudomonadales"
  }
  b000 = as.numeric(sapply(tax["Order"], FUN = function(x) {
      x %!in% b
    }
  )) * 2 + 1
  pseudo000 = as.numeric(tax["Order"] == p)
  pcgcode = b000 - pseudo000

  if (together){
    pcgcode <- factor(as.numeric(pcgcode,ordered = TRUE))
  } else {
    x000 <- as.numeric(tax["Order"] == x)
    pcgcode[x000 == 1] <- 4
    pcgcode <- factor(as.numeric(pcgcode, ordered = TRUE))
  }
  return(pcgcode)
}
