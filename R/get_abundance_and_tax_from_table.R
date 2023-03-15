#' This function reads a .csv/.tsv abundance table, where the last
#' column (or row) is the taxonomy
#'
#' @param exp_file .csv/.tsv abundance table; required
#' @param species_are_rows
#' @param sep
#' @param row.names
#' @param skip
#' @param taxa_fields
#' @param tax_col_name
#' @param tax_sep 4 chars or longer "breaks" renamer (will include this as if it were a real name)
#' @param NA_option
#' @param check.names
#'
#' @return two separate data.frames: abundances and taxa data
#' @importFrom tidyr separate
#' @importFrom stats setNames
#' @importFrom utils read.csv
#' @importFrom gsubfn list
#' @export
get_abundance_and_tax_from_table <- function (exp_file,
                                              species_are_rows=TRUE,
                                              sep="\t",
                                              row.names=1,
                                              skip=1,
                                              taxa_fields=NULL,
                                              tax_col_name="taxonomy",
                                              tax_sep=";",
                                              NA_option="___",
                                              check.names=FALSE) {
  if (is.null(taxa_fields)) {taxa_fields = c("Kingdom","Phylum","Class","Order","Family","Genus","Species")}

  exp <-
    read.csv(
      exp_file,
      sep = sep,
      skip = skip,
      row.names = row.names,
      check.names = check.names
    )

  if (!species_are_rows) {
    exp <- .my_transpose(exp)
  }

  tax <- exp["taxonomy"]
  exp <- exp[1:dim(exp)[2] - 1]
  tax <- tax %>% separate("taxonomy", sep = tax_sep, taxa_fields)
  tax[is.na(tax)] <- NA_option # avoid na-related errors

  return(list(exp, tax))
}
