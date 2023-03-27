# ==============================================================================
# This test plots the results from test-000, showing growth over 10 transfers.
# Only shows the abundances before each dilution - that is, at the end of each
# growth cycle.
# "old" plot --> non-logistic growth, no separation by group
# "new" plot --> logistic growth, each group has a different carrying capacity
#                as specified in the PCG table
# ==============================================================================

# load required packages
library(dplyr)
library(ggplot2)
library(tidyr)
# setwd("~/repos/dilgrowth/tests")

# function for reading files
read_all_files <- function(filename, transfers) {
  # create empty data frame to store results
  df <- data.frame()

  for (i in 0:transfers) {
    temp_df <- read.csv(paste0(filename, i, ".csv"), header = TRUE, check.names = F, row.names = 1)

    # add a column to indicate the file number
    temp_df$File <- i

    # add the data to the main data frame
    df <- bind_rows(df, temp_df)
  }
  return(df)
}

for (func in c('new', 'old')) {
  for (sample in c("sa1", "sa2")) {
    # loop through csv files and add them to the data frame
    df <- read_all_files(
      paste0("testresults/simuls_test_", sample, "/", func, "/simul_result_X", sample, "_t_"),
      transfers = 10
      )
    # remove transfers with NA values
    df <- remove_missing(df)

    # reshape the data into a long format
    df_long <- pivot_longer(df, -c(File), names_to = "Column", values_to = "Value")

    # read PCG table
    pcg_table <- read.csv("testdata/pcgdata.txt", sep="\t")
    pcg_table <- pcg_table[1:(nrow(pcg_table)-1),] # remove last row (general info, not core info)
    pcg_table <- pcg_table[c("Core", "Average", "Leaves")]

    leaf_to_pcg <- c()
    for (group in 1:nrow(pcg_table)) {
      leaves <- strsplit(pcg_table$Leaves[group], ";")[[1]]
      temp <- rep(pcg_table$Core[group], length(leaves))
      names(temp) <- leaves
      leaf_to_pcg <- c(leaf_to_pcg, temp)
    }

    # add a column to indicate the PCG name
    df_long$PCG <- leaf_to_pcg[df_long$Column]
    df_long$PCG[is.na(df_long$PCG)] <- "others"

    # create plot
    p <- ggplot(df_long, aes(x = File, y = Value, color = PCG, group = Column)) +
      geom_line() +
      labs(x = "File Number", y = "Value", title = func)
    print(p)

    # show fixed OTUs
    print(sort(df_long$Value[df_long$File==max(df_long$File)], decreasing = TRUE)[1:3])
  }
}
