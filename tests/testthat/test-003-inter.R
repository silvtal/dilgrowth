# ==============================================================================
# This test creates simulation data similarly to test-000, with growth over 10
# transfers. Also plots the results similarly to test-002-plot-growth.
# In this case, we'll include a table of species interactions. We also compare
# with a percentage growth_step.
# ==============================================================================
pdf("../testresults/test-003.pdf", width = 8.5)

library("dilgrowth")
library("tidyverse")
library("gridExtra")

# run simuls

# ==============================================================================
## Input data
### ABUNTABLE: abundance table, output from BacterialCore.py
ABUNTABLE="../testdata/table_glucosa.txt"
### PCGTABLE: table with information with each PCG, output from BacterialCore.py
PCGTABLE="../testdata/pcgdata.txt"
### INTERTABLE: table with information about inter-species interactions
INTERTABLE="../testdata/interactions.txt"
### SAMPLENAMES: list/subset of the samples (column headers of $ABUNTABLE) we want to run simulations for
SAMPLENAMES=c("sa1")

## Simulation parameters
FIXATION_THRESHOLD="1" # relative abundance an OTU needs to have to be considered "fixed"
DILUTION=0.1
GROWTH_STEP=1
GD_PERC=FALSE

# ==============================================================================
# read abundance table
exp <- read.csv(
  ABUNTABLE,
  sep = "\t",
  skip = 1,
  row.names = 1,
  check.names = FALSE
)
exp <- exp[colnames(exp) != "taxonomy"] # remove taxonomy column if present
colnames(exp) <- as.character(colnames(exp)) # in case sample names are numbers
counts <- exp[SAMPLENAMES]
abun_total <- sum(counts)

# dilute
diluted_counts <- .my_transpose(counts)
diluted_counts <- table(names(diluted_counts)[.Internal(sample(
  x = length(diluted_counts),
  size = sum(diluted_counts)*DILUTION,
  replace = TRUE,
  prob = diluted_counts))])

# parse PCG table if everything's OK
pcg_table <- read.csv(PCGTABLE, sep="\t")
pcg_table <- pcg_table[1:(nrow(pcg_table)-1),] # remove last row (general info, not core info)
pcg_table <- pcg_table[c("Core", "Average", "Leaves")]
pcg_table$Average <- as.numeric(pcg_table$Average)

# carrying capacities for each PCG
abun_others <- abun_total * (1 - sum(pcg_table$Average))
carrying_capacities <- rep(round(abun_others), nrow(counts))
names(carrying_capacities) <- rep("others", nrow(counts))
for (group in 1:nrow(pcg_table)) {
  leaves <- strsplit(pcg_table$Leaves[group], ";")[[1]]
  names(carrying_capacities)[rownames(counts) %in% leaves] <- pcg_table$Core[group]
  carrying_capacities[rownames(counts) %in% leaves]        <- (pcg_table$Average[group] * abun_total) %>% round
}

# parse interactions table
interactions = read.csv(file = INTERTABLE, row.names = 1, check.names = F)


# ==============================================================================
# PLOT 4 PLOTS
plot_list <- list()

for (i in 1:2) {

  if (i == 1) {
    # with INTERACTIONS
    # ----------------------------------------------------------------------------
    # parse interactions table
    interactions = read.csv(file = INTERTABLE, row.names = 1, check.names = F)
    title = "with interactions"
  } else {
    # without INTERACTIONS
    # ----------------------------------------------------------------------------
    interactions[interactions != 0] <- 0
    title = "NO interactions"
  }

  # ==============================================================================
  # GROUPS, NON-LOGISTIC

  this_timestep <- as.vector(diluted_counts)
  this_timestep <- as.numeric(this_timestep)

  all_timesteps <- diluted_counts
  new_carrying_capacities <- carrying_capacities[rownames(counts) %in% names(diluted_counts)]
  new_interactions <- interactions[rownames(counts) %in% names(diluted_counts),
                                   rownames(counts) %in% names(diluted_counts)]

  # Check if any group has no abundance (optional step)
  sum_by_group <- c()
  groups <- unique(names(carrying_capacities))
  for (group in groups) {
    sum_by_group <- c(sum_by_group, sum(.my_transpose(counts)[names(carrying_capacities)==group]))
  }
  zero_groups <-  groups[!(groups %in% unique(names(new_carrying_capacities)))]

  # If there are zero groups, redefine abun_total
  if (length(zero_groups) > 0) {
    for (zg in zero_groups) {
      abun_total <- abun_total - carrying_capacities[zg]
      message(paste0("WARNING: The following group has no abundance and was excluded: ", zg, ". Adjusting abun_total to ", abun_total, " accordingly."))
    }
  }

  while (sum(this_timestep) < abun_total) {
    step          <- check_step(this_timestep, abun_total, GROWTH_STEP, GD_PERC)
    this_timestep <- growth(x = this_timestep,
                            growth_step = step,
                            carrying_capacities = new_carrying_capacities,
                            interactions = as.matrix(new_interactions))
    all_timesteps <- rbind(all_timesteps, this_timestep)
  }

  # PLOT

  leaf_to_pcg <- names(new_carrying_capacities)
  names(leaf_to_pcg) <- names(diluted_counts)

  # reshape the data into a long format
  all_timesteps <- data.frame(all_timesteps, check.names = F)
  all_timesteps$time <- 1:nrow(all_timesteps)
  df_long <- tidyr::gather(all_timesteps, key = "variable", value = "value", -time)

  # add a column to indicate the PCG name
  df_long$PCG <- leaf_to_pcg[df_long$variable]
  df_long$PCG[is.na(df_long$PCG)] <- "others"

  # create plot
  plot_list[[i]] <- ggplot(df_long, aes(x = time, y = value, group = variable, color = PCG)) +
    geom_line() +
    labs(x = "Time", y = "value", title = title) +
    geom_text(data = subset(df_long, time == max(time)),
              aes(label = variable, color = PCG),
              nudge_x = 40, nudge_y = 40, size = 2)

  # plot_list[[i]] <- ggplot(df_long, aes(x = time, y = value, group = variable, color = PCG)) +
  #   geom_line() +
  #   labs(x = "Time", y = "value", title = title)
  #


  # ==============================================================================
  # GROUPS, NON-LOGISTIC, GROWTH_STEP 1%

  GROWTH_STEP <- 0.01
  GD_PERC <- TRUE

  this_timestep <- as.vector(diluted_counts)
  this_timestep <- as.numeric(this_timestep)

  all_timesteps <- diluted_counts
  new_carrying_capacities <- carrying_capacities[rownames(counts) %in% names(diluted_counts)]
  new_interactions <- interactions[rownames(counts) %in% names(diluted_counts),
                                   rownames(counts) %in% names(diluted_counts)]

  # Check if any group has no abundance (optional step)
  sum_by_group <- c()
  groups <- unique(names(carrying_capacities))
  for (group in groups) {
    sum_by_group <- c(sum_by_group, sum(.my_transpose(counts)[names(carrying_capacities)==group]))
  }
  zero_groups <-  groups[!(groups %in% unique(names(new_carrying_capacities)))]

  # If there are zero groups, redefine abun_total
  if (length(zero_groups) > 0) {
    for (zg in zero_groups) {
      abun_total <- abun_total - carrying_capacities[zg]
      message(paste0("WARNING: The following group has no abundance and was excluded: ", zg, ". Adjusting abun_total to ", abun_total, " accordingly."))
    }
  }

  while (sum(this_timestep) < abun_total) {
    step          <- check_step(this_timestep, abun_total, GROWTH_STEP, GD_PERC)
    this_timestep <- growth(x = this_timestep,
                            growth_step = step,
                            carrying_capacities = new_carrying_capacities,
                            interactions = as.matrix(new_interactions))
    all_timesteps <- rbind(all_timesteps, this_timestep)
  }

  # PLOT

  leaf_to_pcg <- names(new_carrying_capacities)
  names(leaf_to_pcg) <- names(diluted_counts)

  # reshape the data into a long format
  all_timesteps <- data.frame(all_timesteps, check.names = F)
  all_timesteps$time <- 1:nrow(all_timesteps)
  df_long <- tidyr::gather(all_timesteps, key = "variable", value = "value", -time)

  # add a column to indicate the PCG name
  df_long$PCG <- leaf_to_pcg[df_long$variable]
  df_long$PCG[is.na(df_long$PCG)] <- "others"

  # create plot
  plot_list[[i+2]] <- ggplot(df_long, aes(x = time, y = value, group = variable, color = PCG)) +
    geom_line() +
    labs(x = "Time", y = "value", title = paste(title, "and 1%growth")) +
    geom_text(data = subset(df_long, time == max(time)),
              aes(label = variable, color = PCG),
              nudge_x = 40, nudge_y = 40, size = 2)
}
grid.arrange(
  grobs = plot_list,
  ncol = 2
)

dev.off()
