# ==============================================================================
# This test creates simulation data similarly to test-000, with growth over 10
# transfers. Also plots the results similarly to test-002-plot-growth.
# In this case, we'll include a table of species interactions. We also compare
# with a percentage growth_step.
# ==============================================================================
library("dilgrowth")
library("tidyverse")
library("grid")
library("gridExtra")
library("patchwork")

# run simuls

# ==============================================================================
## Input data (custom)
my_sample <- data.frame('otu1' = 3000,
                        'otu2' = 3000,
                        'otu3' = 1000,

                        'otu4' = 4000,
                        'otu5' = 4000,
                        'otu6' = 2000,
                        'otu7' = 500 )/10
carrying_capacities <- c('gr_1' = 7000,
                         'gr_1' = 7000,
                         'gr_1' = 7000,

                         'gr_2' = 10500,
                         'gr_2' = 10500,
                         'gr_2' = 10500,
                         'gr_2' = 10500)/10
my_colors <- c("otu1" = "#009900",
               "otu2" = "#466D1D",
               "otu3" = "#66FF66",

               "otu4" = "#B80F0A",
               "otu5" = "#FF0000",
               "otu6" = "#FF6666",
               "otu7" = "darkred")
# total final community size
abun_total <- sum(my_sample)
# set interactions
no_interactions <- data.frame(matrix(0,
                                  nrow = 7,
                                  ncol = 7),
                           row.names = names(my_sample))
colnames(no_interactions) <- names(my_sample)
my_interactions <- no_interactions
# - [intergroup] otu1 grows worse in the presence of otu4
# - [intergroup] otu4 grows more in the presence of otu1
# - [intragroup] otu6 and otu7 help each other
my_interactions["otu1", "otu4"] <- -1
my_interactions["otu4", "otu1"] <- +1
my_interactions["otu6", "otu7"] <- +1
my_interactions["otu7", "otu6"] <- +1

## Simulation parameters
FIXATION_THRESHOLD="1" # relative abundance an OTU needs to have to be considered "fixed"
DILUTION=0.1

abun_total <- sum(my_sample)

# dilute
diluted_counts <- table(names(my_sample)[.Internal(sample(
  x = length(my_sample),
  size = sum(my_sample)*DILUTION,
  replace = TRUE,
  prob = my_sample))])

# ==============================================================================
# PLOT 4 PLOTS
plot_list <- list()

for (i in 1:2) {

  if (i == 2) {
    # with INTERACTIONS
    # ----------------------------------------------------------------------------
    # parse interactions table
    interactions <- my_interactions
    title <- "with interactions"
  } else {
    # without INTERACTIONS
    # ----------------------------------------------------------------------------
    interactions <- no_interactions
    title = "NO interactions"
  }

  # ==============================================================================
  # GROUPS, NON-LOGISTIC

  GROWTH_STEP=1
  GD_PERC=FALSE

  this_timestep <- as.vector(diluted_counts)
  this_timestep <- as.numeric(this_timestep)

  all_timesteps <- diluted_counts
  new_carrying_capacities <- carrying_capacities[names(my_sample) %in% names(diluted_counts)]
  new_interactions <- interactions[names(my_sample) %in% names(diluted_counts),
                                   names(my_sample) %in% names(diluted_counts)]

  # Check if any group has no abundance (optional step)
  sum_by_group <- c()
  groups <- unique(names(carrying_capacities))
  for (group in groups) {
    sum_by_group <- c(sum_by_group, sum(my_sample[names(carrying_capacities)==group]))
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

  # reshape the data into a long format
  all_timesteps <- data.frame(all_timesteps, check.names = F)
  all_timesteps$time <- 1:nrow(all_timesteps)
  df_long <- tidyr::gather(all_timesteps, key = "variable", value = "value", -time)

  # add a column to indicate the PCG name
  leaf_to_pcg <- names(carrying_capacities)
  names(leaf_to_pcg) <- names(diluted_counts)

  df_long$PCG <- leaf_to_pcg[df_long$variable]

  # create plot
  plot_list[[i]] <- ggplot(df_long, aes(x = time, y = value, group = variable, linetype = PCG, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "value", title = title) +
    geom_text(data = subset(df_long, time == max(time)),
    aes(label = variable, color = variable),
    nudge_x = 40, nudge_y = 40, size = 3) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    scale_color_manual(values = my_colors)

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
  new_carrying_capacities <- carrying_capacities[names(my_sample) %in% names(diluted_counts)]
  new_interactions <- interactions[names(my_sample) %in% names(diluted_counts),
                                   names(my_sample) %in% names(diluted_counts)]

  # Check if any group has no abundance (optional step)
  sum_by_group <- c()
  groups <- unique(names(carrying_capacities))
  for (group in groups) {
    sum_by_group <- c(sum_by_group, sum(my_sample[names(carrying_capacities)==group]))
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

  # create plot
  plot_list[[i+2]] <- ggplot(df_long, aes(x = time, y = value, group = variable, linetype = PCG, color = variable)) +
    geom_line() +
    labs(x = "Time", y = "value", title = paste(title, "and 1%growth")) +
    geom_text(data = subset(df_long, time == max(time)),
              aes(label = variable, color = variable),
              nudge_x = 40, nudge_y = 40, size = 3) +
    scale_linetype_manual(values = c("solid", "dotted")) +
    scale_color_manual(values = my_colors)
}

design_panel <- "
  11111112222222
  11111112222222
  11111112222222
  33333334444444
  33333334444444
  33333334444444
  55555555555555
"
pdf("../testresults/test-005.pdf", width = 8.5)

plot_list[[5]] <- textGrob("Interactions are like follows:
- [intergroup] otu1 grows worse in the presence of otu4
- [intergroup] otu4 grows more in the presence of otu1
- [intragroup] otu6 and otu7 help each other", hjust = 0, gp =  gpar(fontsize = 8))

print(
  wrap_plots(
    grobs = plot_list,
    ncol = 2
  ) + plot_layout(design = design_panel)
)

dev.off()
