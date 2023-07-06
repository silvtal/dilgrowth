# ==============================================================================
# This test creates simulation plots to show drift by showing the estochasticity
# of multiple instances of the same simulation. Some simulations will have
# interactions and others won't. Some will have a fixed growth_step value (1) and
# others will have a proggressive one (1%).
# ==============================================================================

library("dilgrowth")
library("tidyverse")
library("patchwork")
library("grid")
library("gridExtra")

# ==============================================================================
## Simulation parameters
FIXATION_THRESHOLD="1" # relative abundance an OTU needs to have to be considered "fixed"
DILUTION=0.15
DIL_NUMBER=10

# ==============================================================================
# set initial abundance vector and carrying capacities
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

# total final community size
abun_total <- sum(my_sample)

# set interactions
interactions <- data.frame(matrix(0,
                                  nrow = 7,
                                  ncol = 7),
                           row.names = names(my_sample))
colnames(interactions) <- names(my_sample)
# - [intergroup] otu1 grows worse in the presence of otu4
# - [intergroup] otu4 grows more in the presence of otu1
# - [intragroup] otu6 and otu7 help each other
interactions["otu1", "otu4"] <- -.5
interactions["otu4", "otu1"] <- +.5
interactions["otu6", "otu7"] <- +1
interactions["otu6", "otu7"] <- +1



# ==============================================================================
# PLOT 4 GROUPS OF PLOTS
grid_list <- list()
colors1 <- c("coral", "lightgreen")
colors2 <- c("darkred", "darkgreen")

# Loop over all combinations of growth_step and interaction
for (GROWTH_STEP in c(1, 0.01)) {
  if (GROWTH_STEP < 1) {GD_PERC=TRUE} else {GD_PERC=FALSE}

  for (include_interactions in c(FALSE, TRUE)) {
    title <- paste0("Grow_step: ", GROWTH_STEP, "; Interactions: ", ifelse(include_interactions, "YES", "NO"))

    plot_list1 <- list()
    plot_list2 <- list()

    # Create the grid and add it to the list

    ## -----------------------------------------------------------
    ## First, three plots for simulated growth within one dilution
    ## (same diluted starting community!)
    ## -----------------------------------------------------------
    # dilute
    diluted_counts <- table(names(my_sample)[.Internal(sample(
      x = length(my_sample),
      size = sum(my_sample)*DILUTION,
      replace = TRUE,
      prob = my_sample))])

    # Check for possible extinctions after dilution
    new_carrying_capacities <- carrying_capacities[colnames(my_sample) %in% names(diluted_counts)]
    new_interactions <- interactions[colnames(my_sample) %in% names(diluted_counts),
                                     colnames(my_sample) %in% names(diluted_counts)]
    if (!include_interactions) {
      new_interactions[new_interactions != 0] <- 0
    }

    # Check if any group has no abundance
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

    for (i in 1:3) {

      this_timestep <- as.vector(diluted_counts)
      this_timestep <- as.numeric(this_timestep)
      all_timesteps <- diluted_counts

      while (sum(this_timestep)+.5 < abun_total) {
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
      df_long_time <- tidyr::gather(all_timesteps, key = "variable", value = "value", -time)

      # add a column to indicate the PCG name
      df_long_time$PCG <- leaf_to_pcg[paste0(df_long_time$variable)]

      # create plot
      plot_list1[[i]] <- ggplot(df_long_time, aes(x = time, y = value, group = variable, color = PCG)) +
        geom_line() +
        scale_x_continuous(limits = c(min(df_long_time$time), max(df_long_time$time))) +
        labs(x = "Time", y = "value") +
        geom_text(data = subset(df_long_time, time == max(time)),
                  aes(label = variable, color = PCG),
                  nudge_x = 40, nudge_y = 40, size = 2) +
        scale_color_manual(values = colors1) +
        xlab("Time step") +
        ylab("Abundance") +
        theme(legend.position = "bottom")
      }




    ## ---------------------------------------------------------------
    ## Second, three plots for simulated growth over several dilutions
    ## ---------------------------------------------------------------
    for (i in 1:3) {
      diluted_counts <- my_sample
      all_dilutions  <- my_sample
      for (dil in 1:DIL_NUMBER) {
        # dilute
        diluted_counts <- table(names(diluted_counts)[.Internal(sample(
          x = length(my_sample),
          size = sum(my_sample)*DILUTION,
          replace = TRUE,
          prob = my_sample))])

        this_timestep <- as.vector(diluted_counts)
        this_timestep <- as.numeric(this_timestep)

        # Check for possible extinctions after dilution
        new_carrying_capacities <- carrying_capacities[colnames(my_sample) %in% names(diluted_counts)]
        new_interactions <- interactions[colnames(my_sample) %in% names(diluted_counts),
                                         colnames(my_sample) %in% names(diluted_counts)]
        if (!include_interactions) {
          new_interactions[new_interactions != 0] <- 0
        }

        # Check if any group has no abundance
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

        while (sum(this_timestep)+.5 < abun_total) {
          step          <- check_step(this_timestep, abun_total, GROWTH_STEP, GD_PERC)

          this_timestep <- growth(x = this_timestep,
                                  growth_step = step,
                                  carrying_capacities = new_carrying_capacities,
                                  interactions = as.matrix(new_interactions))
        }
        all_dilutions <- rbind(all_dilutions, this_timestep)

        # PLOT

        leaf_to_pcg <- names(new_carrying_capacities)
        names(leaf_to_pcg) <- names(diluted_counts)

        # reshape the data into a long format
        all_dilutions <- data.frame(all_dilutions, check.names = F)
        all_dilutions$time <- 1:nrow(all_dilutions)
        df_long <- tidyr::gather(all_dilutions, key = "variable", value = "value", -time)

        # add a column to indicate the PCG name
        df_long$PCG <- leaf_to_pcg[df_long$variable]

        # create plot
        plot_list2[[i]] <- ggplot(df_long, aes(x = time, y = value, group = variable, color = PCG)) +
          geom_line() +
          scale_x_continuous(limits = c(min(df_long$time), max(df_long$time))) +
          labs(x = "Time", y = "value") +
          geom_text(data = subset(df_long, time == max(time)),
                    aes(label = variable, color = PCG),
                    nudge_x = 40, nudge_y = 40, size = 2) +
          scale_color_manual(values = colors2) +
          xlab("Dilution number") +
          ylab("Abundance") +
          theme(legend.position = "bottom")

      }
    }

    print(title)

    plot_time <- wrap_plots(
      grobs = plot_list1,
      ncol = 1,
      guides = "collect"
    ) + plot_annotation(theme = theme(legend.position = "bottom"))

    plot_dilutions <- wrap_plots(
      grobs = plot_list2,
      ncol = 1,
      guides = "collect"
    ) + plot_annotation(theme = theme(legend.position = "bottom"))

    design_panel <- "
  11111111111111
  22222223333333
  22222223333333
  22222223333333
  22222223333333
  22222223333333
  22222223333333
  22222223333333
  22222223333333
  22222223333333
  22222223333333
"
    grid_list[[paste0(GROWTH_STEP, include_interactions)]] <- wrap_plots(
      grobs = list(textGrob(title, x=.3, y=.5,), plot_time, plot_dilutions),
      ncol = 2) + plot_layout(design = design_panel)
  }
}
# save plot
all_grids <- wrap_plots(
  grobs = grid_list,
  ncol = 2,
  convention = "keep"
)

pdf("../testresults/test-004.pdf", height = 22, width = 12)
all_grids
dev.off()
