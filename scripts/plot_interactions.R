# Short script to generate a simple plot showing how interactions work
library(dilgrowth)
library(ggplot2)
library(tidyr)

# Auxiliary function: generates interaction matrix
create_interaction_matrix <- function(n, value) {
  mat <- matrix(0, n, n)
  mat[1,2] <- value  # Interacción de A sobre B
  return(mat)
}

# Options
set.seed(123)
initial_abundance <- c(10, 10)  # A y B con abundancia inicial 10
names(initial_abundance) <- c("A", "B")

# 1. No interactions
sim_no_interaction <- simulate_timeseries(
  counts_data = initial_abundance,
  no_of_dil = 10,
  keep_all_timesteps = TRUE
)

# 2. Negative interactions
sim_negative <- simulate_timeseries(
  counts_data = initial_abundance,
  interactions = create_interaction_matrix(2, -1),
  no_of_dil = 10,
  keep_all_timesteps = TRUE
)

# 3. Positive interactions
sim_positive <- simulate_timeseries(
  counts_data = initial_abundance,
  interactions = create_interaction_matrix(2, 0.5),
  no_of_dil = 10,
  keep_all_timesteps = TRUE
)

# Preparar datos para ggplot
prepare_data <- function(sim_data, title) {
  data.frame(
    time = 0:10,
    A = sim_data$A,
    B = sim_data$B
  ) %>%
    pivot_longer(
      cols = c(A, B),
      names_to = "Species",
      values_to = "Abundance"
    ) %>%
    mutate(Scenario = title)
}

# Combinar todos los datos
all_data <- rbind(
  prepare_data(sim_no_interaction, "No interactions"),
  prepare_data(sim_negative, "Negative interactions (A→B = -1)"),
  prepare_data(sim_positive, "Positive interactions (A→B = 0.5)")
)

# Crear las figuras
p <- ggplot(all_data, aes(x = time, y = Abundance, color = Species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~Scenario) +
  theme_minimal() +
  labs(
    x = "Number of dilutions",
    y = "Abundance",
    title = "Simulation of growth with different types of interactions"
  ) +
  scale_color_brewer(palette = "Set1") +
  theme(
    strip.text = element_text(size = 10),
    plot.title = element_text(hjust = 0.5)
  )

# Guardar la figura
ggsave("figures/interaction_simulations.png", p, width = 10, height = 4)
