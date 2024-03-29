% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_timeseries.R
\name{simulate_timeseries}
\alias{simulate_timeseries}
\title{simulate_timeseries}
\usage{
simulate_timeseries(
  counts_data,
  carrying_capacities = NULL,
  interactions = NULL,
  logistic = FALSE,
  dilution = 8 * 10^(-3),
  no_of_dil = 12,
  fixation_at = 1,
  abun_total = NULL,
  growth_step = 1,
  is_growth_step_a_perc = FALSE,
  keep_all_timesteps = FALSE,
  allow_group_extinctions = TRUE,
  force_continue = FALSE
)
}
\arguments{
\item{growth_step}{Number of individuals that grow each timestep (Default: 1)}

\item{allow_group_extinctions}{If TRUE, simulations will continue even if one
or more groups go extinct, and the function will try to reach fixation in all
groups. Only applicable when carrying_capacities is not NULL (when there are
multiple functional groups) Also, this being FALSE does NOT affect groups
that were not in the community from the start (if there are missing groups
from the start, there will be a warning).}

\item{force_continue}{If TRUE, continue when ALL bugs go extinct because of a
strong dilution. This takes one random bug in order to continue, even when it
was actually diluted out. Affects the entire community, not each individual
functional group.}
}
\value{

}
\description{
Main function. Simulates the changes in abundances for a given initial
abundance data. The model has a dilution and growth system. After each
dilution, one or more random organism is duplicated (unless growth is set to
be logistic, in which case growth is not limited). Depends on c++ functions.
}
