% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{growth_one_group}
\alias{growth_one_group}
\title{growth_one_group}
\usage{
growth_one_group(this_timestep, growth_step, interactions = NULL)
}
\arguments{
\item{this_timestep}{A numeric vector representing the current abundance of
each organism in the population.}

\item{growth_step}{An integer representing the growth step for organisms in
the community, 1 by default.}

\item{interactions}{An optional numeric matrix representing the interaction
between organisms in the population. This argument is set to R_NilValue by
default.}
}
\description{
This Rcpp function simulates the growth of a community of organisms over a
specified time step. It takes the current community abundances, a fixed
growth step, and an optional interactions matrix as inputs, and returns
an updated abundances vector that reflects the growth of the community. The
function calculates the probability of growth for each organism based on
their current abundance and any specified interactions between organisms. The
function then randomly selects organisms to grow based on these probabilities
and increases their abundance by growth_step (1 by default).
The function takes four arguments:
}
