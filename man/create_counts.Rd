% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_counts.R
\name{create_counts}
\alias{create_counts}
\title{create_counts}
\usage{
create_counts(
  exp,
  map,
  m_inic,
  allow_empty_transfers = FALSE,
  select_transfers = NULL,
  orig = "ORIG",
  transfer = "T",
  sa = "SA"
)
}
\arguments{
\item{exp}{Abundance table}

\item{map}{Metadata about m_inic, replicates, sample names...}

\item{m_inic}{Vector of samples to pick (by sample name)}

\item{select_transfers}{Vector of transfers to pick from the abundance table}

\item{orig}{Original sample tag column name in map}

\item{transfer}{Transfer tag column name in map}

\item{sa}{Sample name tag column name in map}
}
\value{

}
\description{
This function parses and puts together in a list of data.frames all the
count data for time series abundance data. For each timestep / transfer, it
saves the mean value of all available replicates.
}
\details{
If there are no available replicates for a given transfer, it saves the data
for the previous transfer (in other words, assumes no abundance changes took
place). This assumption can be overridden with the option
"allow_empty_transfers", that will set "NA" instead.

If select_transfers is not NULL, is has to be a vector. If you want to pick
only transfers 0, 2 and 4, you must input c(0, 2, 4)
}
