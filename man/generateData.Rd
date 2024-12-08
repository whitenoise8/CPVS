% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gen_data.R
\name{generateData}
\alias{generateData}
\title{Generate data based on the functions used in the simulation study in the paper}
\usage{
generateData(n, d, ncp, minLen, sigma = 1, fix_seed = 0)
}
\arguments{
\item{n}{length of the series}

\item{d}{number of series}

\item{ncp}{number of change points}

\item{minLen}{minimum number of observations between consecutive change points}

\item{sigma}{standard deviation of the errors: if scalar, the same value is applied for all segments}

\item{fix_seed}{seed can be fixed for reproducibility. If 0, no seed is fixed}
}
\value{
\code{y}: the simulated data (a matrix \code{n x d} if multivariate case)

\code{signal}: the true simulated signal (a matrix \code{n x d} if multivariate case)

\code{cp}: the change point locations
}
\description{
Generate data based on the functions used in the simulation study in the paper
}
