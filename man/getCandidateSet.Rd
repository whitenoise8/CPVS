% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/candidate_cp.R
\name{getCandidateSet}
\alias{getCandidateSet}
\title{Get the candidate set of change points using Optimistic Search and gain function based on Marginal Likelihood/ELBO (see paper for details)}
\usage{
getCandidateSet(
  y,
  Xsin,
  Xcos,
  Ne,
  sigma2 = 0,
  nu = 1/2,
  thresh = 0,
  rate_c = 0.1,
  Trace = 0
)
}
\arguments{
\item{y}{vector of observations or matrix (\code{T x d}) if multivariate case}

\item{Xsin}{matrix (\code{T x p}) or array (\code{T x p x d}), if multivariate, of sine functions X^(1) based on the grid}

\item{Xcos}{matrix (\code{T x p}) or array (\code{T x p x d}), if multivariate, of cosine functions X^(2) based on the grid}

\item{Ne}{number of effects in the SuSiE algorithm}

\item{sigma2}{default is 0, meaning that the variance of the error term is estimated. It can be fixed to a specific value.}

\item{nu}{default is 1/2. \nu parameter in the optimistic search algorithm}

\item{thresh}{default is equal to 0. Threshold in the gain function to define a change point}

\item{rate_c}{defines the minimum distance from the boundaries for a change point rate_c*T}

\item{Trace}{if 1, print progress of the search}
}
\value{
\code{cp_list}: list of the candidate change point locations

\code{delbo}: list of the gain function values for each element in \code{cp_list}
}
\description{
Get the candidate set of change points using Optimistic Search and gain function based on Marginal Likelihood/ELBO (see paper for details)
}
