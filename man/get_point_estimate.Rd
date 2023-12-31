% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_regimes.R
\name{get_point_estimate}
\alias{get_point_estimate}
\title{Create an estimate on posterior distribution of change-points.}
\usage{
get_point_estimate(regime_fit_object, prob_cutoff)
}
\arguments{
\item{regime_fit_object}{bayesWatch object. Fit with the bayesWatch method.}

\item{prob_cutoff}{float in (0,1). Posterior probabilities above this cutoff will be considered changepoints.}
}
\value{
vector. Indicator values corresponding to change-point locations.
}
\description{
Given a bayesWatch object and a probability cutoff, finds change-points.
}
