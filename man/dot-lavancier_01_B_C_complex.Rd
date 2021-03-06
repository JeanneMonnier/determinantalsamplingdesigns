% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lavancier.R
\name{.lavancier_01_B_C_complex}
\alias{.lavancier_01_B_C_complex}
\title{.lavancier_01_B_C_complex}
\usage{
.lavancier_01_B_C_complex(v, B = TRUE, C = rep(FALSE, nrow(v)))
}
\arguments{
\item{v}{a matrix.}

\item{B}{logical indicating which model of vector should be returned.}

\item{C}{a logical vector of length \code{n} (the number of elements sampled and the number of columns of v), indicating for each draw if the probability weights for the next draw should be plotted.}
}
\value{
A sample from DSD(vv*)
}
\description{
This function is called by \code{.lavancier_mult_complex}.
}
