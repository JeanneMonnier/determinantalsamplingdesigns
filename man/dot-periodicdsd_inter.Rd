% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/periodicdsd.R
\name{.periodicdsd_inter}
\alias{.periodicdsd_inter}
\title{periodicdsd_inter}
\usage{
.periodicdsd_inter(n, r, N, i)
}
\arguments{
\item{n}{a positive number, the number of items to be drawn, such that \code{n <= N}.}

\item{r}{a positive number, such that \code{r} and \code{N} are two relatively prime integers}

\item{N}{a positive number, the number of items to choose from.}

\item{i}{a positive integer, between 0 and n. It is the increment of the recursion.}
}
\value{
The ith step of the construction of the matrix v.
}
\description{
This function is called by \code{periodicdsd}.
}
