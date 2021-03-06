% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Smoothing.R
\name{RRM}
\alias{RRM}
\title{Repeated running medians}
\usage{
RRM(x, y, l = 5, eps = 1e-06, h.seq = c(7, 5, 3, 3))
}
\arguments{
\item{x}{- Raw x tracked coordinates}

\item{y}{- Raw y tracked coordinates}

\item{l}{- Minimum number of frames to be considered as arrest}

\item{eps}{- Cutoff value for determining arrests (maximum distance allowed to travel during arrest)}
}
\value{
Table of specfying beginning and ends of arrests
}
\description{
Repeated running medians
}
