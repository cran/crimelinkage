% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/compareCrimes.R
\name{compareSpatial}
\alias{compareSpatial}
\title{Make spatial evidence variables}
\usage{
compareSpatial(C1, C2, longlat = FALSE)
}
\arguments{
\item{C1}{(n x 2) matrix of coordinates for the crimes}

\item{C2}{(n x 2) matrix of coordinates for the crimes}

\item{longlat}{(logical) if true, the the coordinates are in (Long,Lat), else}
}
\value{
numeric vector of distances between the crimes (in km)
}
\description{
Calculates spatial distance between crimes (in km)
}
\keyword{internal}

