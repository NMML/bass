\name{GetBASSExtent}
\alias{GetBASSExtent}
\title{Return an \code{extent} object from the \code{raster} package
cooresponding the BASS spatial grid.}
\usage{
  GetBASSExtent()
}
\value{
  extent object
}
\description{
  This function returns an object of type \code{extent}
  that specifies the coordinate boundaries of the BASS
  spatial grid. The following values describe the extent
  boundary within the BASS projection. x_min <- 0 x_max <-
  1.4e+06 y_min <- -3.8e+06 y_max <- -1.6e+06
}
\author{
  Josh London \email{josh.london@noaa.gov}
}

